// event time-ordering - scrambles original event ordering
//    => midas event/bank headers no longer connected to their events
//         (but there is nothing useful there anyway)
// we used to use midas timestamps, but event timestamps are more useful
// ----------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "grif-replay.h"
#include "midas-format.h"

#define MAX_GRIFC 16

unsigned event_buffer[EVENTBUFSIZE+128];
volatile unsigned long eventbuf_rdpos;
volatile unsigned long eventbuf_wrpos;
int reorder_insert_event_check(int grifc);

volatile long tsevents_in;
int reorder_events_read;

// output ...
#define TS_EVENT_BUFFER 8192
volatile unsigned long output_ts;
pthread_mutex_t nxtlock;
long tsevents_out;

extern Sort_metrics diagnostics;
void reorder_status(int current_time)
{
   int sum, *err = diagnostics.reorder_error;
   sum = err[0] + err[2] + err[3] + err[4] + err[5] + err[6] + err[7];// err[8] is in addition to other
                                                                      // errors - do not add to total
   printf("Reorder: in:%10ld out:%10ld err:%10d[%5.1f%%] [Desync:%d]\n         ",
          tsevents_in, tsevents_out, sum, (100.0*sum)/tsevents_in, err[8] );
   printf("[init:%d format:%d length:%d early:%d late:%d,%d, unk:%d]\n",
          err[2], err[0], err[3], err[5], err[4], err[6], err[7] );
}

//////////////////////////////////////////////////////////////////////////////
//this type of sort is a bucket-sort with insertion-sort into each time-bucket
//  it is planned to optimize the insertion sort with one skip-list if needed
//////////////////////////////////////////////////////////////////////////////
#define TOO_EARLY_CUTOFF 3000000000
#define REORDER_EVENTS    65536  // *0.5 => max 32k events stored in buffers
#define REORDER_TSLOTS  1000000  // 1 million [~1.25 seconds]
#define BUCKET_SIZE_BITS      7  // 128 timestamps per slot -> 1us
#define OVERFULL_FRACTION   0.5
#define OUTPUT_FRACTION    0.25
#define INIT_WAIT           250 // allow # junk events per grifc at run start
#define REORDER_MAXEVENTSIZE 20 // max 20 words - 80 bytes
typedef struct reorderbuf_struct Tsbuf;
struct reorderbuf_struct {
   volatile Tsbuf *next; unsigned long ts;
   volatile int in_use; int evlen;
   int event[REORDER_MAXEVENTSIZE];
};
Tsbuf  reorder_buffer[REORDER_EVENTS]; // ~100bytes ea -> 100Mbytes
volatile Tsbuf *tslot[REORDER_TSLOTS];
int reorder_bufpos; // next slot to be written

#define MAX_GRIFC 16
int reorder_init[MAX_GRIFC];

pthread_mutex_t nxtlock;
void reorder_main(Sort_status *arg)
{
   int ev_done, ts_slot, len, startup, *err = diagnostics.reorder_error;
   int i,  grifc,  type,  rd_avail,  ts_stat,  err_format;
   unsigned int usecs=1000, *evptr, *evstart, *bufend;
   volatile Tsbuf *bufptr, *nxtptr, *newptr;
   unsigned long ts, tslo;

   startup = 1;
   memset(err,           0, REORDER_ERRTYPES*sizeof(int));
   //for(i=0; i<MAX_GRIFC; i++){ reorder_init[i] = INIT_WAIT; }
   printf("starting reorder input ...\n");
   for(i=0; i<REORDER_EVENTS; i++){
      reorder_buffer[i].in_use = 0;  reorder_buffer[i].next = NULL;
   }
   for(i=0; i<REORDER_TSLOTS; i++){ tslot[i] = NULL; }
   tsevents_in = reorder_bufpos = len = ev_done = ts_stat = err_format = 0;
   evstart = NULL; evptr = bankbuf; bufend = bankbuf + BANK_BUFSIZE;
   while(1){
      if( tsevents_in - tsevents_out > OVERFULL_FRACTION*REORDER_EVENTS ){
         usleep(usecs); continue; // do not over-fill buffer
      }
      if( (rd_avail = bankbuf_wrpos - bankbuf_rdpos) < 1 ){
         if( arg->end_of_data ){ break; }
         usleep(usecs); continue;
      }
      type = (((*evptr) >> 28) & 0xf);
      if( evstart == NULL && type != 0x8 ){ err_format=1; }
      switch(type){
      case 0x8:
         if( evstart != NULL ){ err_format=1; } else { evstart=evptr; }
         grifc = ((*evptr) >> 16) & 0xF;  break;
      case 0xE:
         if( ts_stat != 2 ){ err_format=1; } ev_done = 1; break;
      case 0xA:
         if( ts_stat != 0 ){ err_format=1; }
         tslo = *evptr & 0xFFFFFFF; ts_stat = 1; break;
      case 0xB:
         if( ts_stat != 1 ){ err_format=1; }
         ts   = *evptr &    0x3FFF;
         ts <<= 28; ts += tslo; ts_stat = 2; break;
      default:  break;
      }
      ++len;
      ++bankbuf_rdpos; ++evptr; if( evptr >= bufend ){ evptr -= BANK_BUFSIZE; }
      if( !ev_done ){ continue; }
      ++reorder_events_read;
      // now have full event
      if( err_format ){
         ++err[ERR_FORMAT_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      if( len > REORDER_MAXEVENTSIZE ){
         ++err[ERR_LENGTH_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      // often, at start of run, get junk data from prev run(large timestamps)
      // try and discard this (timestamps larger than 2.6 seconds)
      if( reorder_init[grifc] != 0 ){
         ++err[ERR_INIT_IN]; err[ERR_WORDS_IN] += len;
         if( (ts>>28) == 0 ){ reorder_init[grifc] = 0; }
         else             { --reorder_init[grifc]; }
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      // now have full event without any obvious format errors
      
      // events that are so late, that we have already output their timeslot
      // can no longer be made to be in order
      // - drop them for now (may include anyway - mark as just for singles?)
      if( ts <= output_ts ){
         ++err[ERR_LATE_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      // events that are Extremely early (firmware bugs/data corruption)
      // these would stop reuse of their and subsequent buffer slots until
      // their timeslot comes (reuse would be expected ~64k events later)
      //   i.e. buffer would become blocked for a long time
      //   discard for now, until implement a way of avoiding blocked buf
      if( ts > (output_ts + TOO_EARLY_CUTOFF) ){
         ++err[ERR_EARLY_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      bufptr = NULL;  ts_slot = (ts >> BUCKET_SIZE_BITS) % REORDER_TSLOTS;
      // check if next buffer slot is available to store this event
      // LOOP over buffer (starting at next)
      //   NOTE: looping, instead of just waiting for nxt to become free
      //          - avoids single early event blocking buffer (see above)
      i=0; do {
         ++i; newptr = reorder_buffer + reorder_bufpos++;
         if( reorder_bufpos >= REORDER_EVENTS ){ reorder_bufpos = 0; } // wrap
         if( i >= REORDER_EVENTS ){ // ALL BUSY
            i=0; usleep(100*usecs);
            printf("################### REORDER IN_WAIT\n");
         }
      } while( newptr->in_use );
      // copy event to buf slot found above
      newptr->in_use =  1; newptr->evlen =  len;
      newptr->ts     = ts; newptr->next  = NULL;
      if( evptr >= evstart ){
         memcpy((char *)(newptr->event),(char *)(evstart),4*len);
      } else {  // wrapped in middle of event
         i = bufend - evstart;
         memcpy((char *)(newptr->event),  (char *)(evstart),i);
         memcpy((char *)(newptr->event)+i,(char *)(bankbuf),4*len-i);
      }
      // update linked list ...
      bufptr = NULL;
      pthread_mutex_lock(&nxtlock);
      if( (nxtptr = tslot[ts_slot]) == NULL ){ // currently empty
         tslot[ts_slot] = newptr;
      } else {                              // insert into list IN TIME ORDER
         while( 1 ){ // at this point - bufptr is previous, nxtptr is current
            if( ts <= nxtptr->ts || nxtptr->in_use == 0 ){ // insert here
               if( bufptr == NULL ){
                  tslot[ts_slot] = newptr; newptr->next = nxtptr;
               } else {
                  bufptr->next = newptr; newptr->next=nxtptr;
               }
               break;
            }
            if( nxtptr->next == NULL ){ // insert at end of chain
               nxtptr->next = newptr; break;
            }
            bufptr = nxtptr; nxtptr = nxtptr->next;
         }
      }
      pthread_mutex_unlock(&nxtlock);  ++tsevents_in;
      evstart = NULL;  len = ev_done = ts_stat = 0; err_format=0;
   }
   arg->reorder_in_done = 1;
   printf("Reorder: end input thread\n");
   return;
}

#define MAX_DATA_GAP 4000000000 // 40 seconds
void reorder_out(Sort_status *arg)
{
   int i, j, wr_avail, wrpos, ts_slot, *err;
   int bucket_length = (1<<BUCKET_SIZE_BITS);
   volatile Tsbuf *buf, *nxt;
   unsigned long ts, prev_ts;
   unsigned int usecs=100;

   usleep(100*usecs);
   err = diagnostics.reorder_error;
   printf("starting reorder output ...\n");
   tsevents_out = prev_ts = output_ts = 0;
   eventbuf_rdpos = eventbuf_wrpos = 0;
   ts = 0 - bucket_length;
   while(1){
      i = tsevents_in - tsevents_out;
      if( i < OUTPUT_FRACTION * REORDER_EVENTS && !arg->reorder_in_done ){
         usleep(usecs); continue;
      }
      while(1){ ts += bucket_length; // check slots, in order, for next event
         if( ts - prev_ts >= MAX_DATA_GAP ){ // what is this?
            if( arg->reorder_in_done ){
               arg->reorder_out_done = 1;
               printf("Reorder: end output thread [ts:%3ds] ...\n",
                        (int)(ts/100000000) );
               reorder_status(0);
               return;
            } else { //  error?
               printf("REORDER ERROR-VERY-LONG-DATA-GAP\n");
            }
         }
         ts_slot = (ts >> BUCKET_SIZE_BITS) % REORDER_TSLOTS;
         if( (buf = tslot[ts_slot]) == NULL ){ continue; }
         if( buf->ts > ts+bucket_length ){ continue; }
         break;
      }
      prev_ts = ts;
      // have found event to output (ts is now - or earlier)
      //    iterate over this in-use linked list
      //    NOTE: as this list is sorted in timestamp order, this iteration 
      //    will only be over any equal timestamp events at start of list
      while(1){
         nxt = buf->next;
         if( buf->ts > ts+bucket_length ){ break; } //  done with this list
         else {                                     //  still doing this list
            if( buf->ts < ts ){
               // this is the common ordering error - where events are delayed
               //   due to congestion, and the filter puts them in amongst
               //   much later events, so their timestamp is anomalously early
               ++err[ERR_LATE_OUT];
            }
            // copy event
            while( 1 ){ // wait for space in output buffer
               wr_avail = EVENTBUFSIZE + (eventbuf_rdpos - eventbuf_wrpos);
               if( wr_avail >= buf->evlen ){ break; }
               if( arg->grif_sort_done ){ break; } else { usleep(usecs); }
            }
            ++tsevents_out;
            for(i=0; i<buf->evlen; i++){
               wrpos = eventbuf_wrpos++ % EVENTBUFSIZE;
               event_buffer[wrpos] = buf->event[i];
            }
            pthread_mutex_lock(&nxtlock); // remove from list
            buf->in_use = 0;  buf->next = NULL;
            // since we're taking from start of list, prev is always list-head
            tslot[ts_slot] = nxt;                    // update previous nxtptr
            pthread_mutex_unlock(&nxtlock);
         }
         if( (buf=nxt) == NULL ){ break; }
      }
      output_ts = ts;
   }
   arg->reorder_out_done = 1;
   printf("Reorder: end output thread [ts:%3ds] ...\n", (int)(ts/100000000) );
   reorder_status(0);
   return;
}
