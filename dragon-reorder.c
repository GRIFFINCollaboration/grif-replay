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
#include "dragon-format.h"
#include "midas-format.h"

int reorder_events_read;
volatile long tsevents_in;
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
#define REORDER_MAXEVENTSIZE (DRAGON_EVENTWORDS+8)
typedef struct reorderbuf_struct Tsbuf;
struct reorderbuf_struct {
   volatile Tsbuf *next; unsigned long ts;
   volatile int in_use;
   int event[REORDER_MAXEVENTSIZE];
};
Tsbuf  reorder_buffer[REORDER_EVENTS]; // ~100bytes ea -> 100Mbytes
volatile Tsbuf *tslot[REORDER_TSLOTS];
int reorder_bufpos; // next slot to be written

pthread_mutex_t nxtlock;
void reorder_main(Sort_status *arg)
{
   int i, rd_avail, ts_stat, ev_done, ts_slot, *err = diagnostics.reorder_error;
   unsigned int usecs=1000, *inptr, *bufend;
   volatile Tsbuf *bufptr, *nxtptr, *newptr;
   Dragon_event *evptr;
   FILE *log_fp=fopen("reorder.log","w");

   memset(err,           0, REORDER_ERRTYPES*sizeof(int));
   printf("starting reorder input ...\n");
   for(i=0; i<REORDER_EVENTS; i++){
      reorder_buffer[i].in_use = 0;  reorder_buffer[i].next = NULL;
   }
   for(i=0; i<REORDER_TSLOTS; i++){ tslot[i] = NULL; }
   tsevents_in = reorder_bufpos = ev_done = ts_stat = 0;
   inptr = bankbuf;  bufend = bankbuf + BANK_BUFSIZE;
   while(1){    evptr = (Dragon_event *)inptr;
      if( tsevents_in - tsevents_out > OVERFULL_FRACTION*REORDER_EVENTS ){
         usleep(usecs); continue; // do not over-fill buffer
      }
      if( (rd_avail = bankbuf_wrpos - bankbuf_rdpos) < DRAGON_EVENTWORDS ){
         if( arg->end_of_data ){
            break;
         }
         usleep(usecs); continue;
      }
      fprintf(log_fp, "Reorder #%6d buf[%8d] [wrpos:%8d] ... ", tsevents_in, bankbuf_rdpos, bankbuf_wrpos);
      if( evptr->begin_marker != DRAGON_EVENT_MARK ){ // wait for next start-mark
         ++bankbuf_rdpos; printf("Reorder BAD MARK\n");
         if( ++inptr >= bufend ){ inptr = bankbuf; }
         continue;
      }
         
 // bankbuf_rdpos += DRAGON_EVENTWORDS;
 // inptr += DRAGON_EVENTWORDS; if( inptr >= bufend ){ inptr -= BANK_BUFSIZE; }
      
      ++reorder_events_read;
      // events that are so late, that we have already output their timeslot
      // can no longer be made to be in order
      // - drop them for now (may include anyway - mark as just for singles?)
      if( evptr->ts < output_ts ){
         ++err[ERR_LATE_IN]; err[ERR_WORDS_IN] += DRAGON_EVENTWORDS;
         bankbuf_rdpos += DRAGON_EVENTWORDS;
         fprintf(log_fp, " LATE [ts, out_ts:%ld %ld]\n", evptr->ts, output_ts);
         ts_stat = 0; continue;
      }
      // events that are Extremely early (firmware bugs/data corruption)
      // these would stop reuse of their and subsequent buffer slots until
      // their timeslot comes (reuse would be expected ~64k events later)
      //   i.e. buffer would become blocked for a long time
      //   discard for now, until implement a way of avoiding blocked buf
      if( evptr->ts > (output_ts + TOO_EARLY_CUTOFF) ){
         ++err[ERR_EARLY_IN]; err[ERR_WORDS_IN] += DRAGON_EVENTWORDS;
         bankbuf_rdpos += DRAGON_EVENTWORDS;
         fprintf(log_fp, " EARLY [ts, out_ts:%ld %ld]\n", evptr->ts, output_ts);
         ts_stat = 0; continue;
      }
      fprintf(log_fp, " OK [ts, out_ts:%ld %ld]\n", evptr->ts, output_ts);
      bufptr = NULL;  ts_slot = (evptr->ts >> BUCKET_SIZE_BITS) % REORDER_TSLOTS;
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
      newptr->in_use =  1;
      newptr->ts     = evptr->ts; newptr->next  = NULL;
      if( inptr + DRAGON_EVENTWORDS < bufend ){
         memcpy((char *)(newptr->event),(char *)(inptr),sizeof(Dragon_event));
      } else {  // wrapped in middle of event
         i = bufend - inptr;
         memcpy((char *)(newptr->event),  (char *)(inptr),i);
         memcpy((char *)(newptr->event)+i,(char *)(bankbuf),sizeof(Dragon_event)-i);
      }
      // ONLY AFTER DATA COPIED - update readpos
      bankbuf_rdpos += DRAGON_EVENTWORDS;
      inptr += DRAGON_EVENTWORDS; if( inptr >= bufend ){ inptr -= BANK_BUFSIZE; }
      // update linked list ...
      bufptr = NULL;
      pthread_mutex_lock(&nxtlock);
      if( (nxtptr = tslot[ts_slot]) == NULL ){ // currently empty
         tslot[ts_slot] = newptr;
      } else {                              // insert into list IN TIME ORDER
         while( 1 ){ // at this point - bufptr is previous, nxtptr is current
            if( nxtptr == NULL ){
               printf("IMPOSSIBLE ERROR\n"); break;
            }
            if( evptr->ts <= nxtptr->ts || nxtptr->in_use == 0 ){ // insert here
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
      ev_done = ts_stat = 0;
   }
   arg->reorder_in_done = 1;
   printf("Reorder: end input thread\n");
   fclose(log_fp); return;
}

#define EVT_BUFSIZE         16384
Dragon_event evbuf[EVT_BUFSIZE];
volatile long evbuf_count;
volatile unsigned long evbuf_wrpos;
extern volatile unsigned long evbuf_rdpos;

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
   evbuf_count = evbuf_rdpos = evbuf_wrpos = 0;
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
               //printf("REORDER ERROR-VERY-LONG-DATA-GAP\n");
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
               wr_avail = EVT_BUFSIZE + (evbuf_rdpos - evbuf_wrpos);
               if( wr_avail >= 1 ){ break; }
               if( arg->grif_sort_done ){ break; } else { usleep(usecs); }
            }
            ++tsevents_out;
            wrpos = evbuf_wrpos % EVT_BUFSIZE;
            memcpy((char *)&evbuf[wrpos], (char *)buf->event,
                   sizeof(Dragon_event) );
            ++evbuf_wrpos; // Do this *AFTER* above copy
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
