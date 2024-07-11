// event time-ordering - scrambles original event ordering
//    => midas event/bank headers no longer connected to their events
//         (but there is nothing useful there anyway)
// we used to use midas timestamps, but event timestamps are more useful
// ----------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "grif-replay.h"
#include "midas-format.h"

#define MAX_GRIFC 16

unsigned event_buffer[EVENTBUFSIZE+128];
volatile unsigned long eventbuf_rdpos;
volatile unsigned long eventbuf_wrpos;

// input ...
//    minimum space is #grifcs * #events per grifc buffer
//    max packet size is 9k = 300 per grifc times 4 grifcs = 1200 MINIMUM
// try 16k events and 16k slots
// #define TIMESLOT_SLOTS     262144
// #define TIMESLOT_EVENTS     81920
#define TIMESLOT_SLOTS     16384
#define TIMESLOT_EVENTS    16384
#define TIMESLOT_EVENTSIZE    20 // max 20 words - 80 bytes
typedef struct reorderbuf_struct Tsbuf;
struct reorderbuf_struct {
   volatile Tsbuf *next; unsigned long ts;
   volatile int in_use; int evlen;
   int event[TIMESLOT_EVENTSIZE];
};
Tsbuf timeslot_buffer[TIMESLOT_EVENTS];
volatile Tsbuf *tslot[TIMESLOT_SLOTS];
int tsbufpos; // next slot to be written
volatile long tsevents_in;
int reorder_events_read;

// output ...
#define TS_EVENT_BUFFER 8192
volatile unsigned long output_ts;
pthread_mutex_t nxtlock;
long tsevents_out;

// read events out of bank buffer and copy to appropriate reorder_buffer
//               [each grifc has its own buffer]
//     filters in each grifc keeps events ordered within that module
// [NOTE: grifc filter can fail if any data is delayed more than ~150us]
#define REORDER_TSBUFSZ (REORDER_BUFSIZE>>3)
typedef struct Reorder_buf_struct {
   int     data[REORDER_BUFSIZE]; 
   long   tsbuf[REORDER_TSBUFSZ];
   long evtspos[REORDER_TSBUFSZ]; // position in databuf of each events ts
   long  hdrpos[REORDER_TSBUFSZ]; // position in databuf of each events hdr
   volatile long    wrpos;  volatile long    rdpos;   long tmp_wrpos;  
   volatile long ts_wrpos;  volatile long ts_rdpos;   int pending_events;
   volatile int    events;  volatile int  events_done;
} Reorder_buffer;
Reorder_buffer *reorder_buf[MAX_GRIFC];
volatile int reorder_ready_count;

// bank buffer is 64Mbytes
volatile int reorder_active[MAX_GRIFC];
volatile int reorder_evcount[MAX_GRIFC];
volatile int reorder_discard[MAX_GRIFC];
int reorder_init[MAX_GRIFC];
int reorder_stat[4*MAX_GRIFC];

// single-thread-debug
//   or known file

////////////////////////////////////////////////////////////////////////////
///////////        check-buffer - debugging function              //////////
///////////  (remove when all initial buffer errors are fixed)    //////////
////////////////////////////////////////////////////////////////////////////
#define STATE_HDR  1
#define STATE_TSA  2
#define STATE_TSB  3
#define STATE_TRLR 4
int check_buffer(int grifc)
{
   unsigned long i, ts, tsa, tsb, tmp, evstart, tspos;
   int wcnt, state, err, wrpos, rdpos, tscnt, evcount;
   Reorder_buffer *buf = reorder_buf[grifc];
   unsigned word, type, rel_tspos;

   pthread_mutex_lock(&nxtlock);
   
   state = STATE_HDR;
   err = ts = tsa = tsb = evcount = wcnt = 0;
   tscnt = buf->ts_wrpos - buf->ts_rdpos;
   if( tscnt - buf->pending_events != buf->events ){ err = 1;
      printf("Check Buffer[%d] - Event count wrong\n", grifc);
      printf("   ts_wr:%ld rs_rd:%ld diff:%d pend:%d count:%d\n", buf->ts_wrpos, buf->ts_rdpos, tscnt, buf->pending_events, buf->events);
   }
   tspos = buf->ts_rdpos;
   rel_tspos = tspos     % REORDER_TSBUFSZ;
   rdpos = buf->rdpos    % REORDER_BUFSIZE;
   wrpos = buf->wrpos    % REORDER_BUFSIZE;
   for(i=buf->rdpos; i!=buf->wrpos; i++){
      rdpos = i % REORDER_BUFSIZE; wcnt = i-buf->rdpos;
      word = buf->data[rdpos]; type = word >> 28;
      switch(type){
      case 0: case 1: case 2: case 3: case 0xc: case 0x9:
      case 4: case 5: case 6: case 7: case 0xd: break;
      case 0x8:
         if( state != STATE_HDR ){ err = 1;
            printf("     [%d]  ERR Event %d, word %d[pos:%d][0x%08x]: Repeated Header\n", grifc,
                      evcount, wcnt, rdpos, word);
         }
         tmp = buf->hdrpos[rel_tspos] % REORDER_BUFSIZE;
         if( buf->hdrpos[rel_tspos] != i ){ err = 1;
            printf("     [%d]  ERR Event %d, pos:%ld[%d]: Wrong Hdrpos:%ld[%d]\n", grifc,
                   evcount, i, rdpos, buf->hdrpos[rel_tspos], tmp );
         }
         evstart = i; state = STATE_TSA; break;
      case 0xa:
         if( state != STATE_TSA ){ err = 1;
            printf("     [%d]  ERR Event %d, word %d[pos:%d][0x%08x]: tsa at state %d\n", grifc,
                   evcount, wcnt, rdpos, word, state);
         } tsa = word & 0xfffffff; state = STATE_TSB;
         if( (buf->tsbuf[rel_tspos] & 0xfffffff) != tsa ){ err = 1;
            printf("     [%d]  ERR Event %d: tsa mismatch buf:%x ev:%x\n", grifc, evcount,
                (unsigned)(buf->tsbuf[rel_tspos] & 0xfffffff), (unsigned)tsa );
         } break;
      case 0xb:
         if( state != STATE_TSB ){ err = 1;
            printf("     [%d]  ERR Event %d, word %d[0x%08x]: tsb at state %d\n", grifc,
                   evcount, wcnt, rdpos, word, state);
         }
         if( buf->evtspos[rel_tspos] != i ){ err = 1;
            printf("     [%d]  ERR Event %d, word %d[pos:%d]: wrong tspos:%d[pos:%d]\n", grifc,
                   evcount, wcnt, rdpos, word, buf->evtspos[rel_tspos], buf->evtspos[rel_tspos] % REORDER_BUFSIZE );
         }
         tsb = word & 0x3fff; state = STATE_TRLR;
         tmp = buf->tsbuf[rel_tspos] & 0x1ffffffffff; // mask out-of-order bit
         if( (tmp>>28) != tsb ){ err = 1;
            printf("     [%d]  ERR Event %d: ts mismatch buf:%x ev:%lx\n", grifc, evcount,
                   (unsigned)(buf->tsbuf[rel_tspos]>>28), tsb );
         } break;
      case 0xe:
         if( state != STATE_TRLR ){ err = 1;
            printf("     [%d]  ERR Event %d, word %d[pos:%d][0x%08x]: trlr at state %d\n", grifc,
                   evcount, wcnt, rdpos, word, state);
         }
         ts = tsa + (tsb<<28);
         tmp = buf->tsbuf[rel_tspos] & 0x1ffffffffff; // mask out-of-order bit
         if( tmp != ts ){ err = 1;
            printf("     [%d]  ERR Event %d: ts mismatch buf:%lx ev:%lx\n",
                   grifc, evcount, buf->tsbuf[rel_tspos], ts);
         }
         if( 1 || evcount < 2 || evcount >= buf->events-2 ){
            printf("     [%d]  Event %3d{rdpos:%8ld[%6d] to %8ld[%6d]}{ts_rdpos:%7ld[%6d]}, hdrpos:%7ld[%6d] evtspos:%7ld[%6d]:0x%015x\n", grifc, evcount, evstart, (unsigned)(evstart % REORDER_BUFSIZE), i, (unsigned)(i % REORDER_BUFSIZE), tspos, rel_tspos, buf->hdrpos[rel_tspos],  (unsigned)(buf->hdrpos[rel_tspos] % REORDER_BUFSIZE), buf->evtspos[rel_tspos], (unsigned)(buf->evtspos[rel_tspos] % REORDER_BUFSIZE), tmp );  fflush(stdout);
         }
         ++evcount; ++tspos; rel_tspos = tspos % REORDER_TSBUFSZ;
         state = STATE_HDR; break;
      default:
          printf("     [%d]  ERR Event %d, word %d[0x%08x]: ? word at state %d\n", grifc,
                 evcount, wcnt, rdpos, word, state);  err = 1; break;
      }
   }
   if( evcount != buf->events ){ err = 1;
      printf("     [%d]  ERR Event count mismatch - event-check:%d claimed:%d\n", grifc,
             evcount, buf->events );
   }
   printf("CHKBF[%d] [Completed:%d] DataEvents:%d pending:%d claim:%d rdpos:%ld[%d] wrpos:%ld[%d] tmp_wrpos:%ld[%d] {tsrdpos:%ld[%d] tswrpos:%ld[%d]}\n", grifc, buf->events_done, evcount, buf->pending_events, buf->events, buf->rdpos,  buf->rdpos % REORDER_BUFSIZE, buf->wrpos,  buf->wrpos % REORDER_BUFSIZE, buf->tmp_wrpos,  buf->tmp_wrpos % REORDER_BUFSIZE, buf->ts_rdpos,  buf->ts_rdpos % REORDER_TSBUFSZ, buf->ts_wrpos,  buf->ts_wrpos % REORDER_TSBUFSZ ); fflush(stdout);
   pthread_mutex_unlock(&nxtlock);
   if( err == 0 ){
      return(0);
   }
   return(1);
}

extern Sort_metrics diagnostics;
void reorder_status(int current_time)
{
   int sum, *err = diagnostics.reorder_error;
   sum = err[0] + err[2] + err[3] + err[4] + err[5] + err[6] + err[7];// err[8] is in addition to other
                                                                      // errors - do not add to total
   printf("Reorder: in:%10d out:%10d err:%10d[%5.1f%%] [Desync:%d]\n         ",
          tsevents_in, tsevents_out, sum, (100.0*sum)/tsevents_in, err[8] );
   printf("[init:%d format:%d length:%d early:%d late:%d,%d, unk:%d]\n",
          err[2], err[0], err[3], err[5], err[4], err[6], err[7] );
}

// events are temporarily written into buffer, using local wrpos
// if no errors, global tmp_wrpos is updated to include latest event
// THEN after timestamp sequence check, final wrpos is updated to include evt
void reorder_a_main(Sort_status *arg)
{
   int grifc, evlen, rd_avail, wr_avail, tshi_pos, *err = diagnostics.reorder_error;
   int i, tmp, wrpos, type, bad_blk, bad_cnt, fmterr, *evptr, *bufend;
   unsigned int usecs=100;
   Reorder_buffer *buf;
   long ts1, ts2;
   static int wrap;

   for(i=0; i<MAX_GRIFC; i++){
      if( (reorder_buf[i] = (Reorder_buffer *)calloc(1,sizeof(Reorder_buffer))) == NULL ){
         fprintf(stderr,"Can't alloc space for reorder buffer\n");
         arg->reorder_in_done = 1; return;
      }
   } 
   //memset((char *)reorder_buf,     0, MAX_GRIFC*sizeof(Reorder_buffer));

   printf("starting reorder input ...\n");
   reorder_ready_count = 0;
   evptr = bankbuf; bufend = bankbuf + BANK_BUFSIZE;
   memset((char *)reorder_active,  0, MAX_GRIFC*sizeof(int));
   memset((char *)reorder_discard, 0, MAX_GRIFC*sizeof(int));
   memset((char *)reorder_evcount, 0, MAX_GRIFC*sizeof(int));
   memset((char *)reorder_init,    0, MAX_GRIFC*sizeof(int));
   memset((char *)reorder_stat,    0, MAX_GRIFC*sizeof(int)*4);
   bad_blk = bad_cnt = fmterr = evlen = wrap = 0;  grifc = -1;
   // to avoid annoying complications, all events are checked for correct
   // format here, and obviously bad events are removed
   //
   // events -> output buffer, but if event is bad, do not confirm writes
   // and then overwrite bad event with next one
   //
   // wrpos is temp write-position, reorder_wrpos[grifc] is confirmed posn
   while(1){
      if( arg->shutdown_midas != 0 ){ break; }
      rd_avail = bankbuf_wrpos - bankbuf_rdpos;
      if( rd_avail < 1 ){ usleep(usecs);
         if( arg->end_of_data ){ break; } else { continue; }
      }
      if( (type = ((*evptr >> 28) & 0xf)) == 0x8 ){
         if( grifc != -1 ){ // still writing prev event - missing trailer
            // dump old event (i.e. just start a new one here)
            ++reorder_events_read;    fmterr = wrap = evlen = 0;
            ++reorder_evcount[grifc]; grifc = -1;  ++tsevents_in; buf = NULL;
            // printf("DUMP [?] ts=0x%04x,0x%08x[@????] TRLR:0x%08x[@??] *No-Trailer/Rpt-Header*\n", ts2, ts1, *(int *)evptr );
            ++err[ERR_FORMAT_IN];
         }
         ts1 = ts2 = 0; tshi_pos = -1;
         grifc = (*evptr>>16) & 0xf; buf = reorder_buf[grifc];
         // if( grifc == 15 ){ printf("GRIFC:0x%08x\n", *(int *)evptr ); }
         //   ^ These are the grifc filter/link-status scalar events
         wrpos = buf->ts_wrpos %  REORDER_TSBUFSZ;
         buf->hdrpos[wrpos] = buf->tmp_wrpos;
         wrpos = buf->tmp_wrpos %  REORDER_BUFSIZE;
      } else if( grifc == -1 ){ // no header[unknown grifc] - have to discard
         ++bad_blk;
         ++bankbuf_rdpos;
         ++evptr; if( evptr >= bufend){ evptr -= BANK_BUFSIZE; }
         while(1){ // REMAIN IN THIS LOOP UNTILL NEXT EVENT_HEADER
            rd_avail = bankbuf_wrpos - bankbuf_rdpos;
            if( rd_avail < 1 ){ usleep(usecs);
               if( arg->end_of_data ){ break; } else { continue; }
            }
            if( (type = ((*(evptr) >> 28) & 0xf)) == 0x8 ){
               grifc = ((*(evptr) >> 16) & 0xf); buf = reorder_buf[grifc];
               ++reorder_events_read; break;
            }
            if( type == 0xA ){ ts1 = *(evptr) & 0xfffffff; 
            } else if( type == 0xB ){ ts2 = *(evptr) & 0x3fff; 
            } else if( type == 0xE ){
               //printf("DUMP [?] ts=0x%04x,0x%08x[@????] TRLR:0x%08x[@??] missing header\n", ts2, ts1, *(int *)evptr );
               ++err[ERR_FORMAT_IN];
            }
            ++bankbuf_rdpos; ++bad_cnt;
            ++evptr; if( evptr >= bufend){ evptr -= BANK_BUFSIZE; }
         }
      } else if( type == 0xA ){ fmterr += ts1!=0; ts1 = *(evptr) & 0xfffffff; 
      } else if( type == 0xB ){ fmterr += ts2!=0; ts2 = *(evptr) &    0x3fff;
         tshi_pos = buf->tmp_wrpos + evlen;
      }
      while( 1 ){ // wait for space in reorder buffer
         wr_avail = REORDER_BUFSIZE - evlen -
                   (buf->tmp_wrpos - buf->rdpos);
         if( wr_avail != 0 ){ break; }
         usleep(usecs);
      }
      if( type != 0xC ){ // strip waveforms - buffers too small for these atm.
         buf->data[wrpos] = *evptr; ++evlen;
         if( ++wrpos >= REORDER_BUFSIZE ){ wrpos -= REORDER_BUFSIZE; wrap=1; }
      }
      if( ++evptr >= bufend){ evptr -= BANK_BUFSIZE; }
      ++bankbuf_rdpos; reorder_active[grifc] = 1;
      if( type == 0xE ){
         // due to junk data at BOR - dump up to 1k events per grifc
         // *AND* temp during development - minimum of 50 events
         if( reorder_init[grifc] == 0 ){
            if( (ts2 == 0 && reorder_discard[grifc] > 50) ||
                         ++reorder_discard[grifc] >= 1000 ){
               reorder_init[grifc] = 1;
            } else {
               ++reorder_stat[grifc+2*MAX_GRIFC];
               reorder_stat[grifc+3*MAX_GRIFC] += evlen;
               ++err[ERR_INIT_IN];
               fmterr = 99;
            }
         }
         if( fmterr == 0 && ts1 != 0 && tshi_pos != -1 && ts1 != 2 ){
            buf->tmp_wrpos += evlen;         // update tmp_wrpos
            tmp = (buf->tmp_wrpos-1) % REORDER_BUFSIZE;
            //pthread_mutex_lock(&nxtlock);
            //printf("STORE[%d] ts=0x%04x,0x%08x[Into:%4ld] TRLR:0x%08x[@%6ld] len=%2d  ", grifc, ts2, ts1, buf->ts_wrpos, *(int *)(evptr-1), tmp, evlen ); fflush(stdout);
            //for(i=0; i<evlen; i++){
            //  tmp = (buf->tmp_wrpos - evlen + i) % REORDER_BUFSIZE;
            //   printf("0x%08x ", buf->data[tmp] );
            //} printf("\n"); fflush(stdout);
            //pthread_mutex_unlock(&nxtlock);
            wrpos = buf->ts_wrpos % REORDER_TSBUFSZ; // store ts
            buf->tsbuf[wrpos] = ts1 + (ts2<<28);
            buf->evtspos[wrpos] = tshi_pos;
            
            pthread_mutex_lock(&nxtlock);
            ++buf->ts_wrpos; ++buf->pending_events; // do these together - do not allow to be split
            pthread_mutex_unlock(&nxtlock);

            reorder_insert_event_check(grifc);
         } else {
            if( fmterr == 0 ){ ++err[ERR_FORMAT_IN]; } // error not already counted
            // leave reorder_wrpos as is ( => do not keep this event)
            //printf("DUMP-[%d] ts=0x%04x,0x%08x[Into:%4ld] TRLR:0x%08x[@%6ld] len=%2d %s\n", grifc, ts2, ts1, buf->ts_wrpos, *(int *)(evptr-1), wrpos-1, evlen, fmterr==99?"INIT":"fmt-error" );
         }
         ++reorder_events_read;    fmterr = wrap = evlen = 0;
         grifc = -1;  ++tsevents_in;  buf = NULL;
      }
   }
   arg->reorder_in_done = 1;
   printf("END reorder input ... %d words discarded in %d blocks\n", bad_cnt, bad_blk);
   return;
}
// switch chk to read thread: inc pending in wr_thread and check in rd_thread

#define OUT_OF_ORDER_BIT (0x20000000000)

// each event inserted into reorder buffer - check timestamp sequence
//   for out-of-order events, and mark them
//   also finalize event write
#define NUM_ORDER_EVENT 7
#define EVENT_CHCK ( (unsigned)0xffffffff )
static int insert_event;
// if have sequence of N in-order events, finalize write of first of them
int reorder_insert_event_check(int grifc) // ### ALL CURRENTLY LOCKED ###
{
   int i, pending, tmp, ins_err, evrdpos, tsrdpos, wrpos, wr_avail, evlen, first;
   int early_ts[NUM_ORDER_EVENT], *err = diagnostics.reorder_error;
   Reorder_buffer *buf = reorder_buf[grifc];
   long ts_list[NUM_ORDER_EVENT];

   // currently this is called by reorder_in (which modifies *wrpos)
   // => only rdpos can be changed by other thread during this
   //    record pending and all rdpos in locked section
   //    and use the local copies after that
   // if rdpos* change during this, does not matter,
   // as wrpos updates (possibly overwriting data) are controlled from here
   //
   // ABOVE IRRELEVANT ...
   //    the only processing done here is between wrpos and tmp_wrpos
   //    (do not even touch rdpos, except to check if buf now half full)
   // 
   pthread_mutex_lock(&nxtlock);
   pending = buf->pending_events;
   pthread_mutex_unlock(&nxtlock);
   //if( ++insert_event > EVENT_CHCK ){ check_buffer(grifc); }
   if( pending < NUM_ORDER_EVENT ){ /*check_buffer(grifc); */ return(0); }
   
   ins_err = 0; first = -1;
   for(i=0; i<NUM_ORDER_EVENT; i++){ early_ts[i] = 0;
      wrpos = (buf->ts_wrpos - pending + i) % REORDER_TSBUFSZ;
      ts_list[i] = buf->tsbuf[wrpos] & 0x1ffffffffff;
      if( i > 0 && ts_list[i] < ts_list[i-1] ){
         ins_err = 1; early_ts[i] = 1; if( first == -1 ){ first = i; }
         if( first < 3 ){
            ins_err = 1;
         }
      }
   }
   if( ins_err == 0 || first > 3 ){
      // seven events in buffer => nxt events hdrpos is already set
      wrpos = (buf->ts_wrpos - pending) % REORDER_TSBUFSZ;
      evlen = buf->hdrpos[(wrpos+1) % REORDER_TSBUFSZ] - buf->hdrpos[wrpos];
      
      pthread_mutex_lock(&nxtlock); // use current rdpos in case updated
      wr_avail =  REORDER_BUFSIZE - (buf->wrpos - buf->rdpos);
      if( wr_avail>=REORDER_BUFSIZE/2 && wr_avail-evlen < REORDER_BUFSIZE/2 ){
         --reorder_ready_count;  // just passed half-way mark
      }
      buf->wrpos += evlen; --buf->pending_events; ++buf->events; 
      ++reorder_evcount[grifc];
      pthread_mutex_unlock(&nxtlock);
      
      if( insert_event > EVENT_CHCK ){
         pthread_mutex_lock(&nxtlock);
         tmp = (buf->wrpos-1) % REORDER_BUFSIZE;
         printf("FINLZ[%d] ts=0x%015x[In  :%4ld] TRLR:0x%08x[@%6ld] len=%2d  ", grifc, buf->tsbuf[wrpos], wrpos, buf->data[tmp], tmp, evlen );
         for(i=0; i<evlen; i++){
            tmp = (buf->wrpos - evlen + i) % REORDER_BUFSIZE;
            printf("0x%08x ", buf->data[tmp] );
         } printf("\n"); fflush(stdout);
         pthread_mutex_unlock(&nxtlock);
      }
   } else {
      // figure out which ones are wrong, mark them,
      // and finalize up to and including all wrong events
      // two cases ...
      //    usual, where event is delayed and timestamp is far in past ...
      //       ts       100 101 102 3 103 104 105
      //       early     0   0   0  1  0   0   0
      //       first seen when ts=3 appears
      // 
      //    problem, where event comes early and timestamp is far in future
      //       ts:       100 101 203 102 103 104 105
      //       early:     0   0   0   1   0   0   0
      //       first noticed when event AFTER the wrong one appears
      //
      // now have 7 events, with first error at position 3 (4th item)

      // just 4th item is delayed ...
      if( ts_list[1] >= ts_list[0] && ts_list[2] >= ts_list[1] &&
          /* ts_list[3] is wrong */   ts_list[4] >= ts_list[2] &&
          ts_list[5] >= ts_list[4] && ts_list[6] >= ts_list[5] ){
         ++err[ERR_LATE_IN];
         wrpos = (buf->ts_wrpos - pending) % REORDER_TSBUFSZ;
         buf->         tsbuf[(wrpos+3) % REORDER_TSBUFSZ] |= OUT_OF_ORDER_BIT;
         evlen = buf->hdrpos[(wrpos+4) % REORDER_TSBUFSZ] - buf->hdrpos[wrpos];
         
         pthread_mutex_lock(&nxtlock); // do b4 changing wrpos
         wr_avail = REORDER_BUFSIZE - (buf->wrpos - buf->rdpos);
         if( wr_avail >= REORDER_BUFSIZE/2 && wr_avail-evlen <
             REORDER_BUFSIZE/2 ){ --reorder_ready_count; }
         buf->wrpos += evlen; buf->pending_events -= 4; buf->events += 4; 
         reorder_evcount[grifc] += 4;
         pthread_mutex_unlock(&nxtlock);
         
         if( insert_event > EVENT_CHCK ){ printf("FINLZ[%d] 4 EVENTS\n", grifc); }
      } else if( 
      // just 3rd item is early ...
            ts_list[1] >= ts_list[0] && /* ts_list[2] is wrong */
            ts_list[3] >= ts_list[1] && ts_list[4] >= ts_list[3] &&
            ts_list[5] >= ts_list[4] && ts_list[6] >= ts_list[5] ){
         ++err[ERR_EARLY_IN];
         wrpos = (buf->ts_wrpos - pending) % REORDER_TSBUFSZ;
         buf->         tsbuf[(wrpos+2) % REORDER_TSBUFSZ] |= OUT_OF_ORDER_BIT;
         evlen = buf->hdrpos[(wrpos+3) % REORDER_TSBUFSZ] - buf->hdrpos[wrpos];
         
         pthread_mutex_lock(&nxtlock); // do b4 changing wrpos
         wr_avail = REORDER_BUFSIZE - (buf->wrpos - buf->rdpos);
         if( wr_avail >= REORDER_BUFSIZE/2 && wr_avail-evlen <
             REORDER_BUFSIZE/2 ){ --reorder_ready_count; }
         buf->wrpos += evlen; buf->pending_events -= 3; buf->events += 3; 
         reorder_evcount[grifc] += 3;
         pthread_mutex_unlock(&nxtlock);
         
         //if( insert_event > EVENT_CHCK ){ printf("FINLZ[%d] 3 EVENTS\n", grifc); }
      } else {
         //if( insert_event > EVENT_CHCK ){ printf("UnHnd[%d] [Ev:%6d]: ", grifc, insert_event );
         //   for(i=0; i<NUM_ORDER_EVENT; i++){ printf("%d ", early_ts[i]); }
         //   printf("     T-T0:["); for(i=0; i<NUM_ORDER_EVENT; i++){
         //      printf("%d ", ts_list[i]-ts_list[0]);
         //   }
         //   printf("]     Dt:[XXXX "); for(i=1; i<NUM_ORDER_EVENT; i++){
         //      printf("%d ", ts_list[i]-ts_list[i-1]);
         //   }
         //   printf("]\n");
         //}
         
         // Currently un-handled case, can easily add handling of other common errors
         err[ERR_UNKNOWN_ORDER_IN] += (buf->pending_events-1);
         
         wrpos = (buf->ts_wrpos - pending) % REORDER_TSBUFSZ;
         for(i=0; i<pending; i++){
            buf->tsbuf[ (wrpos+i) % REORDER_TSBUFSZ ] |= OUT_OF_ORDER_BIT;
         }
         evlen = buf->hdrpos[(wrpos+i-1)%REORDER_TSBUFSZ] - buf->hdrpos[wrpos];
         
         pthread_mutex_lock(&nxtlock); // do b4 changing wrpos
         wr_avail = REORDER_BUFSIZE - (buf->wrpos - buf->rdpos);
         if( wr_avail >= REORDER_BUFSIZE/2 && wr_avail-evlen <
             REORDER_BUFSIZE/2 ){ --reorder_ready_count; }
         buf->wrpos += evlen; buf->events += (buf->pending_events-1);
         reorder_evcount[grifc] += (buf->pending_events-1);
         buf->pending_events = 1;
         pthread_mutex_unlock(&nxtlock);
         
         //if( insert_event > EVENT_CHCK ){ printf("FINLZ[%d] %d EVENTS\n", grifc, buf->pending_events-1 ); }
      }
   }
   //if( insert_event > EVENT_CHCK ){ check_buffer(grifc); }
   //if( insert_event % 10000 == 0 ){ printf("Insert:%d\n", insert_event ); }
   return(0);
}

#define EV_COPY     0 // must be zero
#define ERR_INIT    1
#define ERR_RECOVER 2
#define ERR_TS_HIGH 3
int copy_event(Sort_status *arg, int grifc, int discard)
{
   int type, ts_pos, evlen, mis, *err = diagnostics.reorder_error;
   int i, evwrd, rdpos, wrpos, rd_avail, wr_avail;
   Reorder_buffer *buf = reorder_buf[grifc];
   unsigned long ts, tsa, tsb;
   unsigned int usecs=100;
       
   //printf("## Storing event [grifc#%d]:", grifc);

   if( grifc == 11 ){
     evlen = 1;
   }

   evlen = i = 0;
   while( 1 ){ // copy (or discard) whole event
      rd_avail = buf->wrpos - (buf->rdpos + evlen);
      if( rd_avail == 0 ){
	 if( ++i >= 100 ){
 	    fprintf(stderr,"REORDER BUFFER CORRUPTION\n");
            pthread_mutex_lock(&nxtlock);
	    buf->rdpos = buf->wrpos;  buf->ts_rdpos = buf->ts_wrpos;
            buf->pending_events = 0;
            pthread_mutex_unlock(&nxtlock);
            return(0);
         }
         if( arg->shutdown_midas != 0 || arg->end_of_data ){ break; }
         usleep(usecs); continue;
      }
      i = 0; rdpos = (buf->rdpos + evlen++) % REORDER_BUFSIZE;
      evwrd = buf->data[rdpos];
      if( !discard ){
         while( 1 ){ // wait for space in output buffer
            wr_avail = EVENTBUFSIZE + eventbuf_rdpos - eventbuf_wrpos;
            if( wr_avail != 0 ){ break; }
            usleep(usecs);
         }
         wrpos = eventbuf_wrpos % EVENTBUFSIZE;
         event_buffer[wrpos] = evwrd;
         ++eventbuf_wrpos;
      }
      if( (type = (evwrd >> 28) & 0xF) == 0xE ){
         // check for timestamp buffer misalignment
         ts_pos = buf->ts_rdpos % REORDER_TSBUFSZ;
         // ts = tsa + (tsb<<28); ERROR:tsb is unsigned => upper bits lost
         ts = (tsa + (tsb<<28)) & 0x1ffffffffff;
         if( ( buf->tsbuf[ts_pos] & 0x1ffffffffff ) == ts ){
            mis = 0; if( insert_event > EVENT_CHCK ){ /*printf(".");*/ }
         } else {
            ts_pos = (buf->ts_rdpos+1) % REORDER_TSBUFSZ;
            if( ( buf->tsbuf[ts_pos] & 0x1ffffffffff ) == ts ){
               pthread_mutex_lock(&nxtlock);
               ++buf->ts_rdpos; --buf->events;
               pthread_mutex_unlock(&nxtlock);
            } else {
               ts_pos = (buf->ts_rdpos-1) % REORDER_TSBUFSZ;
               if( ( buf->tsbuf[ts_pos] & 0x1ffffffffff ) == ts ){
                  pthread_mutex_lock(&nxtlock);
                  --buf->ts_rdpos; ++buf->events;
                  pthread_mutex_unlock(&nxtlock);
               } else {
                  mis = 1; ++err[ERR_TS_DESYNC_IN];
               }
            }
            //pthread_mutex_unlock(&nxtlock);
            //check_buffer(grifc);
            //pthread_mutex_lock(&nxtlock);
            // printf("#");
            if( grifc == 1 ){
               mis = 1;
            }
         }
         //if( insert_event > EVENT_CHCK ){
         //   check_buffer(grifc);
         //   pthread_mutex_lock(&nxtlock);
         //   printf("%s%s[%d] ts=0x%015lx[Out :%4ld] TRLR:0x%08x[@%6ld] len=%2d EvTS=0x%015x  ", mis?"BAD":"OK_", discard?"Dp":"Cp", grifc, buf->tsbuf[ts_pos], ts_pos, evwrd, rdpos, evlen, ts);
         //   for(i=0; i<11; i++){
         //      rdpos = (buf->rdpos - evlen + i) % REORDER_BUFSIZE;
         //      printf("0x%08x ", buf->data[rdpos] );
         //   } printf("\n"); fflush(stdout);
         //   pthread_mutex_unlock(&nxtlock);
         //}
         pthread_mutex_lock(&nxtlock);
         ++buf->ts_rdpos; buf->rdpos += evlen;;
         --reorder_evcount[grifc]; --buf->events; ++buf->events_done;
         if( buf->wrpos -  buf->rdpos          <  REORDER_BUFSIZE/2 &&
             buf->wrpos - (buf->rdpos - evlen) >= REORDER_BUFSIZE/2 ){ --reorder_ready_count; }
         pthread_mutex_unlock(&nxtlock);
         break;
      } else if( type == 0xA ){ tsa  = evwrd & 0xfffffff;
      } else if( type == 0xB ){ tsb  = evwrd & 0x3fff;
      }
   }
   return(evlen);
}

// allow config to specify smaller value of MAX_GRIFC - often only three
void reorder_a_out(Sort_status *arg)
{
   int i, j, k, ts_id,max_id,inc,type,words,junk,grifc,rdpos,event;
   unsigned long ts, prev_ts, buf_ts, max_ts, diff;
   unsigned int usecs=100, lim=10000000;
   Reorder_buffer *buf;
   
   printf("starting reorder output ...\n");
   while( reorder_buf[MAX_GRIFC-1] == NULL ){ usleep(usecs);  }
   junk = 1; prev_ts = grifc = event = 0;
   while(1){
      if( reorder_ready_count == 0 ){ // require >=1 half-full buffer
         if( !arg->reorder_in_done ){ // while waiting - check count is good
            for(i=0; i<MAX_GRIFC; i++){ buf = reorder_buf[i];
               pthread_mutex_lock(&nxtlock);
               if( buf->wrpos - buf->rdpos >= REORDER_BUFSIZE/2 ){
                  ++reorder_ready_count; }
               pthread_mutex_unlock(&nxtlock);
            }
            if( reorder_ready_count ){ printf("reorder_ready_count BAD\n"); }
            usleep(usecs); continue;
         }
      }
      // find lowest timestamp of non-empty buffers - start where left off
      //      ts = 0; for(i=grifc+1; i!=grifc; i++){
      //         if( i >= MAX_GRIFC ){ i = 0; if( grifc == 0 ){ break; } }
      ts = max_ts = -1;
      for(i=0; i<MAX_GRIFC; i++){ if( !reorder_active[i] ){ continue; }
         buf = reorder_buf[i];
         pthread_mutex_lock(&nxtlock);
         j = buf->ts_wrpos - buf->ts_rdpos - buf->pending_events;
         k = buf->wrpos    - buf->rdpos;
         pthread_mutex_unlock(&nxtlock);
         if( j == 0 || k < 9 ){ reorder_active[i] = 0; continue; }
         rdpos = buf->ts_rdpos % (REORDER_BUFSIZE/8);
         buf_ts = buf->tsbuf[rdpos];
         if( buf_ts & OUT_OF_ORDER_BIT ){ // do not change prev_ts
            copy_event(arg, i, EV_COPY); ++event; ++tsevents_out;
            --i; continue;
         }
         if( ts != -1 ){// check for 60 second or more difference in timestamps
            diff = buf_ts > ts ? buf_ts - ts : ts - buf_ts;
            if( diff > 6000000000 ){ // 6 billion
               printf("@");
               if(  buf_ts > ts ){ // buf_ts is probably wrong - dump it
                  copy_event(arg, i, EV_COPY); ++event; ++tsevents_out;
                  --i; continue;
               } else { // saved ts is probably wrong - dump that
                  // this section will only trigger if only bad saved timestamps
                  // exist - if any good timestamps are saved, the first bad
                  // timestamp will trigger the above (buf_ts > ts) section
                  // AND after saving one or more bad timestamps, the first
                  // good timestamp to appear, will kick them all out here
                  copy_event(arg, ts_id, EV_COPY); ++event; ++tsevents_out;
                  --i; max_ts = ts = -1; continue;
               }
            }
         }
         if( ts     == -1 || buf_ts <     ts ){ ts = buf_ts; ts_id = i; }
         if( max_ts == -1 || buf_ts > max_ts ){ max_ts = buf_ts; max_id = i; }
      }
      if( ts == -1 ){ usleep(usecs);  // no full events yet
         if( arg->reorder_in_done ){ break; } else { continue; }
      }
      copy_event(arg, ts_id, EV_COPY); prev_ts = ts; ++event; ++tsevents_out;
   }
   arg->reorder_out_done = 1;
   printf("END reorder output ...\n");
   for(i=0; i<MAX_GRIFC; i++){
      printf("Grifc%2d - Early[%6d] Late[%6d] BadFormat[%6d:%6dWords]\n", i,
             reorder_stat[i+0*MAX_GRIFC], reorder_stat[i+1*MAX_GRIFC],
             reorder_stat[i+2*MAX_GRIFC], reorder_stat[i+3*MAX_GRIFC]  );
   }
   return;
}

// typical event ... 10 or 11 words
//      0x82a00000 0xd000da7f 0x00010001 0x0000ef97
//      0x90007938 0xa907c044 0xb01e0068
//      0x0405a8a9 0x2f3c03a0 0x00010000 0xee18f938

// bad event caused immediate hang ..
//     0x82a0a067	0xdec193e5	0x00010001	0x76f65618
//     0xefffffff	0x82a1e0b7	0xdec193e6	0x00010001
//     0x00016e04	0x282f1bc0	0x00010000	0xe3943060

// truncated event caused hang in dump-event - stuck waiting for more data
//   midas file was on different grifc, stuck on full buffer, so deadlock
//   0x8280b011  0x00000001  0x9345504d  0xa7832878  0x5ef28746
//   0x00010000	 0x8280c0a7  0x00000001  0x909b81fa  0xa78320a9
//   0x283208a6  0x00010000
// now check trailer is already present before dumping
// add idle flag to buffer - to prevent reapeat checking

///////////////////////////////////////////////////////////////////////////
///////////////////       TIME INDEXING METHOD        /////////////////////
///////////////////////////////////////////////////////////////////////////
// firmware uses 160us reordering window - 16k slots [14bits]
// here will start with 650us              65k slots [16bits]
//
// divide into 8 sections - enforce gaps between read and write
//  [write must be >= 2 sections ahead of read]
//   => read pauses if reach section before write
//      Also write pauses if would get ahead of read
// read increases continually, but write jumps around
//   at high rates - can clear write_busy after few thousand events
//   at low rates - e.g. 100hz - 10k us gap between events
// *** => section scheme will not work ***
//
// can write multiple wraps of buffer, as store whole events,
// and can check full timestamps
//
// also used buffer slots can record if have been emptied - can pause write
// until read catches up => Write control simple
//
// how to control read?
//    keep event counts - if write 2kevents ahead or end-of-data - continue
//    [make buffer ~10k events]

#define INIT_WAIT 250 // allow this many junk events per grifc at run start
void reorder_b_main(Sort_status *arg)
{
   int ev_done, ts_slot, len, startup, *err = diagnostics.reorder_error;
   int i,  grifc,  type,  rd_avail,  ts_stat,  err_format;
   unsigned int usecs=1000, *evptr, *evstart, *bufend;
   volatile Tsbuf *bufptr, *nxtptr, *newptr;
   unsigned long ts, tslo;

   startup = 1;
   memset(err,           0, REORDER_ERRTYPES*sizeof(int));
   for(i=0; i<MAX_GRIFC; i++){ reorder_init[i] = INIT_WAIT; }
   printf("starting reorder input ...\n");
   for(i=0; i<TIMESLOT_EVENTS; i++){
      timeslot_buffer[i].in_use = 0;  timeslot_buffer[i].next = NULL;
   }
   for(i=0; i<TIMESLOT_SLOTS; i++){ tslot[i] = NULL; }
   tsevents_in = tsbufpos = len = ev_done = ts_stat = err_format = 0;
   evstart = NULL; evptr = bankbuf; bufend = bankbuf + BANK_BUFSIZE;
   while(1){
      if( tsevents_in - tsevents_out > 0.75*TIMESLOT_EVENTS ){
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
      if( len > TIMESLOT_EVENTSIZE ){
         ++err[ERR_LENGTH_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      if( reorder_init[grifc] != 0 ){
         ++err[ERR_INIT_IN]; err[ERR_WORDS_IN] += len;
         if( (ts>>28) == 0 ){ reorder_init[grifc] = 0; }
         else             { --reorder_init[grifc]; }
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      // now have full event without any obvious format errors
      if( ts <= output_ts ){
         //fprintf(stdout, "Order ERROR[in:%d]: ts=%lx[slot0x%04x] out_ts=%lx\n", tsevents_in, ts, ts_slot, output_ts);
         ++err[ERR_LATE_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      if( ts > (output_ts + 3000000000) ){
         //fprintf(stdout, "LARGE GAP[in:%d]: ts=%lx[slot0x%04x] out_ts=%lx [Gap=%ds] DROPPING\n", tsevents_in, ts, ts_slot, output_ts, (ts - output_ts)/100000000 );
         ++err[ERR_EARLY_IN]; err[ERR_WORDS_IN] += len;
         evstart = NULL; err_format = len = ev_done = ts_stat = 0; continue;
      }
      bufptr = NULL;  ts_slot = ts % TIMESLOT_SLOTS;
      // check if any buffer slots available to store this event
      // events are removed partially randomly - need to search
      i=0; do {
         ++i; newptr = timeslot_buffer + tsbufpos++;
         if( tsbufpos >= TIMESLOT_EVENTS ){ tsbufpos = 0; } // wrap
         if( i >= TIMESLOT_EVENTS ){ // ALL BUSY
            i=0; usleep(100*usecs);
            printf("################### REORDER IN_WAIT\n");
            //if( arg->shutdown_midas != 0 ){  break; }
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
            if( ts < nxtptr->ts || nxtptr->in_use == 0 ){ // insert here
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

// reorderbuf: event[evlen] ts, inuse, *next
//    there are TIMESLOT_EVENT total slots
//    there is also a TIMESLOT_SLOTS length linked list
//        of currently stored events
//
// events are stored in reorder buffer above
//    - make sure enough events are in it 8192 (do not read to empty)
//
// count up ts from zero
//    if slot[ts] is occupied - iterate over list
//       typically most ts will be in future
//
//       - would like to skip ahead if ALL ts are in future
void reorder_b_out(Sort_status *arg)
{
   int i, j, wr_avail, wrpos, ts_slot, min_slot, *err;
   volatile Tsbuf *buf, *nxt;
   unsigned int usecs=100;
   unsigned long ts, min_ts, start_ts, loop, skip;
   
   usleep(100*usecs);
   err = diagnostics.reorder_error;
   printf("starting reorder output ...\n");
   tsevents_out = ts = output_ts = skip = loop = start_ts = 0;
   eventbuf_rdpos = eventbuf_wrpos = 0; min_slot = -1;
   while(1){
      // if( arg->shutdown_midas != 0 ){ break; }
      i = tsevents_in - tsevents_out;
      if( i < TS_EVENT_BUFFER ){
         if( !arg->reorder_in_done ){ usleep(usecs); continue; }
         else { break; }
      }
      // continue counting up with ts, and find next in-use timeslot
      while(1){
         //if( (++loop % 500000000) == 0 ){ // 100 million - 1 second realtime
         //   printf("REORDER_LOOP %2d SECONDS [SKIPPED %5.1fs]\n",
         //          (int)((loop+1)/100000000), skip/100000000.0
         //      );
         //}
   // if pass through whole buffer without any event-ts < ts
   // can immediately skip aheaad to first ts in buffer
      // reorder_in keeps linked lists in time order
      //  => they are frequently reordered
      //  => this optimization, storing the earliest timestamp is NOT valid
      //   unless events currently being added by reorder_in are far enough
      //   in the future, as to not insert an earlier event than our "min_ts"
      //
      // Also require same condition to avoid having to lock tslot[ts_slot]
      //  - if events are far future - any update of tslot[ts_slot], after
      //    we stored it in buf just above, will be to add a far future event
      //    
         if( ts - start_ts == TIMESLOT_SLOTS ){
            if( min_slot == -1 ){ printf("REORDER IMPOSSIBLE-ERROR\n");}else{
               skip += min_ts-ts;
               ts = start_ts = min_ts; min_slot = -1;
            }
         }
         ts_slot = ts % TIMESLOT_SLOTS;
         buf = tslot[ts_slot];
         if( buf == NULL ){ ++ts; continue; }
         // if( buf->in_use == 0 ){ printf("NOTINUSE\n"); ++ts; continue; }
         if( buf->ts > ts ){
            if( min_slot == -1 || buf->ts < min_ts ){
               min_slot = ts_slot; min_ts = buf->ts;
            }
            ++ts; continue;
         }
         break;
      }
      min_slot = -1; start_ts = ts; // found event - reset minimum_ts search
      
      // iterate over this in-use linked list tslot[ts % TIMESLOT_SLOTS]
      //   will only be iterating over any equal timestamp events
      while(1){ 
         nxt = buf->next;
         if( buf->ts > ts ){ break; } //  done with this list
         else {                       //  timestamps equal 
            if( buf->ts < ts ){       // (or less than => error)
               //fprintf(stdout, "Order Error[Out:%d] at ts=%ld [buf=%ld]\n", tsevents_out, ts, buf->ts);
               // this is the usual ordering error - where events are delayed
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
            if( wr_avail < buf->evlen ){ break; } // exit sort
            ++tsevents_out;
            //if( (++tsevents_out % 500000) == 0 ||
            //    (ts % 1000000000) == 0 ){ reorder_status(); }
            for(i=0; i<buf->evlen; i++){
               wrpos = eventbuf_wrpos++ % EVENTBUFSIZE;
               event_buffer[wrpos] = buf->event[i];
            }
            pthread_mutex_lock(&nxtlock); // remove from list
            buf->in_use = 0;  buf->next = NULL;
            tslot[ts_slot] = nxt; // update previous nxtptr
            pthread_mutex_unlock(&nxtlock);
         }
         if( (buf=nxt) == NULL ){ break; }
      }
      ++ts; output_ts = ts;
   }
   arg->reorder_out_done = 1;
   printf("Reorder: end output thread [ts:%3ds skip:%3ds] ...\n",
          (int)(ts/100000000), (int)(skip/100000000) );
   reorder_status(0);
   return;
}

void show_status_b()
{
   static int ref[TIMESLOT_EVENTS];
   int i, j, bufidx;
   volatile Tsbuf *buf;
   
   memset(ref, 0, TIMESLOT_EVENTS*sizeof(int));
   printf("TIMESLOT LIST ...\n");
   j=0; for(i=0; i<TIMESLOT_SLOTS; i++){
      if( tslot[i] == NULL ){ ++j; continue; }
      printf(" #%05d[+%4d] -", i, j); buf = tslot[i]; j = 0;
      while(1){
         bufidx =  buf-timeslot_buffer;  ++ref[bufidx];
         printf("  buf#%05d:%s[bufts=%07ld]", bufidx, buf->in_use ? "busy" : "idle", buf->ts);
         if( (buf = buf->next) == NULL){ printf("\n"); break; }
         if( (++j % 2) == 0 ){ printf("\n                 "); }
      }
      j=0;
   }
   printf("EVENT BUFFERS ...\n");
   for(i=0; i<TIMESLOT_EVENTS; i++){
      buf = &timeslot_buffer[i];
      if( buf->in_use == 0 && buf->next == NULL ){ continue; }
      printf("  buf#%05d:%s[bufts%07ld] Nxt=%05d\n", i, buf->in_use ? "busy" : "idle", buf->ts, buf->next==NULL ? 0 : buf->next-timeslot_buffer);
   }
}
