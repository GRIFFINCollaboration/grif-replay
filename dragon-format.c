#include <stdio.h>
#include <time.h>
#include <string.h>

#include "grif-replay.h"
#include "dragon-format.h"
#include "midas-format.h"

#define IO32_BANKLEN 9 // this bank is a fixed length

///////////////////////////////////////////////////////////////////////////
// griffin scheme, midas   -> bankbuf{_wr/rdpos}[4M-Words] RAW BANK DATA
//                 reorder -> reorder 65k-events{ts,next,active,len,data}
//                                    1M-slot-ptrs  @   1us per slot
//                             wait till > 32k events before taking events out
//                         -> event_buffer{_wr/rdpos}[8*1M-words]
//                  unpack -> grif_event[1k-Events] circ buffer for coincs
// 
// *Why only 65k events* ??
///////////////////////////////////////////////////////////////////////////
// dragon/iris, midas   -> unpack -> bankbuf{_wr/rdpos}[4M-Words] UNPACKED DATA
//              reorder -> reorder 65k-events{ts,next,active,len,data}
//                               1M-slot-ptrs  @   1us per slot
//                         wait till > 32k events before taking events out
//                   
//                         -> event_buffer{_wr/rdpos}[8*1M-words]
//
///////////////////////////////////////////////////////////////////////////

// some data in this bank is redundant - check unneeded values, then drop
int unpack_io32_bank(Dragon_event *evt, int *ptr, int len, int type)
{
   Io32_event *io_ptr = &evt->io32;
   int i, j, k;
   if( evt->type == 0 ){ evt->type = type; }
   else if( evt->type != type ){
      printf("Mixture of head.tail banks in single event\n"); return(-1);
   }
   if(len != IO32_BANKLEN){
      printf("IO32 - Wrong bank Length [%d]\n", len); return(-1);
   }
   io_ptr->fw_version    = ptr[0];
   io_ptr->trig_count    = ptr[1]; // since BOR
   io_ptr->trig_ts       = ptr[2];
   io_ptr->readout_begin = ptr[3];
   io_ptr->readout_end   = ptr[4];
   io_ptr->trig_bitmask  = ptr[8];
   if( ptr[3] - ptr[2] != ptr[5] ){
      printf("IO32 - inconsistent Trigger Latency [%d-%d != %d]\n",
             ptr[3], ptr[2], ptr[5] ); return(-1);
   }
   if( ptr[4] - ptr[3] != ptr[6] ){
      printf("IO32 - inconsistent Readout Time [%d-%d != %d]\n",
             ptr[4], ptr[3], ptr[6] ); return(-1);
   }
   if( ptr[4] - ptr[2] != ptr[7] ){
      printf("IO32 - inconsistent Busy Time [%d-%d != %d]\n",
            ptr[4], ptr[2], ptr[7]  ); return(-1);
   }
   return(0);
}
       
int unpack_io32_fifobank(Dragon_event *evt, int *ptr, int len, int type)
{
   Io32_fifo *io_ptr = &evt->io32_ts;
   static unsigned long max_ts;
   unsigned long tshi;
   int i;
   if( evt->type == 0 ){ evt->type = type; }
   else if( evt->type != type ){
      printf("Mixture of head.tail banks in single event\n"); return(-1);
   }
   io_ptr->fw_version = ptr[0];
   io_ptr->timestamp  = ptr[1]; // since BOR
   io_ptr->ts_route   = ptr[2];
   io_ptr->ts_ctrl    = ptr[3];
   io_ptr->rollover   = ptr[4];
   for(i=0; i<len-5; i++){ // could just memcpy - code not much shorter
      if( evt->ts == 0 && ((ptr[i+5] >> 30) & 0x3) == 0 ){ // 1st timestamp-chan-word
         // there is some overlap between upper and lower bits
         // need to examine long datafile (so far rollover always been zero)
         // tshi = ((ts_ctrl >> 15) & 7F) | (rollover << 8);
         //tshi = ((io_ptr->ts_ctrl >> 15) & 0xFF) >> 3; // seems 3 bits at start
         tshi = ((io_ptr->ts_ctrl >> 18) & 0xFF);
         evt->ts =    (tshi << 30) | (ptr[i+5] & 0x3FFFFFFF);
         if( evt->ts < 0 ){
            printf("YY\n");
         }
         max_ts = (evt->ts > max_ts) ? evt->ts : max_ts;
         if( max_ts - evt->ts > 10000000000 ){
            printf("XX\n");
         }
      }
      if( i >= MAX_IO32_FIFOWORDS ){
         printf("IO32_FIFO event too long [%d] - truncated\n", len); break;
      }
      ++io_ptr->nwords;  io_ptr->data[i] = ptr[i+5];
   }
   if( i != len-5 ){ return(-1); } // bank contains errors, unpacking was aborted
   return(0);
}

#define V1190_EVHDR 8
#define V1190_BKHDR 1
#define V1190_EVTRL 16
#define V1190_BKTRL 3
#define V1190_TTAG  17
#define V1190_DATA  0

int unpack_v1190_bank(Dragon_event *evt, int *ptr, int len, int type)
{
   int i, dtype, dcnt, tdc, evid, id, tmp, words, everr;
   V1190_event *ev_ptr = &evt->v1190;
   if( evt->type == 0 ){ evt->type = type; }
   else if( evt->type != type ){
      printf("Mixture of head.tail banks in single event\n"); return(-1);
   }
   everr=0;
   for(i=0; i<len; i++){
      if( everr ){ break; }
      dtype = (ptr[i] >> 27) & 0x1F;
      if( i == 0 && dtype != V1190_EVHDR ){
         printf("V1190 missing Event-Header[dtype=%d]\n", dtype); break;
      }
      if( i != 0 && dtype == V1190_EVHDR ){
         printf("V1190 Event-Header in wrong position[word %d]\n", i); break;
      }
      if( i == len-1 && dtype != V1190_EVTRL ){
         printf("V1190 missing Event-Trailer[dtype=%d]\n", dtype); break;
      }
      if( i != len-1 && dtype == V1190_EVTRL ){
         printf("V1190 Event-Trailer in wrong position[word %d]\n", i); break;
      }
      if( i != len-2 && dtype == V1190_TTAG ){
         printf("V1190 Time-Tag in wrong position[word %d]\n", i); break;
      }
      if( i == 1 && dtype != V1190_BKHDR ){
         printf("V1190 missing Bank-Header[dtype=%d]\n", dtype); break;
      }
      switch( dtype ){
      case V1190_EVHDR: evid = ev_ptr->ev_count = (ptr[i] >> 5) & 0x3FFFF;  break;
      case V1190_BKHDR: tdc = (ptr[i] >> 24) & 0x3;  id = (ptr[i] >> 12) & 0xFFF;
                        ev_ptr->bunch_id = ptr[i] & 0xFFF; dcnt = 0;
                        if( id != (evid & 0xFFF) ){ everr=1;
                           printf("V1190 event-id mismatch [%d vs %d]\n", id, evid);
                        }
                        if( tdc != 0 && tdc != 1 ){
                           everr=1; printf("V1190 unknown tdc bank[%d]\n", tdc);
                        } break;
      case V1190_DATA:  // move LT/CHAN to MSByte for simpler decoding later
                        if( ev_ptr->d_count >= MAX_V1190_WORDS ){
                           ++dcnt; break; // discard, but keep count
                        }
                        tmp = ((ptr[i] & 0x07F80000) << 5) | (ptr[i] & 0x7FFFF);
                        ++dcnt; ev_ptr->data[ev_ptr->d_count++] = tmp; break;
      case V1190_BKTRL: tmp = (ptr[i] >> 24) & 0x3;  id = (ptr[i] >> 12) & 0xFFF;
                        words = ptr[i] & 0xFFF;
                        if( id != (evid & 0xFFF) ){ everr=1;
                           printf("V1190 event-id mismatch [%d vs %d]\n", id, evid);
                        }
                        if( words != dcnt+2 ){ everr=1;
                           printf("V1190 bad word-count[%d vs %d]\n",words,dcnt+2);
                        }
                        if( tmp != tdc ){ everr=1;
                           printf("V1190 tdc bank mismatch[%d bs %d]\n", tmp, tdc);
                        } break;
      case V1190_TTAG:  ev_ptr->time_tag = ptr[i] & 0x7FFFFFF; break;
      case V1190_EVTRL: ev_ptr->status = (ptr[i] >> 24) & 0x7;
                        words          = (ptr[i] >>  5) & 0xFFFF;
                        if( words != len ){ everr=1;
                           printf("V1190 bad word-countB[%d vs %d]\n",words, len);
                        } break;
      default: everr=1; printf("V1190 unknown word type[%d]\n", dtype); break;
      }
   }
   if( dcnt >= MAX_V1190_WORDS ){
      printf("V1190[%s]:Too many words[%d] in event %d\n", type==1?"H":"T", dcnt, ev_ptr->ev_count);
   }
   if( i != len ){ return(-1); } // bank contains errors, unpacking was aborted
   return(0);
}

#define V792_HDR   2
#define V792_DATA  0
#define V792_TRL   4

int unpack_v792_bank(Dragon_event *evt, int *ptr, int len, int id, int type)
{
   int i, dtype, evid, tmp, words, everr;
   V792_event *ev_ptr;
   if( evt->type == 0 ){ evt->type = type; }
   else if( evt->type != type ){
      printf("Mixture of head.tail banks in single event\n"); return(-1);
   }
   if( id == 0 ){
      ev_ptr = &evt->v792a;
   } else if(id == 1 ){
      ev_ptr = &evt->v792b;
   } else {
      printf("V792: unknown module-id [%d]\n", id); return(-1);
   }
   everr=0;
   for(i=0; i<len; i++){
      if( everr ){ break; }
      dtype = (ptr[i] >> 24) & 0x7;
      if( i == 0 && dtype != V792_HDR ){
         printf("V792 missing Header[dtype=%d]\n", dtype); break;
      }
      if( i != 0 && dtype == V792_HDR ){
         printf("V792 Header in wrong position[word %d]\n", i); break;
      }
      if( i == len-1 && dtype != V792_TRL ){
         printf("V792 missing Trailer[dtype=%d]\n", dtype); break;
      }
      if( i != len-1 && dtype == V792_TRL ){
         printf("V792 Trailer in wrong position[word %d]\n", i); break;
      }
      switch( dtype ){
      case V792_HDR:  words = (ptr[i] >> 8) & 0x3F;
                      if( words != len-2 ){ everr=1;
                         printf("V792 bad word-countB[%d vs %d]\n", words, len);
                      } break;                      break;
      case V792_DATA: if( ev_ptr->d_count >= MAX_V792_WORDS ){
                         printf("V792 too many words - ignoring\n"); break;
                      }
                      ev_ptr->data[ev_ptr->d_count++] = ptr[i]; break;
      case V792_TRL:  ev_ptr->ev_count = ptr[i] & 0xFFFFFF; break;
      default: everr=1; printf("792 unknown word type[%d]\n", dtype); break;
      }
   }
   if( i != len ){ return(-1); } // bank contains errors, unpacking was aborted
   return(0);   
}
