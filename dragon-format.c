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
            printf("@");
         }
         max_ts = (evt->ts > max_ts) ? evt->ts : max_ts;
         if( max_ts - evt->ts > 10000000000 ){
            printf("$");
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

// Caen V1190 TDC data format (from caen manual) ...
//   each data word starts with a 5bit "word-type", followed by type-dependant fields
//                        TYPE
//           GlobalHdr   0100_0 22bit-EvCnt 5bit-Geo
//             [TdcHdr]  0000_1 X 2bit-TDC 12bit-Event-id 12bit-Bunch-id
//                Data   0000_0 Lead/Trl 7bit-Chan 19bit-Measurement
//            [TdcTrlr]  0001_1 X 2bit-TDC 12bit-Event-id 12bit-Word-Cnt
//             [ERROR]   0010_0 X 2bit-TDC 9*X 15bit-Error-Flags
//  *ALSO* ExtTrigTime   1000_1  27bit ext-trigger-time
//           [trig arrival rel to cnt_reset [lsb 800ns]]
//        GlobalFooter   1000_0  3bit-STAT 3*X 16bit-WrdCnd 5bit-Geo
//
// [ => dataval is [18:0]  datachan# is 7bits [26:20] ** MOVED BY US to [30:24] ** ]
//
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
      case V1190_DATA:  if( ev_ptr->d_count >= MAX_V1190_WORDS ){
                           ++dcnt; break; // discard, but keep count
                        }
                        // move LT/CHAN to MSByte for simpler decoding later
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

// Caen V79s ADC/QDC data format (from caen manual) ...
//      ADC  Hdr   5bit-Geo[Unused] 0 1 0 8bit-Crate[unused] 0 0 6bit WrdCnt 8*X
//      QDC  Data  5bit-Geo[Unused] 0 0 0 3*X 5bit chan XX Undr Ovr 12bit-ADC
//           ...
//         Footer  5bit-Geo[Unused] 1 0 0 EvtCnt[24bits]
//
// [ => dataval is [11:0]  datachan# is 5bits [20:16]  ]
//
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

// TODO - add single items [eg dssd_tdc_front, odd tacs/tdcs
//   ** ALSO ** multiple v792 modules in tail [map 2nd modules to channels 16-31]

int unpack_head_event(Dragon_event *evt)
{
   Head_data *head = &evt->head_tail_data.head_data;
   int i, cnt, val, edge, badval, srcchan, dstchan;

   // v792a adc data [single module only for head events]
   if( (cnt = evt->v792a.d_count) > 0 ){
      for(i=0; i<cnt; i++){
         val = evt->v792a.data[i];
         srcchan = (val >> 16) & 0x1f;
         badval = (val >> 12) & 0x3; // under/overflow bits
         val &= 0xfff; 
         if( (dstchan = head_adc_dstchan[srcchan]) == -1 ){ continue; } // unassigned
         switch( head_adc_dettype[srcchan] ){
         case SUBSYS_BGO:
            if( dstchan >= 0 && dstchan < BGO_MAXCHAN ){ 
               head->bgo_energy[dstchan] = badval ? 0 :
                  bgo_adc_offset[dstchan] + bgo_adc_slope[dstchan] * val;
            } else { printf("UnpckHead:CodingErr1\n"); } break;
         case SUBSYS_SB:
            if( dstchan >= 0 && dstchan < SB_MAXCHAN ){ 
               head->sb_energy[dstchan] = badval ? 0 :
                  sb_adc_offset[dstchan] + sb_adc_slope[dstchan] * val;
            } else {
               printf("UnpckHead:CodingErr2\n"); } break;
         case SUBSYS_NAI:
            if( dstchan >= 0 && dstchan < NAI_MAXCHAN ){ 
               head->nai_energy[dstchan] = badval ? 0 :
                  nai_adc_offset[dstchan] + nai_adc_slope[dstchan] * val;
            } else { printf("UnpckHead:CodingErr3\n"); } break;
         case SUBSYS_GE:
            if( dstchan == 0 ){ // only a single channel
               head->ge_energy = badval ? 0 :ge_adc_offset + ge_adc_slope * val;
            } else { printf("UnpckHead:CodingErr4\n"); } break;
         default: printf("UnpckHead:CodingErr5\n"); break;
         }
      }
   }
   // v1190 tdc data words ...
   if( (cnt = evt->v1190.d_count) > 0 ){
      for(i=0; i<cnt; i++){
         val = evt->v1190.data[i];
         srcchan = (val >> 24) & 0x7f; // NOTE: chan was moved to MSB from original position
         edge    = (val >> 31) & 0x1;
         val = (val & 0xffffff) | (edge<<31);
         if( (dstchan = head_tdc_dstchan[srcchan]) == -1 ){ continue; } // unassigned
         switch( head_tdc_dettype[srcchan] ){
         case SUBSYS_BGO:
            if( dstchan >= 0 && dstchan < BGO_MAXCHAN ){ 
               head->bgo_time[dstchan] = bgo_tdc_offset[dstchan] +
                  bgo_tdc_slope[dstchan] * val;
            } else { printf("UnpckHead:CodingErr6\n"); } break;
         case SUBSYS_XTDC:
            if( dstchan == 0 ){ head->xtdc = head_xtdc_offset + head_xtdc_slope * val; }
            else { printf("UnpckHead:CodingErr7\n"); } break;
         case SUBSYS_RFTDC:
            if( dstchan == 0 ){ head->rftdc = head_rftdc_offset + head_rftdc_slope * val; }
            else { printf("UnpckHead:CodingErr8\n"); } break;
         case SUBSYS_TDC0:
            if( dstchan == 0 ){ head->tdc0 = head_tdc0_offset + head_tdc0_slope * val; }
            else { printf("UnpckHead:CodingErr9\n"); } break;
         default: printf("UnpckHead:CodingErr10\n"); break;
         }
      }
   }
   return(0);
}

int unpack_tail_event(Dragon_event *evt)
{
   Tail_data *tail = &evt->head_tail_data.tail_data;
   int i, j, cnt, val, edge, badval, srcchan, dstchan;

   for(j=0; j<2; j++){  // v792a adc data (two modules)
      if( (cnt = j ? evt->v792a.d_count : evt->v792b.d_count) <= 0 ){ continue; }
      for(i=0; i<cnt; i++){
         val = j ? evt->v792a.data[i] : evt->v792b.data[i];
         srcchan = (V792_MAXCHAN * j) + (val >> 16) & 0x1f;
         badval = (val >> 12) & 0x3; // under/overflow bits
         val &= 0xfff; 
         if( (dstchan = tail_adc_dstchan[srcchan]) == -1 ){ continue; } // unassigned
         switch( tail_adc_dettype[srcchan] ){
         case SUBSYS_DSSD:
            if( dstchan >= 0 && dstchan < DSSD_MAXCHAN ){ 
               tail->dssd_energy[dstchan] = badval ? 0 :
                  dssd_adc_offset[dstchan] + dssd_adc_slope[dstchan] * val;
            } else { printf("UnpckTail:CodingErr1\n"); } break;
         case SUBSYS_IC:
            if( dstchan >= 0 && dstchan < IC_MAXCHAN ){ 
               tail->ic_energy[dstchan] = badval ? 0 :
                  ic_adc_offset[dstchan] + ic_adc_slope[dstchan] * val;
            } else { printf("UnpckTail:CodingErr2\n"); } break;
         case SUBSYS_MCP:
            if( dstchan >= 0 && dstchan < MCP_MAXCHAN ){ 
               tail->mcp_energy[dstchan] = badval ? 0 :
                  mcp_adc_offset[dstchan] + mcp_adc_slope[dstchan] * val;
            } else { printf("UnpckTail:CodingErr3\n"); } break;
         case SUBSYS_MCPTAC:
            if( dstchan == 0 ){ 
               tail->mcptac_energy = badval ? 0 : mcptac_adc_offset + mcptac_adc_slope * val;
            } else { printf("UnpckTail:CodingErr4\n"); } break;
         default: printf("UnpckTail:CodingErr5\n"); break;
         }
      }
   }
   // v1190 tdc data words ...
   if( (cnt = evt->v1190.d_count) > 0 ){
      for(i=0; i<cnt; i++){
         val = evt->v1190.data[i];
         srcchan = (val >> 24) & 0x7f; // NOTE: chan was moved to MSB from original position
         edge    = (val >> 31) & 0x1;
         val = (val & 0xffffff) | (edge<<31);
         if( (dstchan = tail_tdc_dstchan[srcchan]) == -1 ){ continue; } // unassigned
         switch( tail_tdc_dettype[srcchan] ){
         case SUBSYS_DSSD:
            if( dstchan >= 0 && dstchan < DSSD_TDCCHAN ){ 
               tail->dssd_time[dstchan] = dssd_tdc_offset[dstchan] +
                  dssd_tdc_slope[dstchan] * val;
            } else { printf("UnpckTail:CodingErr6\n"); } break;
         case SUBSYS_IC:
            if( dstchan >= 0 && dstchan < IC_MAXCHAN ){ 
               tail->ic_time[dstchan] = ic_tdc_offset[dstchan] +
                  ic_tdc_slope[dstchan] * val;
            } else { printf("UnpckTail:CodingErr7\n"); } break;
         case SUBSYS_MCP:
            if( dstchan >= 0 && dstchan < MCP_TDCCHAN ){ 
               tail->mcp_time[dstchan] = mcp_tdc_offset[dstchan] +
                  mcp_tdc_slope[dstchan] * val;
            } else { printf("UnpckTail:CodingErr8\n"); } break;
         case SUBSYS_XTDC:
            if( dstchan == 0 ){ tail->xtdc = tail_xtdc_offset + tail_xtdc_slope * val; }
            else { printf("UnpckTail:CodingErr9\n"); } break;
         case SUBSYS_RFTDC:
            if( dstchan == 0 ){ tail->rftdc = tail_rftdc_offset + tail_rftdc_slope * val; }
            else { printf("UnpckTail:CodingErr10\n"); } break;
         case SUBSYS_TDC0:
            if( dstchan == 0 ){ tail->tdc0 = tail_tdc0_offset + tail_tdc0_slope * val; }
            else { printf("UnpckTail:CodingErr11\n"); } break;
         default: printf("UnpckTail:CodingErr12\n"); break;
         }
      }
   }
   return(0);
}
