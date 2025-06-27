#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "grif-replay.h"
#include "grif-format.h"
#include "midas-format.h"

void dbg_dump_event(unsigned *buf, int len){
   int i; for(i=0; i<len; i++){ printf("0x%08x ", *buf++); } printf("\n");
}
void dbg_grifbuf(unsigned *evstrt, int reorder){
   if( reorder ){
   printf("EVBUF[size:%d]: wrpos:%8ld rdpos:%8ld [%8ld avail]\n",
      EVENTBUFSIZE, eventbuf_wrpos%EVENTBUFSIZE, eventbuf_rdpos%EVENTBUFSIZE,
      eventbuf_wrpos-eventbuf_rdpos );
   printf("                                  evstart:%8ld\n", evstrt-event_buffer);
   } else {
   printf("BANKBUF[size:%d]: wrpos:%8ld rdpos:%8ld [%8ld avail]\n",
      BANK_BUFSIZE, bankbuf_wrpos%BANK_BUFSIZE, bankbuf_rdpos%BANK_BUFSIZE,
      bankbuf_wrpos-bankbuf_rdpos );
   printf("                                  evstart:%8ld\n", evstrt-bankbuf);
   }
}

int debug;
Grif_event grif_event[PTR_BUFSIZE];
short waveform[MAX_SAMPLE_LEN];
int     scalar[MAX_SCALAR_LEN];
int waveforms = 1; // process [1] or ignore them

unsigned tmp_buf[256];
volatile long grif_evcount;
volatile long grifevent_wrpos;
extern volatile long grifevent_rdpos;
static int unpack_grif3_event(unsigned *buf, int len, Grif_event *, int);
extern Grif_event grif_event[PTR_BUFSIZE];
// for a single coincidence window - can mark start of window wrt latest frag
// but for multiple windows, may be better to look backwards, as far as req.
//    (subject to buffer size)
extern volatile long grifevent_rdpos, grifevent_wrpos;

#define GRIF_ERR_TYPES      5
#define GRIF_ERR_HDR        0
#define GRIF_ERR_TRLR       1
#define GRIF_ERR_ADDR       2
#define GRIF_ERR_PHWORDS    3
#define GRIF_ERR_TRIGMATCH  4

static int grif_err[GRIF_ERR_TYPES];
void grif_status(int current_time)
{
   static int prev_time, prev_evt;
   int dt = current_time - prev_time, de = grif_evcount - prev_evt;
   int sum = grif_err[GRIF_ERR_ADDR] + grif_err[GRIF_ERR_HDR] + grif_err[GRIF_ERR_TRLR];

   prev_time = current_time; if( dt == 0 ){ dt = 1; }
   prev_evt  = grif_evcount;

   printf("GrifEvt: in:%10ld err:%10d[%5.1f%%]         %6.3f Mevt/s\n         ",
          grif_evcount, sum, (100.0*sum)/grif_evcount, de/(1000000.0*dt) );
   printf("[Hdr:%d Trlr:%d addr:%d phwrds:%d trigMatch%d]\n",
          grif_err[GRIF_ERR_HDR], grif_err[GRIF_ERR_TRLR], grif_err[GRIF_ERR_ADDR],
                           grif_err[GRIF_ERR_PHWORDS], grif_err[GRIF_ERR_TRIGMATCH] );
}

void grif_main(Sort_status *arg)
{
   int i, wrpos, len, rd_avail, wr_avail, wcnt, wrap;
   unsigned bufsize, *bufstart, *bufend, *evstart, *evptr;
   unsigned int usecs=100;
   Grif_event *ptr;

   memset(grif_err, 0, GRIF_ERR_TYPES*sizeof(int));
   bufsize = EVENTBUFSIZE;
   evptr = evstart = bufstart = event_buffer;
   bufend = event_buffer + bufsize;
   grif_evcount = grifevent_wrpos = wrpos = wcnt = wrap = 0;
   printf("starting grif thread ...\n");
   while(1){ ptr = &grif_event[wrpos];
      // if( arg->shutdown_midas != 0 ){  break; }
      rd_avail =eventbuf_wrpos - eventbuf_rdpos;
      if( rd_avail < 10 ){  // leave some margin
         if( arg->reorder_out_done ){ break; }
          usleep(usecs); continue;
      }
      if( (( (*evptr)>>28 )&0xf) != 0xE ){
	if( (++wcnt % 256) == 0 ){
            printf("SEGFAULT:wcnt=%d\n", wcnt);
         }
         ++evptr;       // DO NOT consume input data here or event
         if( evptr >= bufend){  // start can be overwritten before we use it
            evptr -= bufsize;
            if( (wrap = wcnt) > 255 ){
               printf("SEGFAULT:wcnt=%d\n", wcnt);
            }
            memcpy(tmp_buf, (char *)evstart, wcnt*sizeof(int));
         }
         continue;
      }
      // *evptr == 0xE......
      if( wrap ){
         memcpy(&tmp_buf[wrap+1], bufstart, (wcnt+1-wrap)*sizeof(int));
         evstart = tmp_buf; wrap = 0;
      }
      len = wcnt+1; // include word#0
      while(1){ // wait for space to store event ..
         rd_avail = grifevent_wrpos - grifevent_rdpos;
         if( (wr_avail = PTR_BUFSIZE - rd_avail) != 0 ){ break; }
         usleep(usecs);
      }
      if( unpack_grif3_event(evstart, len, ptr, waveforms) == 0 ){
         //if( !arg->sort_thread ){ process_event(ptr, wrpos); }
         //if( ( ++grif_evcount % 1000000) == 0 ){
         //   fprintf(stderr, "%10d events ...\n", grif_evcount);
         //}
      } else { --grifevent_wrpos; } // dump this event
      if( ++evptr >= bufend ){ evptr -= bufsize; }
      ++wcnt; evstart = evptr;  ++grifevent_wrpos;
      eventbuf_rdpos += wcnt; // consume input data here
      wcnt = 0;  wrpos = grifevent_wrpos %  PTR_BUFSIZE;  ++grif_evcount;
   }
   arg->grif_sort_done = 1;
   printf("grif_ordered thread finished rd:%ld wr:%ld\n",
          eventbuf_rdpos, eventbuf_wrpos );
   return;
}

int process_grif3_bank(unsigned *evntbuf, int length)
{
   unsigned *bufend = evntbuf+length, *evstrt=evntbuf;
   unsigned *ptr = evntbuf;
   Grif_event *evt;
   int wrpos;

   while( ptr < bufend ){
      wrpos = grifevent_wrpos % PTR_BUFSIZE;
      evt = &grif_event[wrpos];
      if( ((*(ptr++)>>28)&0xff) == 0xE ){
         if( unpack_grif3_event(evstrt, ptr-evstrt, evt, 0) == 0 ){
            process_event(evt, wrpos);
         }
         evstrt=ptr; ++grifevent_wrpos;
      }
   }
   return(0);
}

extern short address_chan[MAX_ADDRESS];
extern int reorder_events_read;
// header different - no multi qt's, integral bits changed
int unpack_grif3_event(unsigned *evntbuf, int evlen, Grif_event *ptr, int process_waveforms)
{
   int i, ebad, type, value, qtcount, master_port, grifc_port, done=0;
   unsigned int val32, *evstrt = evntbuf;
   static int savelen, prevtrig, errcount;
   int *wave_ptr = NULL;  int discard = 0;

   if( debug ){ printf("--CLR EVT[%4ld]\n", ptr - grif_event ); }
   memset(ptr, 0, sizeof(Grif_event) );
   ptr->master_id = -1;  // ptr->file_id = grif_evcount;
   if( ((*evntbuf) & 0x80000000) != 0x80000000 ){
      ++grif_err[GRIF_ERR_HDR];
      //fprintf(stderr,"griffin_decode: bad header in event %d\n", grif_evcount );
      //dump_event(evstrt,evntbuflen); return(-1);
   }
   for(i=0; i<evlen; i++){
      val32 = *(evntbuf++);
      type = val32 >> 28; value = val32 & 0x0fffffff;
      switch( type ){
      case 0x8:                                            /*  Event header */
         if( i != 0 ){
            ++grif_err[GRIF_ERR_HDR];
            //fprintf(stderr,"Event %d(chan%d) 0x8 not at start of data\n",
            //   grif_evcount, ptr->chan );
	 }
	 qtcount = 0;
         ptr->dtype  = ((value & 0x000000F) >>  0);

  //       if( ptr->dtype == 10 || ptr->dtype==11 ){
	 //   printf("RCS\n");
	// }
         ptr->address= ((value & 0xFFFF0) >>  4);
	 // ptr->address >= 0x8000 - currently this will be caen data events
	 //   which have had their address altered to allow reordering
         //   now the address should be changed back to what it was
         if( ptr->address == 0xFFFF ){
            val32 -= 0x80000000;
         }
         if( (unsigned)(ptr->address) >= 0x8000 ){
	    extern int grifc_to_boardid[16];
	    int board_id, grifc;
            grifc = (ptr->address >> 12 ) & 0xF;
	    if( (board_id = grifc_to_boardid[grifc]) != -1 ){
	       ptr->address = 0x8000 + (board_id * 0x100) + (ptr->address & 0xFF);
            }
	 }
         // fprintf(stdout,"%d\n",ptr->address);
         ptr->chan = address_chan[(unsigned short)ptr->address];
         if( ptr->dtype != 0xF && (ptr->chan < 0) && ptr->address != 0xFFFF ){
            ++grif_err[GRIF_ERR_ADDR];
            //if( ++errcount < 100 || (errcount % 1000 == 0) ){
             fprintf(stderr,"Ignoring Event - Unknown address [0x%04x] returns chan %d\n", ptr->address, ptr->chan);
            //}
            return(-1);
         }
         ptr->alt_chan = -1; // initialize as -1, used in addback
         //wave_ptr  = &ptr->waveform_length;     /* should be zero here */
 	 break;
      case 0x9:                      /* Channel Trigger Counter [AND MODE] */
         ptr->trig_req =  value & 0x0fffffff;
         break;
      case 0xa:                                           /*  Time Stamp Lo */
         ptr->ts   = ( value & 0x0fffffff );
         break;
      case 0xb:                               /* Time Stamp Hi and deadtime */
	 ptr->ts   |= ( (long)(value & 0x0003fff) << 28);
	 ptr->deadtime     = ( (value & 0xfffc000) >> 14);
         ptr->ts = ptr->ts;
	 break;
      case 0xc:                                             /* waveform data */
         if( wave_ptr == NULL){
          //  fprintf(stderr,"griffin_decode: no memory for waveform\n");
         } else if( process_waveforms == 1 ){/* + 14->16bit sign extension */
	    waveform[(*wave_ptr)  ]   = value & 0x3fff;
            waveform[(*wave_ptr)++] |= ((value>>13) & 1) ? 0xC000 : 0;
	    waveform[(*wave_ptr)  ]   =(value & 0xfffc000) >> 14;
            waveform[(*wave_ptr)++] |= ((value>>27) & 1) ? 0xC000 : 0;
	 }
	 break;
      case 0xd:                                   /* network packet counter */
         ptr->net_id = val32;
	 break;
      case 0xe:      // 14bit acc 14bit req                /* Event Trailer */
         if( (i+1) != evlen ){
            ++grif_err[GRIF_ERR_TRLR];
            //fprintf(stderr,"Event %d(chan%d) 0xE before End of data\n",
            //   grif_evcount, ptr->chan );
	 }
	 ptr->trig_acc = (val32 & 0xfffc000) >> 14;
         if( ptr->dtype < 12 && (val32 & 0x3fff) != (ptr->trig_req & 0x3fff)
                                                    && val32 != 0xefffffff ){
            ++grif_err[GRIF_ERR_TRIGMATCH];
            //fprintf(stderr,"Event 0x%x(chan%d) trig_req mismatch [%d]!=[%d]\n",
            // grif_evcount, ptr->chan, (ptr->trig_req&0x3fff), (val32&0x3fff) );
            // dbg_dump_event(evstrt, evlen);
         }
         break;
      case 0x0: case 0x1: case 0x2: case 0x3:
      case 0x4: case 0x5: case 0x6: case 0x7:
         if( ptr->dtype == 0xF ){ // scalar events have different format
            if( ptr->address != 0xFFFF ){ discard = 1; }
            //ptr->scl_present = 1;
            //scalar[ptr->scalar_length++]  = val32;
            break;
         }
	 if( i < 3 ){ // header words - either [mstpat/ppg mstid] or ppg
            if( ptr->address == 0xFFFF ){
               ptr->master_pattern = val32;
            } else {
 	       //ptr->wf_present = (val32 & 0x8000) >> 15;
	       ptr->pileup     = (val32 & 0x001F);
	       if( ( *(evntbuf) >> 31 ) == 0 ){
                  ptr->master_id   = *(evntbuf++);
                  ++i;
               }
            }
	    break;
	 } else { // if dtype=6, maybe RF - extend sign from 30 to 32bits
	    if( ptr->dtype == 6 && (val32 & (1<<29)) ){ val32 |= 0xC0000000; }
            switch(++qtcount){
            case 1:  /* Energy */
            ptr->q1  = (ptr->dtype==6) ? val32 : val32 & 0x01ffffff;
            ebad = (value >> 25) & 0x1;
            ptr-> integ1 |= ((val32 & 0x7c000000) >> 17); ptr->nhit = 1;
            break;
            case 2: /* CFD Time */
	       if(ptr->dtype==6){
	          ptr->cfd      = val32 & 0x3ff;
	          ptr->cc_short = ((val32>>10) & 0x7fff);
	       } else {
	          ptr->cfd = val32 & 0x003fffff;
	       }
               ptr-> integ1 |= ((val32 & 0x7FC00000) >> 22);
               break;
            case 3:  /* descant long*/
               if(ptr->dtype==6){
                 //ptr->cc_long  = val32; // Not used? cc_long was changed to psd in grif-format.h
               } else { ptr->integ2 =   val32 & 0x003FFF;
                 ptr->nhit   = ((val32 & 0xFF0000) >> 16);
               }
               break;
            case 4:  /* descant short*/
               if(ptr->dtype==6){ ptr->cc_short  = val32; }
               else { ptr->q2 =  val32 & 0x3FFFFF;
                 ebad  = (val32 >> 25) & 0x1;
               }
               break;
            case 5:
               ptr->integ3 =   val32 & 0x00003FFF;
               ptr->integ4 = ((val32 & 0x3FFF0000) >> 16);
               break;
            case 6: ptr->q3 =  val32 & 0x3FFFFF;
               ebad  = (val32 >> 25) & 0x1;
               break;
            case 7: ptr->q4 =  val32 & 0x3FFFFF;
               ebad  = (val32 >> 25) & 0x1;
               break;
            default:
               ++grif_err[GRIF_ERR_PHWORDS];
               //fprintf(stderr,"Event %d(chan%d) excess PH words\n",
               //   grif_evcount, ptr->chan );
               break;
            }
         }
         break;
      case 0xf: fprintf(stderr,"griffin_decode: 0xF.......\n");
                /* Unassigned packet identifier */ return(-1);
      default:  fprintf(stderr,"griffin_decode: default case\n"); return(-1);
      }
   }

   return( discard ); // discard scalars other than ppg-"scalar"
}
