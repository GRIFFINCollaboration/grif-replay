#include <stdio.h>
#include <stdlib.h> // abort()
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "grif-replay.h"
#include "midas-format.h"

static char recbuf[RECORDSIZE];
static int recbufpos, recordlen;

static int errcount;

static int bank_buffer[MAX_BANK_SIZE];
static int bankpos, bank_len;

static Midas_event_header          ev_head;
static Midas_allbank_header   allbank_head;
static Midas_bank_header         bank_head;

static int swap_required;
static int first_bank;

static void swapInt (char* b, int len)
{ char t; while( len > 0 ){ t=b[0]; b[0]=b[3]; b[3]=t; t=b[1]; b[1]=b[2]; b[2]=t; b+=4; len-=4; } }
static void swapShort(char* b, int len)
{ char t; while( len > 0 ){ t=b[0]; b[0]=b[1]; b[1]=t; b+=2; len-=2; } }
static void swapWords (char* b, int len)
{ char t; while( len > 0 ){ t=b[0]; b[0]=b[2]; b[2]=t; t=b[1]; b[1]=b[3]; b[3]=t; b+=4; len-=4; } }

extern Sort_metrics diagnostics;
void midas_status(int current_time)
{
   int dt = diagnostics.midas_last_timestamp - diagnostics.midas_run_start;
   float rate = diagnostics.midas_file_bytes / ((dt == 0) ? 1 : dt);
   printf("MIDAS  : in %10ld runtime:%3ds =>    [%8.4f Mbytes/s]\n",
      diagnostics.midas_file_bytes, dt, rate/1024/1024
   );
   dt = current_time - diagnostics.run_sort_start;
   rate = diagnostics.midas_file_bytes / ((dt == 0) ? 1 : dt);
   printf("%22s sorttime:%3ds => [%8.4f Mbytes/s]\n", "", dt, rate/1024/1024);
}
//////////////////////////////////////////////////////////////////////////
///////       Short section to translate caen data to grif fmt     ///////
//////////////////////////////////////////////////////////////////////////

#define TMP_BANKBUFLEN 65536
#define CAEN_ID            5
#define CAEN_WORDS         8
#define DESCANT_DTYPE      6        // ### Read from ODB - what if missing?
unsigned tmp_bankbuf[TMP_BANKBUFLEN];

static int caen_event_id;

// 0x973e90 <recbuf+953744>: ...
//    0xa0000066           0x60000008       0x00000000 0x0007ceef
//  HDR[0xA],0x66Wrds=102  BrdID=12,mask=8    COUNTER     TIME
//    0x80000062           0x7200_0000      0x0000e5ef   0x00000327
//  HDR[0x8],0x62Wrds    0111_0010 #samp    [MSB],TSTAMP   Extra[CFD]
//              0:Dual-Trace 11:?? 1:extra 0:Wvfm  010:ex_fmt
//    0x00cc00c5 0x00014ae0 0x000003ff 0x00a90001
//     E,CC    [MSB],TSTAMP  Extra[CFD]    E,CC
//    0x00020a17 0x0000026f 0x002f0035 0x00024ceb
//    0x00000187 0x00a800b6 0x00027607 0x000003ff
//    0x01020004 0x8003bfd3 0x0000027f 0x02930133
//    0x0003de8a 0x000001bf 0x00a500aa 0x00057332
//    0x00000233 0x00e700ce 0x00058a8a 0x0000037f
//    0x00270013 0x00058b8c 0x00000357 0x00bf00a9
//    0x00081f71 0x0000008d 0x0090008c 0x0008c7fd
//    0x000001cf 0x009f00a7 0x00000091 0x00c700ba
//    0x800c1e0e 0x0000033f 0x004d0056 0x000e98d3
//    0x000000cb 0x0054005c 0x000f31e6 0x0000031f
//    0x00840099 0x00108028 0x000001ab 0x008f0095
//    0x00109d7f 0x0000036f 0x00760043 0x00124c35
//    0x0000011f 0x0095008f 0x001293f7 0x000003b7
//    0x00f600dd 0x00155500 0x000002af 0x00fc00e1
//    0x001624a6 0x000000bd 0x00c600c5 0x0016e1c1
//    0x0000034f 0x005c0075 0x8019adf9 0x0000000f
//    0x00240023 0x8019db7d 0x0000011f 0x00ad009a
//    0x001a69aa 0x00000337 0x007a0077 0x801c9553
//    0x00000121 0x052e04f8 0x801d8b6c 0x00000073
//    0x004b0047 0x001ddab4 0x000000d5 0x00b400c2
//    0x001ef587 0x0000032b 0x00ae00ba 0x001f3a8f
//    0x000000e3 0x00d200c8 0x00000003 0x000000bf
//    0x6704b5d1 0x000001a8 0x000001a0 0x00000001
//    0x4e454143 0x01980006 0xa0000066 0x60000020
//    0x00000001 0x000cd7a2 0x80000062 0x72000000
//    0x80013d8e 0x00000177 0x00f200e0 0x00072e42
//    0x000000df 0x00280019 0x00078b89 0x000000dd
//    0x00a500b3 0x0008c690 0x0000032f 0x002e002e


// assign caen board_id to grifc mapping in order of channel definition in odb
int board_to_grifc[32];   // board-id is 5bits -> 32 values
int grifc_to_boardid[16]; // to allow translating back to original address
int read_caen_odb_addresses(int odb_daqsize, unsigned short *addr_table)
{
   int i, addr, board_id, curr_grifc=8;
   memset(board_to_grifc,   -1, 32*sizeof(int));
   memset(grifc_to_boardid, -1, 16*sizeof(int));
   for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
      if( ((addr = addr_table[i]) & 0x8000 ) == 0 ){ continue; }
      board_id = (addr - 0x8000) >> 8;
      if( board_to_grifc[board_id] != -1 ){ continue; }
      board_to_grifc[board_id] = curr_grifc;
      grifc_to_boardid[curr_grifc++] = board_id;
      if( curr_grifc >= 16 ){
         printf("midas-format.c:read_caen_odb_addresses - more than 8 caen board id's - can't sort this data\n"); abort();
      }
   }
   return(0);
}
int translate_caen_bank(unsigned *ptr, int len)
{
   int i, j, k, m, outpos, tmp_outpos, out_start, addr, msb, board_id;
   int e, cc, cfd, ovr, lost, samples, dual_trace, ext_format, use_waveform;
   int ev_size, digital, chan_words, chan_pair, chan_mask, board_words;
   unsigned long timestamp, grif_ts;

   i = outpos = 0;
   while( i < len ){ // bank can contain several board structures ?
      if( ( ptr[i] >> 28 ) != 0xA ){ return(-1); }// board header
      board_words = ptr[i  ] & 0xFFFFFFF;
      if( i++ + board_words > len ){ return(-1); }
      board_id  =  ptr[i  ] >>  27;   // 5-bit id
      chan_mask =  ptr[i  ] & 0xFF; // 8-bit mask
      i+=3; // SKIP AGG/AGG-TIME
      for(chan_pair=0; chan_pair<8; chan_pair++){
         if( ! ( (chan_mask >> chan_pair) & 0x1) ){ continue; }
         if( ( ptr[i  ] >> 31 ) != 0x1 ){ return(-1); } // chan header bit
         chan_words = ptr[i  ] & 0x3FFFFF;              // 22-bit size
         if( i++ + chan_words > len ){ return(-1); }
         if( ( ptr[i  ] >> 29 ) != 0x3 ){ return(-1); } // next header bits?
         dual_trace  =  (ptr[i]>>31 ) & 0x1;
         ext_format  = ((ptr[i]>>28 ) & 0x1) ? ((ptr[i] >> 24 )& 0x7)+1 : 0;
         use_waveform = (ptr[i]>>27 ) & 0x1;
         samples   = 4*( ptr[i++] & 0xFFFF );
         ev_size = samples + 2 + (ext_format != 0);
         if( !(ev_size==2 && (chan_words % ev_size != 0)) && // UNCLEAR
                             (chan_words % ev_size != 2) ){ return(-1); }
         // event data for this channel-pair follows ...
         for(j=0; j<(chan_words-2)/ev_size; j++){
         out_start = outpos;
            msb       = ptr[i  ] >> 31; // msb set => odd-channel
            timestamp = ptr[i++] & 0x7FFFFFFF;
	  //addr = (0x8000 + (board_id * 0x100) + 2*chan_pair + msb);
            addr = (0x0000 + (board_to_grifc[board_id] * 0x1000) +
		                               2*chan_pair + msb);

   /* 0 */  tmp_bankbuf[outpos++] = (0x8<<28) + ((CAEN_ID)<<25) +
               ((CAEN_WORDS)<<20) + (addr<<4) + (DESCANT_DTYPE);
   /* 1 */  tmp_bankbuf[outpos++] = (use_waveform<<15) | dual_trace;
   /* 2 */  tmp_bankbuf[outpos++] = caen_event_id++ & 0x7FFFFFFF;
   /* 3 */  tmp_bankbuf[outpos++] = (0x9<<28) + 0;
   /* 4 */  tmp_bankbuf[outpos++] = 0; // fill in timestamp later
   /* 5 */  tmp_bankbuf[outpos++] = 0; // fill in timestamp later

            if( use_waveform ){
               // waveform uses all 32bits of each word - no room for grif hdr
               // store digital waveform after main waveform
               tmp_outpos = outpos + samples;  digital = 0;
               for(k=0; k<samples; k++){
                  tmp_bankbuf[outpos++] = (0xC<<28) + (ptr[i] & 0x3FFF) +
                                                 ((ptr[i] & 0x3FFF0000)>>2);
                  // digitals use 8 of the 14bits in upper and lower samples
                  m = ((ptr[i  ]>>14) & 0x3) + ((ptr[i  ]>>18) & 0xC000);
                  digital |= (m << 2*(k%4));
                  if( k%4 == 3 ){
                     tmp_bankbuf[tmp_outpos++] = (0xC<<28) + digital;
                     digital = 0;
                  }
                  ++i;
               }
               outpos = tmp_outpos;
            }
            e = cfd = cc = ovr = lost = 0;
            switch( ext_format ){
            case 0:  break; // no extra word present
            case 1:  timestamp |= ((long)(ptr[i] & 0xffff0000))<<15;
                    /* 16 bits of baseline */
                     ++i; break;
            case 2:  timestamp |= ((long)(ptr[i] & 0xffff0000))<<15;
                    /* [15..12} - 4 bits of flags */
                     ++i; break;
            case 3:  timestamp |= ((long)(ptr[i] & 0xffff0000))<<15;
                    /* [15..12} - 4 bits of flags */
                     cfd = ptr[i] & 0x3FF;  ++i; break;
            case 4:  ++i; break; // missing from analyzer code?
            case 5: tmp_bankbuf[out_start+3] |= (ptr[i]>>16);
                    lost = ptr[i] & 0xffff; ++i; break;
            case 6:  ++i; break;   // cfd something ???
            case 7:  ++i; break;   // missing from analyzer code?
            case 8:  ++i; break;   // debug word
            default: ++i; break;
            }
            grif_ts = (unsigned long)(ceil(timestamp / 5)); // (timestamp was in 2ns units)
            tmp_bankbuf[out_start+4] = (0xA<<28) +  (grif_ts & 0xFFFFFFF);
            tmp_bankbuf[out_start+5] = (0xB<<28) + ((grif_ts>>28) & 0x3fff);
            e   =  ptr[i  ] >> 16; ovr = (ptr[i] >> 15) & 0x1;
            cc  = (ptr[i++]      ) & 0x7fff; // 15bits
            tmp_bankbuf[outpos++] = (ovr<<25) + e;
          //  tmp_bankbuf[outpos++] = ((cc<<10) | cfd); // cfd:10 bits
            tmp_bankbuf[outpos++] = ((cc<<10) | (timestamp&0x3ff)); // Save lower 10 bits of 2ns timestamp. Largest value is 2048ns
            tmp_bankbuf[outpos++] = (0xE<<28) + ((lost&0x3fff)<<14)
                              + (tmp_bankbuf[out_start+3] & 0x3fff);
         }
      }
   }
   return( outpos );
}
//       mod=1=>grif16  [words No-including-waveforms]
//          XXXX   XXX X XXXX   XXXX     XXXX   XXXX   XXXX   XXXX
// Header  <-8->  |Mod]#Words| <----------------ADDR------->  Dtype
// NetPkt    D    ------------------------------------------------
//  Mst1    00 14bit-filter-pattern  Wvfm          5bits:NumPileup
//  Mst2    0  31bit ID
// Trig      9     28bit trig_req
// Trig      A     28bit timestamp
// Trig      B     14bit deadtime  14bit timestamp
//  No waveform samples counter exists
//           1'b0, 5bits:integ[13:9], 26bits:pulseheight[25:0]
//           1'b0, 9bits:integ[ 8:0], 22bits:cfd[21:0]

//////////////////////////////////////////////////////////////////////////
///////         copy bank contents to buffer for other thread      ///////
//////////////////////////////////////////////////////////////////////////

unsigned bankbuf[BANK_BUFSIZE];
volatile unsigned long bankbuf_wrpos;
volatile unsigned long bankbuf_rdpos;
int copy_bank(unsigned *ptr, int size)
{
   int words_used, words_left, wrpos, words_to_end;
   unsigned int usecs=100;

   while(1){
      words_used = bankbuf_wrpos - bankbuf_rdpos;
      words_left = BANK_BUFSIZE - words_used;
      if( size > words_left ){ usleep(usecs); continue; }
      wrpos = bankbuf_wrpos % BANK_BUFSIZE;
      words_to_end = BANK_BUFSIZE - wrpos;
      if( size < words_to_end ){ // no wrap
         memcpy((char *)(bankbuf+wrpos), ptr, 4*size);
      } else { // write at end and wrap to start of bankbuf
         memcpy((char *)(bankbuf+wrpos), ptr,  4*words_to_end);
         memcpy((char *)(bankbuf),       ptr + words_to_end,
                                                    4*(size - words_to_end));
      }
      bankbuf_wrpos += size; return(0);
   }
}

//////////////////////////////////////////////////////////////////////////
///////     midas-main (usually run in separate thread)            ///////
//////////////////////////////////////////////////////////////////////////
// read data banks from midas file (additionally handle odb record at BOR)
//   either copy data to buffer or process it immediately[if single thread]
void midas_main(Sort_status *arg)
{
   int items, len, single_thread = arg->single_thread; //save!
   unsigned int usecs=100, *ptr;
   //static int evcount;
   char *bank_name;
   time_t tstamp;
   extern int process_grif3_bank(unsigned *buf, int len); // if single thread

   bankbuf_wrpos = bankbuf_rdpos = 0;
   recbufpos = recordlen = 0;
   while(1){
      if( arg->shutdown_midas != 0 ){ fprintf(stderr,"SHUTDOWN\n"); break; }
      if( arg->end_of_data == 1 ){ usleep(usecs); continue; }
      // buffer read position is taken care of in next_event and next_bank
      if( next_event(arg) < 0 ){     // midas event
         if( open_next_subrun(arg) == 0 ){
            recbufpos = recordlen = 0; continue;
         }
         arg->end_of_data = 1; continue;
      }
      diagnostics.midas_last_timestamp = ev_head.timestamp;
      if( diagnostics.midas_run_start == 0 ){
         diagnostics.midas_run_start = ev_head.timestamp;
      } else if( ev_head.timestamp != diagnostics.midas_run_start ){
         diagnostics.midas_datarate = diagnostics.midas_file_bytes /
            (ev_head.timestamp - diagnostics.midas_run_start);
      } else {
         ;// no datarate yet
      }
      while( (items = next_bank(arg, &bank_name)) > 0 ){
         len = ( (bank_head.data_size+3) & ~3) / 4; // bytes to ints
         ptr = (unsigned *)(recbuf+recbufpos);
         if( strcmp(bank_name,"CAEN") == 0 ){
            if( (len = translate_caen_bank(ptr, len)) <= 0 ){ continue; }
            ptr = tmp_bankbuf;
         } else if( strcmp(bank_name,"GRF3") != 0 &&
                    strcmp(bank_name,"GRF4") != 0 ){ continue; }//ignore othrs
         if( single_thread ){
            process_grif3_bank(ptr, len);
         } else { copy_bank(ptr, len); }
      }
      if( items == -2 ){
         if( arg->odb_ready == 0 ){ // odb "event" @ start of file
            read_odb_items(*(int *)bank_name, (int *)(recbuf+recbufpos) );
            arg->odb_ready = 1;
         }
         recbufpos += ev_head.data_size;
      } else if( arg->odb_ready == 0 ){
         arg->odb_ready = 2; // missing
      }
   }
   fprintf(stdout,"shutting down midas thread ...\n");
   return;
}

int next_record(Sort_status *arg)
{
   int bytes;

   // any data remaining, copy to start, before reading next block
   if( recordlen > recbufpos ){
      memmove(recbuf, recbuf+recbufpos, recordlen-recbufpos);
      recordlen -= recbufpos;
      recbufpos = 0;
   } else {
      recbufpos = recordlen = 0;
   }
   bytes = fread(recbuf+recordlen, 1, RECORDSIZE-recordlen, arg->data_fp);
   recordlen += bytes; diagnostics.midas_file_bytes += bytes;
   arg->midas_bytes += bytes;
   // printf("Read %ld bytes\n", diagnostics.midas_file_bytes);
   if( bytes <= 0 ){
      fprintf(stderr,"EOF at %ld on data_fp [%ld]\n", ftell(arg->data_fp), diagnostics.midas_file_bytes );
   }
   return(bytes);
}

#define MIDAS_HDRLEN sizeof(Midas_event_header)

static int evbase; // recbufpos at start of event
int next_event(Sort_status *arg)
{
   int bytes, bytes_done, bytes_avail, bytes_remain;

   // check if full event available, if not, grab new record
   while( recordlen-recbufpos <                                 MIDAS_HDRLEN ||
          recordlen-recbufpos < *(int *)(recbuf+recbufpos+12) + MIDAS_HDRLEN ){
      //printf("nextEvent: evlen=%6d bank:[wr:%ld rd:%ld]\n", *(int *)(recbuf+recbufpos+12), bankbuf_wrpos, bankbuf_rdpos );
      if( next_record(arg) <= 0 ){ return(-1); }
   }
   memcpy((char *)&ev_head, recbuf+recbufpos, sizeof(Midas_event_header) );
   recbufpos += sizeof(Midas_event_header);
   arg->midas_timestamp = ev_head.timestamp;
   if( arg->debug ){
      printf("\n\nEvent     id: %d",             ev_head.event_id     );
      printf("    Trigger mask: %d",             ev_head.trigger_mask );
      printf("    Serial   num: %d\ntime",       ev_head.serial_num   );
      printf("stamp   :      %s\n",ctime((time_t*)&ev_head.timestamp) );
      printf("data size   : %d\n",               ev_head.data_size    );
   }
   first_bank = 1; evbase = recbufpos;
   return(0);
}

#define BANK_IS_32BIT (1<<4)

int next_bank(Sort_status *arg, char **bank_name) // loop over banks in event
{
   Midas_bank16_header bank16_head;
   static char bank[5];
   char format[64];
   int items;

   if( first_bank ){ first_bank = 0;
      if( ev_head.event_id & 0x8000 ){
         *(int *)bank = ev_head.data_size; // temp store length here
         *bank_name = bank;
         return(-2);
      } /* "odb" event */
      memcpy((char *)&allbank_head, recbuf+recbufpos, sizeof(Midas_allbank_header) );
      swap_required = allbank_head.flags > 0x10000;
      if( swap_required ){
         swapInt( (char *)&allbank_head, 2*sizeof(int) );
      }
      if( arg->debug ){
         printf("    Allbank Size : %d",   allbank_head.allbanksize );
         printf("    Flags        : %d",   allbank_head.flags );
         printf("    (swap needed): %s\n", swap_required ? "yes" : "no" );
      }
      recbufpos += sizeof(Midas_allbank_header);
      memset((char *)&bank_head, 0, sizeof(Midas_bank_header) );
   }
   /* NOTE: size is rounded up to next multiple of 8 bytes (4 shorts) */
   /* before reading first bank, data_size is zero */
   recbufpos +=  ((bank_head.data_size + 7) & ~7);
   if( recbufpos - evbase >= ev_head.data_size ||
       recbufpos - evbase >= allbank_head.allbanksize ){ // end of event
      if( recbufpos - evbase != ev_head.data_size ){
         fprintf(stderr,"%d wrong event length %d != %d\n",
		 errcount++, recbufpos - evbase, ev_head.data_size );
	 return(-1);
      }
      return(-1);
   }
   if( allbank_head.flags & BANK_IS_32BIT ){
      memcpy((char *)&bank_head, recbuf+recbufpos, sizeof(Midas_bank_header) );
      recbufpos += sizeof(Midas_bank_header);
      if( swap_required ){
         swapInt( (char *)&bank_head.data_type, 2*sizeof(int) );
      }
   } else {
     /* old 16bit items - just copy items to new 32bit header */
      memcpy((char*)&bank16_head,recbuf+recbufpos,sizeof(Midas_bank16_header));
      recbufpos += sizeof(Midas_bank16_header);
      if( swap_required ){
         swapShort( (char *)&bank16_head.data_type, 2*sizeof(short) );
      }
      memcpy((char *)&bank_head, (char *)&bank16_head, 4 );
      bank_head.data_type = bank16_head.data_type;
      bank_head.data_size = bank16_head.data_size;
   }
   if( bank_head.data_size < 0 ){
      fprintf(stderr,"bank length < 0\n"); return(-1);
   }
   if( arg->debug ){
      printf("  Bank Name : %4s", bank_head.name );
      printf("  Type : %d",  bank_head.data_type );
      printf("  size : %d\n",  bank_head.data_size );
   }
   switch( bank_head.data_type ){
   case 4:
      if( swap_required ){ swapShort(recbuf+recbufpos, bank_head.data_size ); }
      if( arg->debug ){ sprintf(format,"    %%s[%%2d] = %%6u (0x%%04x)"); }
      items = bank_head.data_size/sizeof(short); break;
   case 6: /* scaler byte order seems to be wrong */
      if( swap_required  ){ swapInt(recbuf+recbufpos, bank_head.data_size ); }
      if( arg->debug ){ sprintf(format,"    %%s[%%2d] = %%10u (0x%%08x)"); }
      items = bank_head.data_size/sizeof(int); break;
   case 9:
      if( swap_required ){ swapInt(recbuf+recbufpos, bank_head.data_size ); }
      if( arg->debug ){ sprintf(format,"    %%s[%%2d] = %%6.2f (0x%%08x)"); }
      items = bank_head.data_size/sizeof(int); break;
   default:
      fprintf(stderr,"%s bank has unknown type: %d ... ignoring\n",
              bank_head.name, bank_head.data_type );
      return(-1);
   }
   memcpy(bank, bank_head.name, 4); bank[4]='\0';
   *bank_name = bank;
   return( items );
}
