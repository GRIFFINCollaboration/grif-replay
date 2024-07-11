#include <stdio.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

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
   printf("MIDAS  : in %10d runtime:%3ds =>    [%8.4f Mbytes/s]\n",
      diagnostics.midas_file_bytes, dt, rate/1024/1024
   );
   dt = current_time - diagnostics.run_sort_start;
   rate = diagnostics.midas_file_bytes / ((dt == 0) ? 1 : dt);
   printf("%22s sorttime:%3ds => [%8.4f Mbytes/s]\n", "", dt, rate/1024/1024);
}

//////////////////////////////////////////////////////////////////////////
///////         copy bank contents to buffer for other thread      ///////
//////////////////////////////////////////////////////////////////////////

unsigned bankbuf[BANK_BUFSIZE];
volatile unsigned long bankbuf_wrpos; 
volatile unsigned long bankbuf_rdpos;
int copy_bank()
{
   int size, words_used, words_left, wrpos, words_to_end;
   unsigned int usecs=100;

   size = ( (bank_head.data_size + 3) & ~3) / 4; // bytes to ints
   while(1){
      words_used = bankbuf_wrpos - bankbuf_rdpos;
      words_left = BANK_BUFSIZE - words_used;
      if( size > words_left ){ usleep(usecs); continue; }
      wrpos = bankbuf_wrpos % BANK_BUFSIZE;
      words_to_end = BANK_BUFSIZE - wrpos;
      if( size < words_to_end ){ // no wrap
         memcpy((char *)(bankbuf+wrpos), recbuf+recbufpos, 4*size);
      } else { // write at end and wrap to start of bankbuf
         memcpy((char *)(bankbuf+wrpos), recbuf+recbufpos,  4*words_to_end);
         memcpy((char *)(bankbuf),       recbuf+recbufpos + 4*words_to_end,
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
   int items, single_thread = arg->single_thread; //save!
   unsigned int usecs=100;
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
         if( strcmp(bank_name,"GRF3") == 0||strcmp(bank_name,"GRF4") == 0 ){
            if( single_thread ){
               process_grif3_bank((unsigned *)(recbuf+recbufpos), bank_head.data_size/4);
            } else { copy_bank(); }
         }
      }
      if( items == -2 ){
         if( arg->odb_ready == 0 ){ // odb "event" @ start of file
            read_odb_items(*(int *)bank_name, recbuf+recbufpos);
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

