/* read midas data files
 set sources=( grif-replay.c midas-format.c grif-format.c histogram.c web_server.c config.c reorder.c user_sort.c default_sort.c test_config.c )
 gcc -g     -o grif-replay $sources -rdynamic -ldl -lm -lpthread
 gcc -g -O3 -o grif-replay $sources -rdynamic -ldl -lm -lpthread
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "config.h"
#include "histogram.h"
#include "grif-format.h"
#include "midas-format.h"

int sort_next_file(Config *cfg, Sort_status *sort);

int coinc_events_cutoff = 25;
char midas_runtitle[SYS_PATH_LENGTH];
Sort_metrics diagnostics;
static Sort_status sort_status;
volatile int shutdown_server = 0;
static pthread_t web_thread;
static void online_loop(Config *cfg, Sort_status *sort);
int main(int argc, char *argv[])
{
   Sort_status *sort = &sort_status;
   int web_arg=1;  Config *cfg;

   sort->single_thread = 0;
   pthread_create(&web_thread, NULL,(void* (*)(void*))web_main, &web_arg);
   while( !shutdown_server ){ // monitor file queue and sort any added files

      if( sort->current_filenum == sort->final_filenum ){ sleep(1); continue; }
      copy_config(configs[0], configs[1]); // copy config0 to cfg1 for sorting
      cfg = configs[1];
      if( sort->online_mode ){
         online_loop(cfg, sort);
      } else
      if( open_next_sortfiles(sort) == 0 ){
         sort_next_file(cfg, sort);
         fprintf(stdout,"DONE\n");
         close_sortfiles(sort);
      }
      if( ++sort->current_filenum == FILE_QLEN ){ sort->current_filenum = 0; }
   }
   pthread_join(web_thread, NULL);
   exit(0);
}

extern int frag_hist[MAX_COINC_EVENTS];
void show_coinc_stats()
{
   char tmp[64];
   int i, sum;
   for(i=0; i<15; i++){
      sprintf(tmp, "Coinc[%8d]:%d", i, frag_hist[i]);
      printf("%-25s%c", tmp, (((i+1)%3)==0)?'\n':' ');
   }
   sum=0; for(; i<30; i++){ sum +=  frag_hist[i]; }
   printf("Coinc[ 15-  29]:%d\n", sum);
   sum=0; for(; i<100; i++){ sum +=  frag_hist[i]; }
   printf("Coinc[ 30-  99]:%d\n", sum);
   sum=0; for(; i<250; i++){ sum +=  frag_hist[i]; }
   printf("Coinc[100- 249]:%d\n", sum);
   sum=0; for(; i<500; i++){ sum +=  frag_hist[i]; }
   printf("Coinc[250- 499]:%d\n", sum);
   sum=0; for(; i<MAX_COINC_EVENTS; i++){ sum +=  frag_hist[i]; }
   printf("Coinc[500-%4d]:%d\n", MAX_COINC_EVENTS, sum);
}

static int presort_window_start, sort_window_start;
static int done_events;
extern void grif_main(Sort_status *arg);
extern void reorder_main(Sort_status *arg);
extern void reorder_out(Sort_status *arg);
extern void sort_main(Sort_status *arg);
static pthread_t midas_thread, grif_thread, ordthrd, ordthr2;
static int reorder_save, singlethread_save, sortthread_save;
extern int (*midas_module_main)(Sort_status *);
int sort_next_file(Config *cfg, Sort_status *sort)
{
   time_t end, start=time(NULL);
   done_events = 0;
   presort_window_start = sort_window_start = 0;
   memset(&diagnostics, 0, sizeof(Sort_metrics) );
   diagnostics.run_sort_start = start;
   sort->shutdown_midas = sort->end_of_data = 0;
   sort->reorder_in_done = sort->reorder_out_done = 0;
   sort->grif_sort_done = sort->odb_done = sort->odb_ready = 0;
   reorder_save = sort->reorder;
   singlethread_save = sort->single_thread;
   sortthread_save = sort->sort_thread;
   if( singlethread_save == 1 ){
      if( sort->online_mode ){
         midas_module_main(sort);
      } else {
         midas_main(sort);
      }
   } else {
      printf("creating midas thread\n");
      if( !sort->online_mode ){
         pthread_create(&midas_thread, NULL, (void* (*)(void*))midas_main, sort);
      } else {
         pthread_create(&midas_thread, NULL, (void* (*)(void*))midas_module_main, sort);
      }
      printf("creating reorder threads\n");
      pthread_create(&ordthrd,NULL,(void* (*)(void*))reorder_main,sort);
      pthread_create(&ordthr2,NULL,(void* (*)(void*))reorder_out, sort);

      while( !sort->odb_ready ){     // wait for midas thread to read odb event
         usleep(1);
      }
      //user_sort_init();   // user histos already defined, but defaul histos
      init_default_histos(configs[1], sort);     // depend on odb in datafile
      pthread_create(&grif_thread, NULL,(void* (*)(void*))grif_main, sort);

      sort_main(sort); // this exits when sort is done

      sort->shutdown_midas = 1;
      pthread_join(midas_thread, NULL);
      pthread_join(ordthrd, NULL);
      pthread_join(ordthr2, NULL);
      pthread_join(grif_thread, NULL);
   }
   end=time(NULL);
   cfg->midas_start_time = diagnostics.midas_run_start;
   cfg->midas_runtime    = diagnostics.midas_last_timestamp+1;
   cfg->midas_runtime   -= cfg->midas_start_time;
   memcpy(cfg->midas_title, midas_runtitle, SYS_PATH_LENGTH);
   if( !sort->online_mode ){
      write_histofile(cfg, sort->histo_fp);
   } else {
      unload_midas_module();
   }
   printf("File took %ld seconds\n", end-start);
   show_coinc_stats();
   return(0);
}

static void online_loop(Config *cfg, Sort_status *sort)
{
   time_t end, start;
   char tmp[64];
   FILE *fp;
   pthread_create(&midas_thread, NULL, (void* (*)(void*))midas_module_main, sort);
   sort->odb_ready = 0; // only done on module load atm.
   while(1){
      while( !sort->run_in_progress ){ usleep(1); } // wait for run to start
      start =time(NULL);
      done_events = 0;
      memset(&diagnostics, 0, sizeof(Sort_metrics) );
      diagnostics.run_sort_start = start;
      sort->shutdown_midas = sort->end_of_data = 0;
      sort->reorder_in_done = sort->reorder_out_done = 0;
      sort->grif_sort_done = sort->odb_done = 0;
      reorder_save = sort->reorder;
      singlethread_save = sort->single_thread;
      sortthread_save = sort->sort_thread;
      presort_window_start = sort_window_start = 0;

      printf("creating reorder threads\n");
      pthread_create(&ordthrd,NULL,(void* (*)(void*))reorder_main,sort);
      pthread_create(&ordthr2,NULL,(void* (*)(void*))reorder_out, sort);

      while( !sort->odb_ready ){     // wait for midas thread to read odb event
         usleep(1);
      }
      //user_sort_init();   // user histos already defined, but defaul histos
      init_default_histos(configs[1], sort);     // depend on odb in datafile
      pthread_create(&grif_thread, NULL,(void* (*)(void*))grif_main, sort);

      sort_main(sort); // this exits when sort is done

      sort->shutdown_midas = 1;
      // DO NOT SHTUDOWN MIDAS THREAD
      pthread_join(ordthrd, NULL);
      pthread_join(ordthr2, NULL);
      pthread_join(grif_thread, NULL);

      end=time(NULL);
      cfg->midas_start_time = diagnostics.midas_run_start;
      cfg->midas_runtime    = diagnostics.midas_last_timestamp+1;
      cfg->midas_runtime   -= cfg->midas_start_time;
      memcpy(cfg->midas_title, midas_runtitle, SYS_PATH_LENGTH);
      fp = NULL;
      if( strlen(cfg->histo_dir) > 0 ){
         sprintf(tmp,"%s/run%05d.tar", cfg->histo_dir, sort->run_number);
         if( (fp=fopen(tmp,"w")) != NULL ){
            write_histofile(cfg, fp); fclose(fp);
         } else {
            printf("Can't open histo file: %s to write\n", tmp);
         }
      }
      if( fp == NULL ){ // no histo_dir or not writable - use cwd
         sprintf(tmp,"./run%05d.tar", sort->run_number);
         if( (fp=fopen(tmp,"w")) != NULL ){
            write_histofile(cfg, fp); fclose(fp);
         } else {
            printf("Can't open histo file: %s to write\n", tmp);
         }
      }
      copy_config(configs[0], configs[1]); // prepare for next run
      printf("File took %ld seconds\n", end-start);
      show_coinc_stats();
   }
   unload_midas_module();
   return;
}

extern volatile long grifevent_wrpos;
volatile long grifevent_rdpos;
extern Grif_event grif_event[MAX_COINC_EVENTS];

extern volatile unsigned long bankbuf_wrpos;
extern volatile unsigned long bankbuf_rdpos;
extern volatile long tsevents_in;
extern long tsevents_out;
extern volatile unsigned long eventbuf_rdpos;
extern volatile unsigned long eventbuf_wrpos;
extern volatile long grif_evcount;
void show_sort_state()
{
   int val = sort_status.midas_bytes/1000;
   int v2 = bankbuf_wrpos/1000000, v3 = bankbuf_rdpos/1000000;
   int v4 = tsevents_in/1000, v5 = tsevents_out/1000;
   int v6 = eventbuf_wrpos/1000000, v7 = eventbuf_rdpos/1000000;
   int v8 = grifevent_wrpos/1000, v9 = grifevent_rdpos/1000;
   int v10 = grifevent_wrpos - grifevent_rdpos;
   printf("MIDAS:read %d Mbytes [~%d Kevents]\n", val/1000, val/50);
   printf("      BUF in:%dMbytes out:%dMbytes  [Cap:%5.1f%%]\n",
          4*v2, 4*v3, (100000000.0*(v2-v3))/BANK_BUFSIZE );
   printf("REORDER: In:%dKevents Out%dKevents\n", v4, v5 );
   printf("     BUF In:%dMbytes  Out:%dMbytes  [Cap:%5.1f%%]\n",
          4*v6, 4*v7, (100000000.0*(v6-v7))/EVENTBUFSIZE );
   printf("GRIF: Unpacked:%dKevents  ", (int)(grif_evcount/1000) );
   printf("      Sorted  :%dKevents\n", done_events/1000 );
   printf("     BUF In:%dKevents Out%dKevents [=%d][Cap:%5.1f%%]\n\n",
          v8, v9, v10, (100.0*v10)/MAX_COINC_EVENTS);
}
Sort_status *get_sort_status(){ return( &sort_status ); }

void sort_main(Sort_status *arg)
{
   int i, len, nxtpos, rd_avail;
   static long grifevent_nxtpos;
   unsigned int usecs=100;

   printf("starting sort_main ...\n");
   grifevent_rdpos = grifevent_nxtpos = nxtpos = 0;
   while(1){
      // if( arg->shutdown_midas != 0 ){  break; }
      rd_avail = grifevent_wrpos - grifevent_nxtpos;
      if( arg->grif_sort_done && rd_avail < 1 ){ break; }
      if( rd_avail < 1 ){ usleep(usecs); continue; }
      process_event(&grif_event[nxtpos], nxtpos);
      nxtpos = ++grifevent_nxtpos % MAX_COINC_EVENTS;
   }
   printf("sort_main finished\n");
   return;
}

static int proc_calls, sorted, skipped, prefull, sortfull, completed_events;
// called when each new event read into ptr -> list[slot]
int process_event(Grif_event *ptr, int slot)
{
   time_t cur_time = time(NULL);
   static int calls, prv_call;
   static time_t prv_time;
   int dt = cur_time - prv_time, de = calls - prv_call;
   if( prv_time == 0 ){ prv_time = cur_time; }
   prv_call = ++calls;

   if( cur_time-prv_time >= 10 ){
      printf("----------------------------------------------------------------\n");
      midas_status(cur_time); reorder_status(cur_time);  grif_status(cur_time);
      printf("ProcEvt: %10d[Good:%3d%% Skip:%3d%% WinFull:%3d%%] %6.3f Mevt/s\n",
             calls, (int)(100.0*(calls-skipped-prefull)/calls),
             (int)(100.0*skipped/calls), (int)(100.0*prefull/calls), (de/(1000000.0*dt))
      );
      prv_time = cur_time;
   }
   if( singlethread_save == 1 && sort_status.odb_done == 0 ){ calls = 0;
      init_default_histos(configs[1], &sort_status); sort_status.odb_done = 1;
   }
   apply_gains(ptr);
   insert_presort_win(ptr, slot);
   return(0);
}

char *debug_show_ts(long ts)
{
   static char tmp[32];
   int deci_ms = ts/10000;
   int remain = ts - (deci_ms*10000);

   sprintf(tmp,"%5.1fms+%04d", 00000.1*deci_ms, remain);
   return(tmp);
}
char *debug_show_chan(Grif_event *ptr)
{
   static char tmp[32];
   if( ptr->chan != -1 ){
      sprintf(tmp," %3d", ptr->chan);
   } else {
      sprintf(tmp,"%04x", ptr->address);
   }
   return(tmp);
}
extern char chan_name[MAX_DAQSIZE][CHAN_NAMELEN];


// add event to presort window (event has just been read in)
//    recalculate coincwin (sorting any events leaving window)
// => final win of run won't be sorted, as these events will not leave window
int insert_presort_win(Grif_event *ptr, int slot)
{
   int window_width = 500; // 5us to capture all pileup events - MAXIMUM (indiv. gates can be smaller)
   int win_count, win_end;
   Grif_event *alt;
   long dt;

   /*
   static int prv_evt[65535], count;
   ++count;
   i = (slot-1-window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
   printf("PRST:Chan[%4s:prev:%5d][%s] [win:%5d[%05d-%05d:%s]          ",
          debug_show_chan(ptr),
          prv_evt[ptr->address] == 0 ? -1 : count-prv_evt[ptr->address],
          debug_show_ts(ptr->ts),
          i, window_start, slot-1,
          debug_show_ts(grif_event[window_start].ts) );
   prv_evt[ptr->address] = count;
   if( ptr->chan != -1 ){
      printf("%s E=%6d[cal:%8.1f] id:0x%08x\n",
             chan_name[ptr->chan],ptr->energy,ptr->esum, ptr->master_id);
   } else {
      printf("---------- E=%6d[cal:%8.1f] id:0x%08x\n",
             ptr->energy, ptr->esum, ptr->master_id );
   }
   */

   ///////////////// Presort window (used for suppression/addback)
   while( presort_window_start != slot ){ alt = &grif_event[presort_window_start];
   win_count = (slot - presort_window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }

      // should exit while-loop when no more events outside window
      //    *BUT* add error recovery - if window too full, dump events
      if( dt < window_width ){
         //if( win_count < 0.45*MAX_COINC_EVENTS  ){ break; }
         if( win_count < coinc_events_cutoff ){ break; } // LIMIT to ?? events
         else { ++prefull; }
      }

      // event[win_start] is leaving window
      //    ( either because dt > coincwidth OR due to error recovery)
      // NOTE event[slot] is out of window - use slot-1 as window-end
      if( (win_end = slot-1) < 0 ){ win_end = MAX_COINC_EVENTS-1; } // WRAP
      pre_sort(presort_window_start, win_end);
      insert_sort_win(alt, presort_window_start); // add event to next window
      if( ++presort_window_start >= MAX_COINC_EVENTS ){ presort_window_start=0; } // WRAP
   }
   return(0);
}

// add event to main sort window (event has just left presort window)
int insert_sort_win(Grif_event *ptr, int slot)
{
   int window_width = 200; // 2us - MAXIMUM (indiv. gates can be smaller)
   int win_count, win_end;
   Grif_event *alt;
   long dt;

   /* i = (slot-1-window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
   printf("MAIN:Chan[%4s:    :     ][%s] [win:%5d[%05d-%05d:%s]\n",
          debug_show_chan(ptr),

          debug_show_ts(ptr->ts),
          i, window_start, slot-1,
          debug_show_ts(grif_event[window_start].ts) );
   */
   while( sort_window_start != slot ){ alt = &grif_event[sort_window_start];
       win_count = (slot - sort_window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }

      // should exit while-loop when no more events outside window
      //    *BUT* add error recovery - if window too full, dump events
      if( dt < window_width ){
       //if( win_count >= 0.45*MAX_COINC_EVENTS ){ ++sortfull; } else {
         if( win_count > coinc_events_cutoff ){ ++sortfull; } else {
            // now removed all events not in coinc with newly added fragment
            //   so can update coinc window counters with just-added frag
            user_addto_window(sort_window_start, slot);
            break;
         }
      }

      // event[win_start] is leaving window
      //    ( either because dt > coincwidth OR due to error recovery)
      // NOTE event[slot] is out of window - use slot-1 as window-end
      if( (win_end = slot-1) < 0 ){ win_end = MAX_COINC_EVENTS-1; } // WRAP
      if( alt->chan != -1 ){ ++sorted;
         default_sort(sort_window_start, win_end, SORT_ONE);
         //user_sort(window_start, win_end, SORT_ONE);
      } else { ++skipped; }
      if( ++sort_window_start >= MAX_COINC_EVENTS ){ sort_window_start=0; } // WRAP
      ++grifevent_rdpos;  ++completed_events;
   }
   return(0);
}


/////////////////////////////////////////////////////////////////////////////////
///////// Alternate sort method [build complete events, then sort them] /////////
/////////   [faster?] but slightly worse at finding all coincidences    /////////
/////////////////////////////////////////////////////////////////////////////////


// add event to window (event has just been read in)
//    (previously would loop over events already in window to see if they are now OUT)
// NEW method - if first event is OUT ...
//    build an event with all fragments before current fragment (which starts a new event)
int build_event(Grif_event *ptr, int slot)
{
   int window_width = 200; // 2us - MAXIMUM (indiv. gates can be smaller)
   int win_count, win_end;
   static int window_start;
   Grif_event *alt;
   long dt;

   if( window_start == slot ){ return(0); }

   alt = &grif_event[window_start];
   win_count = (slot - window_start+2*MAX_COINC_EVENTS) % MAX_COINC_EVENTS;
   dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }
   if( dt < window_width ){ return(0); }

   if( (win_end = slot-1) < 0 ){ win_end = MAX_COINC_EVENTS-1; } // WRAP
   sort_built_event(window_start, win_end);

   window_start = slot;
   return(0);
}

int sort_built_event(int window_start, int win_end)
{
   pre_sort(window_start, win_end); // only fold of first fragment is set
   default_sort(window_start, win_end, SORT_ALL);
   //user_sort(window_start, win_end, SORT_ALL);
   return(0);
}
