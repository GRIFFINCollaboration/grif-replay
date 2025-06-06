// conpile with -fPIC
// possibly compile rest of program with -rdynamic
//   -rdynamic is needed if the loaded shared object refers to symbols
//             in the rest of the program

// symlink midas.h and libmidas.a to currentdir
// gcc -g -fPIC -rdynamic -fPIC -shared -o midas_module.so midas_module.c libmidas.a -lrt -lz -lutil -lnsl -lpthread

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
#include "grif-format.h"
#include "config.h"
#include "histogram.h"
#include "midas-format.h"
#include "midas.h"

extern Sort_metrics diagnostics;
extern char midas_runtitle[SYS_PATH_LENGTH];

static int midas_local_run_state;
static int single_thread_flag;
static Grif_event grif_event;
int unpack_griffin_bank( unsigned *data, int words ){ }
int process_decoded_fragment( Grif_event *tmp ){ }

void midas_event_callback(HNDLE buf, HNDLE req_id, EVENT_HEADER *pheader, void *pevent)
{
   int size, words, evlen, words_used, words_left, wrpos, words_to_end;
   unsigned int usecs=100, *data;

   if( midas_local_run_state == 0 ){ return; }

   diagnostics.midas_last_timestamp = pheader->time_stamp;
   if( diagnostics.midas_run_start == 0 ){
      diagnostics.midas_run_start = pheader->time_stamp;
   }
   if( (words = bk_locate(pevent, "GRF4", (DWORD *)&data)) <= 0 ){
      return;
   }
   diagnostics.midas_file_bytes += 4*words;
   if( single_thread_flag ){ // sort the event
      // NOT PROPERLY IMPLEMENTED YET
      while(words > 0 ){ // loop over fragments in bank till words all used
         if( (evlen = unpack_griffin_bank( data, words )) < 0 ){ break; }
         data += evlen; words -= evlen;
         if( process_decoded_fragment( &grif_event ) ){
            words = 0; break; // skipping all remaining data in bank
         }
      }
   } else { // copy bank to bank buffer and return
      // bankbuf is unsigned[BANK_BUFSIZE = 4M] 4Mword, 16Mbytes
      size = words;
      while(1){
         words_used = bankbuf_wrpos - bankbuf_rdpos;
         words_left = BANK_BUFSIZE - words_used;
         if( size > words_left ){ usleep(usecs); continue; }
         wrpos = bankbuf_wrpos % BANK_BUFSIZE;
         words_to_end = BANK_BUFSIZE - wrpos;
         if( size < words_to_end ){ // no wrap
            memcpy((char *)(bankbuf+wrpos), data, 4*size);
         } else { // write at end and wrap to start of bankbuf
            memcpy((char *)(bankbuf+wrpos), data,  4*words_to_end);
            memcpy((char *)(bankbuf),       data + 4*words_to_end,
                                          4*(size - words_to_end));
         }
         bankbuf_wrpos += size; return;
      }
   }
}

static INT midas_bor(INT run_number, char *err)
{
   Sort_status *arg = get_sort_status();
   arg->run_number      = run_number;
   arg->run_in_progress = midas_local_run_state = 1;
   diagnostics.midas_file_bytes = 0;
   printf("MIDAS BOR\n");
   return CM_SUCCESS;
}

static INT midas_eor(INT run_number, char *err)
{
   Sort_status *arg = get_sort_status();
   arg->run_in_progress = midas_local_run_state = 0;
   arg->end_of_data = 1;
   printf("MIDAS EOR\n");
   return CM_SUCCESS;
}

// Original midas_main() ...
//    zero bankbuf_rd/wrpos, recbufpos/len
//    while(1){ check shutdown/EOF
//       if nxt_event <0 { open nxt subrun or set EOF }
//       while( items=next_bank() > 0 ){
//         process or copy, also handle odb_event
//       } 
//    }
// this version ...
//    will mostly idle in cm_yield, event callback should just copy to bankbuf
//       (*could* also process event in this thread, if singlethread)
//    directly read tables from odb
//    zero stuff (inc histos) on BOR
//    save histos             on EOR
//
//
extern char midas_host[64];
extern char midas_expt[64];

static char data_dir[256];
static char run_title[256];

extern int debug, odb_daqsize;
extern char chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
extern short   *addrs;
extern int    *dtypes;
extern float   *gains;
extern float *offsets;
extern float   *quads;
void midas_module_main(Sort_status *arg)
{
   int i, runnum, size, state, status, event_id=1, trigger_id=TRIGGER_ALL;
   HNDLE hDB, hKey, hBuf, req_id;
   char key[64];
  
   single_thread_flag = arg->single_thread;
   bankbuf_wrpos = bankbuf_rdpos = 0;
   //debug = 1;

   if( cm_connect_experiment(midas_host, midas_expt, "Grif_Replay", NULL)
                                                     != CM_SUCCESS ){
      fprintf(stderr,"Cant connect to experiment\n"); return;
   }
   if( cm_get_experiment_database(&hDB, NULL) != CM_SUCCESS ){
      fprintf(stderr,"Cant open ODB\n");
      cm_disconnect_experiment(); return;
   }
   
   // read odb tables and set arg->odb_ready
   // (later will re-read on every run start)
   // (will have to monitor for changes, and if so, redo histos)
   // ALSO CALL user_sort_init() on new run
   //    (or at end of run if sorting finished)
   sprintf(key,"/DAQ/PSC/PSC");
   if( (status=db_find_key(hDB, 0, key, &hKey)) != DB_SUCCESS){
     cm_msg(MINFO,"Replay","Key %s not found", key); return;
   }
   size = 0; db_get_record_size(hDB,hKey,0,&size);
   odb_daqsize = size/sizeof(short);
   size = MAX_DAQSIZE*sizeof(short);
   if( (db_get_data(hDB,hKey,addrs,&size,TID_SHORT)) != DB_SUCCESS){
      cm_msg(MINFO,"Replay","Can't get data for Key %s", key); return;
   }
   
   size = MAX_DAQSIZE*CHAN_NAMELEN;
   sprintf(key,"/DAQ/PSC/chan");
   if( (status=db_find_key(hDB, 0, key, &hKey)) != DB_SUCCESS){
     cm_msg(MINFO,"Replay","Key %s not found", key); return;
   }
   if( (db_get_data(hDB,hKey,chan_name,&size,TID_STRING)) != DB_SUCCESS){
      cm_msg(MINFO,"Replay","Can't get data for Key %s", key); return;
   }
   
   size = MAX_DAQSIZE*sizeof(int);
   sprintf(key,"/DAQ/PSC/datatype");
   if( (status=db_find_key(hDB, 0, key, &hKey)) != DB_SUCCESS){
     cm_msg(MINFO,"Replay","Key %s not found", key); return;
   }
   if( (db_get_data(hDB,hKey,dtypes,&size,TID_INT)) != DB_SUCCESS){
      cm_msg(MINFO,"Replay","Can't get data for Key %s", key); return;
   }
   
   sprintf(key,"/DAQ/PSC/gain");
   if( (status=db_find_key(hDB, 0, key, &hKey)) != DB_SUCCESS){
     cm_msg(MINFO,"Replay","Key %s not found", key); return;
   }
   if( (db_get_data(hDB,hKey,gains,&size,TID_FLOAT)) != DB_SUCCESS){
      cm_msg(MINFO,"Replay","Can't get data for Key %s", key); return;
   }
   
   sprintf(key,"/DAQ/PSC/offset");
   if( (status=db_find_key(hDB, 0, key, &hKey)) != DB_SUCCESS){
     cm_msg(MINFO,"Replay","Key %s not found", key); return;
   }
   if( (db_get_data(hDB,hKey,offsets,&size,TID_FLOAT)) != DB_SUCCESS){
      cm_msg(MINFO,"Replay","Can't get data for Key %s", key); return;
   }
   
   sprintf(key,"/DAQ/PSC/quadratic");
   if( (status=db_find_key(hDB, 0, key, &hKey)) != DB_SUCCESS){
     cm_msg(MINFO,"Replay","Key %s not found", key); return;
   }
   if( (db_get_data(hDB,hKey,quads,&size,TID_FLOAT)) != DB_SUCCESS){
      cm_msg(MINFO,"Replay","Can't get data for Key %s", key); return;
   }
   gen_derived_odb_tables();
   arg->odb_ready = 1;

   // last arg is sequencing number[1-1000]
   if( cm_register_transition(TR_START, midas_bor, 998) != CM_SUCCESS||
       cm_register_transition(TR_STOP , midas_eor, 999) != CM_SUCCESS){
      cm_msg(MINFO,"Replay","Can't register run transitions"); return;
   }

   size = 4; sprintf(key,"/Runinfo/Run number");
   db_get_value(hDB, 0, key, &runnum, &size, TID_INT, FALSE);
   size = 4; sprintf(key,"/Runinfo/state");
   db_get_value(hDB, 0, key, &state, &size, TID_INT, FALSE);
   arg->run_number      = runnum;
   arg->run_in_progress = midas_local_run_state = ( state != STATE_STOPPED );

   size = 256; sprintf(key,"/Experiment/Run_parameters/Run_Title");
   db_get_value(hDB, 0, key, &run_title, &size, TID_STRING, FALSE); run_title[255]=0;
   size = strlen(run_title); 
   memcpy( midas_runtitle, run_title, size ); midas_runtitle[size] = 0;
   
   size = 256; sprintf(key,"/Logger/Data dir");
   db_get_value(hDB, 0, key, &data_dir, &size, TID_STRING, FALSE);
   set_directory(configs[0], "Histo", data_dir);
   set_directory(configs[1], "Histo", data_dir);
  
   bm_open_buffer("SYSTEM", (32*1024*1024), &hBuf);
   bm_set_cache_size(hBuf, 100000, 0); // 100kbytes
   bm_request_event(hBuf, event_id, trigger_id, GET_NONBLOCKING,
                                     &req_id, midas_event_callback );
   do {
      status = cm_yield(1000);
   } while( status != RPC_SHUTDOWN && status != SS_ABORT );
   bm_delete_request(req_id); // done automatically on bm_close_buffer
   bm_close_buffer(hBuf);
   cm_disconnect_experiment();
   fprintf(stdout,"shutting down midas_module thread ...\n");
   return;
}
//   cm_get_env(host, expt)
//   cm_connect_experiment1(hostm expt, ana_name, NULL, odbsize, wd_timeout)
//   cm_get_expt_database(&hDB,NULL)
//   register transitions[start,stop,pause,resume]
//   bm_open_buffer
//   bm_set_cache_size() cache 100k bytes (avoid event-by-event calls)
//   bm_request_event() [analyzer.c:ev_id=1,trig=ALL,NONBLOCK,SYSTEM]
//   loop_online()
