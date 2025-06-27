#ifndef GRIF_REPLAY_H
#define GRIF_REPLAY_H

#define SYS_PATH_LENGTH    2048 // allow this much space for filesystem paths
                                // also any other "large" strings

// tar file format limits some string lengths ...
#define HISTO_FOLDER_LENGTH 155 // stored in 155 char field
#define HANDLE_LENGTH       100 // stored in 100 char field
#define TITLE_LENGTH        100 // stored in 100 char field
//#define FILE_BUFSIZ     1048576 // 64k * 4 * 4 [1M]
#define FILE_BUFSIZ     (1024*1024*1024) // allows matrices up to 16k*16k*4
#define STRING_LEN          256

#define MAX_HISTO_GATES     64
#define MAX_SORT_VARS      256
#define MAX_HISTOGRAMS   16384
#define MAX_ADDRESS    0x10000
#define MAX_CONFIGS        256
#define MAX_CONDS         1024
#define MAX_DAQSIZE       1200 // max #channels that can be defined in odb
#define CHAN_NAMELEN        32 // size of midas odb strings
#define MAX_GATE_CONDS      64
#define MAX_GATES          256
#define MAX_GLOBALS        256
#define MAX_CALIB  MAX_DAQSIZE
#define SMALL_HISTO_BINS 65536 // 64k bins

#define SORT_ALL 0     // for built-events, sort all events in window at once
#define SORT_ONE 1     //     not-built, sort events one-by-one as they leave window

extern int debug;

#define FILE_QLEN 256
typedef struct sortfile_struct {
   char  *data_dir; char  *data_name;
   char *histo_dir; char *histo_name;
   char  *conf_dir; char  *conf_name;
   long data_size; int run; int subrun;
   int run_digits; int subrun_digits;
   int num_subruns;  char  *cal_src;
   char file_info[4][256]; char **arg;  int narg;  int carg;
   char recent_cal[256]; // others attempted on open-subrun
} Sortfile;

typedef struct sortstatus_struct {
   volatile int  end_of_data; volatile int odb_done;  int debug;
   volatile int  reorder_out_done;    int single_thread;
   volatile int  reorder_in_done;     int reorder;
   volatile int  current_filenum;     int odb_ready;
   volatile int  final_filenum;       int sort_thread;
   volatile long midas_bytes;         int cal_overwrite;
   volatile int  midas_timestamp;     FILE *data_fp;
   volatile int  shutdown_midas;      FILE *histo_fp;
   volatile int  grif_sort_done;      FILE *cal_fp;
   volatile int  online_mode;         int  sum_mode;
   volatile int  run_in_progress;
   volatile int  run_number;
} Sort_status;

extern void web_main(int *);
extern Sort_status *get_sort_status();

//////////////////////////// diagnostics ////////////////////////////////

#define REORDER_ERRTYPES 9
#define ERR_FORMAT_IN 0
#define ERR_WORDS_IN  1
#define ERR_INIT_IN   2
#define ERR_LENGTH_IN 3
#define ERR_LATE_IN   4
#define ERR_EARLY_IN  5
#define ERR_LATE_OUT  6
#define ERR_UNKNOWN_ORDER_IN  7
#define ERR_TS_DESYNC_IN      8

typedef struct sortmetrics_struct {
   int reorder_error[REORDER_ERRTYPES];
   long midas_file_bytes;
   int midas_run_start;
   int midas_datarate;
   int midas_last_timestamp;
   int run_sort_start;
} Sort_metrics;

extern void midas_status(int), reorder_status(int), grif_status(int);
extern int gen_derived_odb_tables();
extern int open_next_subrun(Sort_status *arg);
extern int read_odb_items(int len, int *bank_data);

#endif

// dynamically change #reorder_slots
