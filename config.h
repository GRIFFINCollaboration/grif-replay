// TODO ..
//   commented "inuse check for remove global
// remove_gate_from_group
//  gate/grouplist ordering may need to be reset on clearing
// config mtime to ALL new functions
//   commented free(sort->global[i]->name) due to crash

#ifndef CONFIG_H
#define CONFIG_H

#include "grif-replay.h"

int send_spectrum_list(char *name, int fd);
int send_spectrum(int num, char urlarg[][STRING_LEN], char *, int fd);
int send_sort_status(int fd);
int most_recent_calib_file(char *data_dir, int data_run, char *result);
int send_datafile_list(char *path, int fd, int type);
int send_histofile_list(char *path, int fd);
int send_configfile_list(char *path, int fd);
int send_file_details(char *path, int fd);
int add_sortfile(char *path, char *histodir, char *confdir, char *calsrc);
int open_next_sortfiles(Sort_status *arg);
int free_sortfile(Sortfile *sort);
int close_sortfiles(Sort_status *arg);
int end_current_sortfile(int fd);
void unload_midas_module();
int user_addto_window(int win_strt, int new_frag);
int default_sort(int win_idx, int frag_idx, int flag);
int sort_built_event(int window_start, int win_end);

// Sorting window sizes. Can be set as Globals
extern int presort_window_width;
extern int sort_window_width;

/////////////////////////////////////////////////////////////////////////
/////////////////////       Histograms       ////////////////////////////
/////////////////////////////////////////////////////////////////////////

#define DEFAULT_CONFIG "last.json"

typedef struct th1i_struct Histogram;

typedef struct global_struct {
   char name[STRING_LEN]; int min; int max;
} Global;

typedef struct cal_coeff_struct {
   char name[CHAN_NAMELEN]; float offset; float gain; float quad;
   float pileupk1[7], pileupk2[7], pileupE1[7];
   float crosstalk0[16], crosstalk1[16], crosstalk2[16];
   short address; short datatype;
} Cal_coeff;

typedef struct sortvar_struct {     // sortvars used by histos AND conditions
   int value;      int offset;    int dtype;     //  but "use_count" only used
   int use_count_x;  int valid;   int local;     //  to refer to histo_use
   Histogram *histo_list_x[MAX_HISTOGRAMS];      // (inc/dec in add/rmv_histo)
   char name[STRING_LEN]; char title[STRING_LEN];
} Sortvar;

typedef struct cond_struct {      // use count inc/dec when un/applying gates
   char name[STRING_LEN]; Sortvar *var; int op; int value;      // (to histos)
   int use_count;  int veto;
   int passed; int pass_count;  // passed only used during sort
} Cond;

typedef struct gate_struct {
   char name[STRING_LEN]; int nconds; Cond *conds[MAX_GATE_CONDS];
   int use_count;  int passed;
} Gate;

typedef struct histo_folder_struct {
   char name[HISTO_FOLDER_LENGTH];
   struct histo_folder_struct *next_subfolder;  // 1 level lower
   struct histo_folder_struct *next_folder;     // same level
   struct th1i_struct *first_histo;             // same level
} Folder;

typedef struct th1i_struct TH1I;
typedef struct th2i_struct TH2I;

#define SYMMETERIZE -1 // When the ybins are set to this for a 2D histogram it will be symmeterized

// float has around 24bits integer precision
struct th1i_struct {  long  file_data_offset;    int data_size;
   int      type;  TH1I    *next;   char    path[HISTO_FOLDER_LENGTH];
   int     xbins;  int     ybins;   char          title[TITLE_LENGTH];
   int     *data;  int valid_bins;  char        handle[HANDLE_LENGTH];
   int underflow;  int   overflow;
   int   entries;  int  num_gates;  char  *gate_names[MAX_HISTO_GATES];
   Sortvar *xvar;  Sortvar  *yvar;  int  *gate_passed[MAX_HISTO_GATES];
   int xmin; int xmax; int ymin; int ymax; int suppress; int user;
   int xrange; int yrange;  int done_flag;  int symm;
   int   (*Reset)(TH1I *);
   int   (*Fill)(TH1I *, int, int);
   int   (*SetBinContent)(TH1I *, int, int);
   int   (*GetBinContent)(TH1I *, int);
   int   (*SetValidLen  )(TH1I *, int);
};
struct th2i_struct {  long  file_data_offset;  int data_size;
   int      type;  TH1I    *next;   char     path[HISTO_FOLDER_LENGTH];
   int     xbins;  int     ybins;   char          title[TITLE_LENGTH];
   int     *data;  int valid_bins;  char        handle[HANDLE_LENGTH];
   int underflow;  int   overflow;
   int   entries;  int  num_gates;  char  *gate_names[MAX_HISTO_GATES];
   Sortvar *xvar;  Sortvar  *yvar;  int  *gate_passed[MAX_HISTO_GATES];
   int xmin; int xmax; int ymin; int ymax; int suppress;  int user;
   int xrange; int yrange;  int done_flag;  int symm;
   int   (*Reset)(TH2I *);
   int   (*Fill)(TH2I *, int, int, int);
   int   (*SetBinContent)(TH2I *, int, int, int);
   int   (*GetBinContent)(TH2I *, int, int);
   int   (*SetValidLen  )(TH2I *, int);
};

#define DISK_CONFIG 0
#define MEM_CONFIG  1

// group config element stuff together, with add/delete/copy fns
//   to allow quick switch between storing data or pointers (malloc/free)
typedef struct config_set_struct { int  type; // memory(live,sort) or disk
   FILE             *histo_fp;
   char name[SYS_PATH_LENGTH];     // "live", "sort", or pathname of tar file
   char data_dir[SYS_PATH_LENGTH];   // most recent datafile directory
   char histo_dir[SYS_PATH_LENGTH];  // most recent histofile directory
   char config_dir[SYS_PATH_LENGTH]; // most recent configfile directory
   char midas_title[SYS_PATH_LENGTH];// title of midas run(histo_config)
   char current_path[HISTO_FOLDER_LENGTH]; // current histogram path
   int  midas_start_time;  int midas_runtime;  int mtime;  int lock;
   int folders_valid;    Folder first_folder;
   int current_depth;    Folder *treepath[HISTO_FOLDER_LENGTH]; // length 2296
   int ncal;             Cal_coeff *calib[MAX_CALIB];
   int nglobal;          Global *globals[MAX_GLOBALS];
   int nconds;           Cond  *condlist[MAX_CONDS];       //   sorted list
   int ngates;           Gate  *gatelist[MAX_GATES];//(inactive at end)           47384
   int nusedvar;         Sortvar *usedvars[MAX_SORT_VARS];
   int nuser;            Histogram *user_histos[MAX_HISTOGRAMS];
   int nhistos;          Histogram *histo_list[MAX_HISTOGRAMS];
   int nsortvar;         Sortvar varlist[MAX_SORT_VARS];                     // 33921336
   Cond cond_array[MAX_GATES];  Gate gate_array[MAX_GATES]; // unsorted
   Global global_array[MAX_GLOBALS];                                         // 34115896
   Histogram histo_array[MAX_HISTOGRAMS];
   Cal_coeff calib_array[MAX_CALIB];  int odb_daqsize;
} Config;

extern Config *configs[MAX_CONFIGS]; // start unallocated

extern int init_default_histos(Config *cfg, Sort_status *arg);
extern int remove_histo(Config *cfg, char *name);
extern int add_histo(Config *cfg, char *path, char *name, char *title, int xbins, char *xvarname, int xmin, int max, int ybins, char *yvarname, int ymin, int ymax);
extern Histogram *find_histo(Config *cfg, char *name);
extern int read_histo_data(Histogram *histo, FILE *fp);

extern Config *add_config(char *name);
extern int remove_config(Config *cfg);
extern int next_condname(Config *cfg);
extern int queue_sum_histos(Config *cfg, int num, char url_args[][STRING_LEN], int fd);
extern int set_calibration(Config *cfg, int num, char url_args[][STRING_LEN], int fd);
extern int set_pileup_correction(Config *cfg, int num, char url_args[][STRING_LEN], int fd);
extern int set_crosstalk_correction(Config *cfg, int num, char url_args[][STRING_LEN], int fd);

extern int sum_histos(Config *cfg, Sortfile *sort);
/////////////////////////////////////////////////////////////////////////
/////////////////////          Gains         ////////////////////////////
/////////////////////////////////////////////////////////////////////////
extern int edit_calibration(Config *cfg, char *name, float offset, float gain, float quad, float pileupk1[7], float pileupk2[7], float pileupE1[7], float ct0[7], float ct1[7], float ct2[7], int address, int type, int overwrite);

/////////////////////////////////////////////////////////////////////////
/////////////////////       Variables        ////////////////////////////
/////////////////////////////////////////////////////////////////////////

/////// sort variables are currently predefined and fixed at compilation time
/////// maybe later will add ability to add/remove them during runtime
//extern int add_variable(char *name);
//extern int remove_variable(char *name);
extern Sortvar *find_sortvar(Config *cfg, char *name);

extern int add_global(Config *cfg, char *name, int value, int val2);
extern int remove_global(Config *cfg, char *name);

/////////////////////////////////////////////////////////////////////////
/////////////////////          Gates         ////////////////////////////
/////////////////////////////////////////////////////////////////////////

#define GATEOP_LT    0
#define GATEOP_LE    1
#define GATEOP_GT    2
#define GATEOP_GE    3
#define GATEOP_EQ    4
#define GATEOP_RA    5
//#define GATEOP(OP) (         ((OP)==GATEOP_EQ) ? "=" :      \
//   ((OP)==GATEOP_LT) ? "<" : ((OP)==GATEOP_LE) ? "<=" :       \
//   ((OP)==GATEOP_GT) ? ">" : ((OP)==GATEOP_GE) ? ">=" : "UNK" )
#define GATEOP(OP) (                                           \
   ((OP)==GATEOP_RA) ? "RA" : ((OP)==GATEOP_EQ) ? "EQ" :       \
   ((OP)==GATEOP_LT) ? "LT" : ((OP)==GATEOP_LE) ? "LE" :       \
   ((OP)==GATEOP_GT) ? "GT" : ((OP)==GATEOP_GE) ? "GE" : "UNK" )

extern int add_cond(Config *cfg, char *name, char *var, char *op, int value);
extern int remove_cond(Config *cfg, char *name);
extern int apply_gate(Config *cfg, char *histoname, char *gatename);
extern int unapply_gate(Config *cfg, char *histoname, char *gatename);
extern int add_gate(Config *cfg, char *name);
extern int remove_gate(Config *cfg, char *name);
extern int add_cond_to_gate(Config *cfg, char *gatename, char *condname);

/////////////////////////////////////////////////////////////////////////
/////////////////          Config Files         /////////////////////////
/////////////////////////////////////////////////////////////////////////

extern int init_config();
extern int init_user_config(Config *cfg);
extern int init_default_config(Config *cfg);
extern int write_config(Config *cfg, FILE *fp);
extern int copy_config(Config *src, Config *dst);
extern int clear_config(Config *cfg);
extern int clear_calibrations(Config *cfg);
extern int delete_config(Config *cfg);
extern int load_config(Config *cfg, char *filename, char *buffer);
extern int save_config(Config *cfg, char *filename, int overwrite);

extern int set_directory(Config *cfg, char *name, char *path);
extern int set_midas_param(Config *cfg, char *name, char *value);

#define HPGeE        0 // HPGE
#define HPGeA        1 //
#define HPGeT        2 //
#define HPGeTS       3 //
#define HPGePH       4 //
#define HPGeC        5 //
#define HPGeCL       6 //
#define HPGePU       7 //
#define HPGeIP       8 //
#define HPGeDT       9 //
#define HPGeEU      10 // #
#define HPGeAU      11 //
#define GRGTHETA    12 //
#define GRGPHI      13 //
#define CLVTHETA    14 //
#define CLVPHI      15 //
#define SEPE        16 // SCEPTAR
#define SEPT        17 //
#define SEPTS       18 //
#define EPPH        19 //
#define SEPPU       20 //
#define SEPTHETA    21 //
#define SEPPH       22 //
#define SEPNUM      23 //
#define PACE        24 // PACES
#define PACT        25 //
#define PACTS       26 //
#define PACPH       27 //
#define PACPU       28 //
#define PACTHETA    29 //
#define PACPHI      30 //
#define PACNUM      31 //
#define LBLE        32 // LaBr3
#define LBLT        33 //
#define LBLTS       34 //
#define LBLPH       35 //
#define LBLPU       36 //
#define LBLTHETA    37 //
#define LBLPHI      38 //
#define LBLNUM      39 //
#define LBT         40 // TACs
#define LBTT        41 //
#define LBTTS       42 //
#define LBTPH       43 //
#define LBTPU       44 //
#define LBTNUM      45 //
#define GRSE        46 // Clover BGO
#define GRST        47 //
#define GRSTS       48 //
#define GRSPH       49 //
#define GRSPU       50 //
#define GRSNUM      51 //
#define GRSPOS      52 //
#define GRSTYPE     53 //
#define LBSE        54 // Ancillary BGO
#define LBST        55 //
#define LBSTS       56 //
#define LBSPH       57 //
#define LBSPU       58 //
#define LBSNUM      59 //
#define LBSPOS      60 //
#define MIDAS_Time  61 // Time Differences
#define TD_GRG_GRG  62 //
#define TD_SEP_SEP  63 //
#define TD_PAC_PAC  64 //
#define TD_LBL_LBL  65 //
#define TSD_GRG_GRG 66 //
#define TSD_SEP_SEP 67 //
#define TSD_PAC_PAC 68 //
#define TSD_LBL_LBL 69 //
#define TD_GRG_SEP  70 //
#define TSD_GRG_SEP 71 //
#define TD_GRG_PAC  72 //
#define TSD_GRG_PAC 73 //
#define TD_SEP_PAC  74 //
#define TSD_SEP_PAC 75 //
#define TD_GRG_LBL  76 //
#define TSD_GRG_LBL 77 //
#define TD_SEP_LBL  78 //
#define TSD_SEP_LBL 79 //
#define ANG_GRG_GRG 80 // Angular Differences [HpGeDistDependant]
#define ANG_CLV_CLV 81 //
#define ANG_SEP_SEP 82 //
#define ANG_PAC_PAC 83 //
#define ANG_LBL_LBL 84 //
#define ANG_GRG_SEP 85 //
#define ANG_GRG_PAC 86 //
#define ANG_GRG_LBL 87 //
#define ANG_PAC_LBL 88 //
#define ANG_SEP_LBL 89 //
#define PPG_NUM     90 // Cycle Timing (PPG events)
#define PPG_TIME    91 //
#define PPG_PAT     92 //

#endif
