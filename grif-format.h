#include "grif-replay.h" // STRING_LEN

#define MAX_SAMPLE_LEN  4096
#define ENERGY_BINS    65536 /* 65536 131072 262144 */
#define NUM_CHAN        4096
#define MAX_SCALAR_LEN   256
#define PTR_BUFSIZE     1024

typedef unsigned short uint_16;
typedef struct griffin_fragment_struct { // was 74 bytes, now ?            
   long        ts;   uint_16    address; short  deadtime;                  
   char      dtype;  char    array_posn; char       nhit;   char  pileup;  
   int          q1;  int         integ1; int          q2;   int   integ2;  
   int          q3;  int         integ3; int          q4;   int   integ4;  
   int         cfd;  int       trig_req; int    trig_acc;   int   net_id;  
   int   master_id;  int master_pattern; int         psd;   int cc_short;  
// items below are derived from items above ...
//    esum:is addback energy, alt_chan is addback-otherCrystal
   float      ecal;  int           chan; int      subsys;   int suppress;  
   float      esum;  int   multiplicity; int     delta_t;   int alt_chan;  
   float  alt_ecal;  int            tof; int    pu_class;   int alt2_chan;
   float alt2_ecal;
} Grif_event;

// NOTES on Grif_event ...
//    midas data only contains the following ...
//       ts, address; deadtime, dtype, array_posn. nhit, pileup, q1-4, integ1-4
//       cfd, trig_req, trig_acc. net_id, master_id, master_pattern
//    The other itmes are derived from the above ...
//       chan, subsys, ecal, crystal, fold, psd, suppress, bank_start,
//       esum_i, esum_supr, e_supr, trigger_num, multiplcity
//
// NOTES on types ...
//    address is 16bits 0-FFFF                               => unsigned 16bits
//    chan is never over 32k, and want to allow -1 for unset => signed 16bits
//
// Q* are the original Q-Sum from the digitizer
// Energy are the Q-Sums divided by integ-times
// Ecal is calibrated energy#1

// ppg_pattern is now wave_expected, num_pileup

// midas-timestamp should be redundant and equal to BOR time+timestamp
//   can just check this in midas part

//#######################################################################
//########         Subsystem and Detector definitions          ##########
//#######################################################################

// do not alter order without also changing subsys_e_vs_e, subsys_dt
//                                              in default_sort
#define MAX_SUBSYS       24
#define SUBSYS_HPGE_A     0
#define SUBSYS_PACES      1
#define SUBSYS_LABR_L     2
#define SUBSYS_RCMP       3
#define SUBSYS_ARIES_A    4 // GRIF16
#define SUBSYS_ZDS_A      5 // GRIF16
#define SUBSYS_TAC_LABR   6
#define SUBSYS_LABR_BGO   7
#define SUBSYS_BGO        8
#define SUBSYS_SCEPTAR    9
#define SUBSYS_DESCANT   10
#define SUBSYS_DESWALL   11
#define SUBSYS_DSG       12
#define SUBSYS_QED_STRIP 13
#define SUBSYS_IGNORE    14
#define SUBSYS_HPGE_B    16
#define SUBSYS_ARIES_B   17 // CAEN
#define SUBSYS_ZDS_B     18 // CAEN
#define SUBSYS_TAC_ZDS   19
#define SUBSYS_TAC_ART   20
#define SUBSYS_COMPTON   21
#define SUBSYS_DCOMPTON  15
#define SUBSYS_QED_PIXEL 22
#define SUBSYS_UNKNOWN   23

extern char subsys_handle[MAX_SUBSYS][8];
extern char subsys_name[MAX_SUBSYS][STRING_LEN];

// #####################################################################

#define N_CLOVER 16
#define N_HPGE 64
#define N_BGO 320
#define N_ARIES 76
#define N_LABR 8
#define N_TACS 12
#define N_RCMP_POS 6
#define N_RCMP_STRIPS 32
#define N_QED_POS 6
#define N_QED_STRIPS 32
#define N_DES_WALL 60

extern int process_event(Grif_event *ptr, int slot);
extern int insert_presort_win(Grif_event *ptr, int slot);
extern int insert_sort_win(Grif_event *ptr, int slot);
extern int pre_sort_enter(int start_idx, int frag_idx);
extern int pre_sort_exit(int frag_idx, int end_idx);
extern int pre_sort_triples(int frag_idx, int end_idx);

// User sort function declarations
extern int calc_coincvars(Grif_event *ptr1, Grif_event *ptr2);
extern int user_removefrom_window(int win_strt, int new_frag);
extern int test_gates(Grif_event *ptr, Grif_event *alt);
