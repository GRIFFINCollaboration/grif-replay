#define MAX_SAMPLE_LEN  4096
#define ENERGY_BINS    65536 /* 65536 131072 262144 */
#define NUM_CHAN        4096
#define MAX_SCALAR_LEN   256
#define PTR_BUFSIZE     1024

typedef unsigned short uint_16;
// **************************************************************************
// MANY VALUES IN THIS STRUCTURE HAVE TO BE ACCESSED USING A HARDCODED OFFSET
// FROM THE START OF THE STRUCTURE - DO NOT CHANGE THE ORDERING OR
// INSERT NEW VALUES, WITHOUT ALSO ADJUSTING THE OFFSETS IN USER_SORT.C
// **************************************************************************
typedef struct griffin_fragment_struct { // was 74 bytes, now ?            //OFFSET
   long        ts;   uint_16    address; short  deadtime;                  //0
   char      dtype;  char    array_posn; char       nhit;   char  pileup;  //3
   int          q1;  int         integ1; int          q2;   int   integ2;  //4
   int          q3;  int         integ3; int          q4;   int   integ4;  //8
   int         cfd;  int       trig_req; int    trig_acc;   int   net_id;  //12
   int   master_id;  int master_pattern; int         psd;   int cc_short;  //16
// items below are derived from items above ...
   float      ecal;  int           chan; int      subsys;   int suppress;  //20
   float      esum;  int   multiplicity; int     delta_t;   int alt_chan;  //24
   float  alt_ecal;  int            tof; int    pu_class;                  //28
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

extern int process_event(Grif_event *ptr, int slot);
extern int insert_presort_win(Grif_event *ptr, int slot);
extern int insert_sort_win(Grif_event *ptr, int slot);
extern int pre_sort_enter(int start_idx, int frag_idx);
extern int pre_sort_exit(int frag_idx, int end_idx);

// User sort function declarations
extern int calc_coincvars(Grif_event *ptr1, Grif_event *ptr2);
extern int user_removefrom_window(int win_strt, int new_frag);
extern int test_gates(Grif_event *ptr, Grif_event *alt);
