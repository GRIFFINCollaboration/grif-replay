#define MAX_SAMPLE_LEN  4096
#define ENERGY_BINS    65536 /* 65536 131072 262144 */
#define NUM_CHAN        4096
#define MAX_SCALAR_LEN   256
#define MAX_COINC_EVENTS 1024

// **************************************************************************
// MANY VALUES IN THIS STRUCTURE HAVE TO BE ACCESSED USING A HARDCODED OFFSET
// FROM THE START OF THE STRUCTURE - DO NOT CHANGE THE ORDERING OR
// INSERT NEW VALUES, WITHOUT ALSO ADJUSTING THE OFFSETS IN USER_SORT.C
// **************************************************************************
typedef struct griffin_fragment_struct { // was 74 bytes, now ?
   int     address;  int chan;          int     dtype;    int array_posn; // 0
   int        ecal;  int e2cal;         int     e3cal;    int    e4cal;   // 4
   int      energy;  int q;             int     e_bad;    int    integ;   // 8
   int     energy2;  int q2;            int    e2_bad;    int   integ2;   //12
   int     energy3;  int q3;            int    e3_bad;    int   integ3;   //16
   int     energy4;  int q4;            int    e4_bad;    int   integ4;   //20
   int         cfd;  int cc_long;       int  cc_short;    int     nhit;   //24
   int    trig_req;  int trig_acc;      int     pileup;   int suppress;   //28
   int   deadtime;   int bank_start;    int    ts_int;    int    esum_i;  //32
   int   master_id;  int master_pattern;int  esum_supr;   int    e_supr;  //36
   int     crystal;  int fold;          int     subsys;   int    dummy4;  //40
   int      angle1;  int angle2;        int     angle3;   int    angle4;  //44
 //int      dummy1;  int dummy2;        int     dummy3;   int    dummy4;  //##
   int      net_id;  int trigger_num;   long timestamp;  long ts;
   int  wf_present;  int waveform_length; int file_id;
   int scl_present;  int scalar_length;  float    esum;
} Grif_event;

// Q* are the original Q-Sum from the digitizer
// Energy are the Q-Sums divided by integ-times
// Ecal is calibrated energy#1

// ppg_pattern is now wave_expected, num_pileup

// midas-timestamp should be redundant and equal to BOR time+timestamp
//   can just check this in midas part

// Global counters for finding the Ge-RCMP coincidences
int presort_ge_events_passed;
int presort_rcmp_events_passed;
int presort_rcmp_fb_events_built;
int singles_ge_events;
int singles_rcmp_events;
int singles_rcmp_fb_events;
int coinc_ge_rcmp_events;
int coinc_rcmp_ge_events;
int coinc_ge_rcmp_fb_events;
int coinc_rcmp_ge_fb_events;
