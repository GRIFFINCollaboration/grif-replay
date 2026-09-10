//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#####        HISTOGRAM DEFINITION AND INITIALIZATION              #####
//#######################################################################

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "config.h"
#include "grif-format.h"
#include "histogram.h"
#include "default_sort_griffin.h"


int DEBUG_OUTPUT=0; // 0 for off, 1 for on

char subsys_handle[MAX_SUBSYS][8] = {
  "GRGA", "PAC",  "LBL",  "RCS",
  "ARTA", "ZDSA", "LBT",  "LBS",
  "BGO",  "SEP",  "DSC",  "DSW",
  "DSG",  "QEDs", "XXX",  "DCSA",
  "GRGB", "ARTB", "ZDSB", "TACZ", // secondary names start after #16
  "TACA", "CS",   "QED",  "DCSB",
  "XXX3", "XXX4", "XXX5", "UNK"
};
char subsys_name[MAX_SUBSYS][STRING_LEN] = {
  "Griffin",   "PACES",    "LaBrX",   "RCMP",     //  0- 3
  "ARIES",     "ZDSA",     "TAC_LBL", "LaBrS",    //  4- 7
  "BGO",       "Sceptar",  "Descant", "DES_WALL", //  8-11
  "Des_Ancil", "QEDs",     "Ignore",  "DCSA",     // 12-15
  "Grif_B",    "ARS_B",    "ZDS_B",   "TAC_ZDS",  // 16-19
  "TAC_ART",   "CS",       "QED",     "DCSB",     // 20-23
  "Ignore3",   "Ignore4",  "Ignore5", "Unknown"   // 24-27
}; // final entry will be used if not found - make sure it is not empty

int          odb_daqsize;// number of daq channels currently defined in the odb
int         subsys_table[MAX_DAQSIZE];
int        crystal_table[MAX_DAQSIZE]; // Ge/BGO have 4 "crystals" per clover
int        element_table[MAX_DAQSIZE]; // BGO have 5 elements per crystal
int       polarity_table[MAX_DAQSIZE]; // 1 is negative, 0 is positive, -1 is unset
short       address_chan[MAX_ADDRESS];
short  addr_table[MAX_DAQSIZE]; short   *addrs = addr_table;
char    chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
int   dtype_table[MAX_DAQSIZE]; int    *dtypes = dtype_table;
float  gain_table[MAX_DAQSIZE]; float   *gains = gain_table;
float  offs_table[MAX_DAQSIZE]; float *offsets = offs_table;
float  quad_table[MAX_DAQSIZE]; float   *quads = quad_table;
float  pileupk1[MAX_DAQSIZE][7];
float  pileupk2[MAX_DAQSIZE][7];
float  pileupE1[MAX_DAQSIZE][7];
float  crosstalk[MAX_DAQSIZE][3][16];
short *chan_address = addr_table;
int subsys_initialized[MAX_SUBSYS];
int subsys_deadtime_count[MAX_SUBSYS];
int subsys_prg_ddtm[MAX_SUBSYS];
int previous_trig_acc[MAX_DAQSIZE];

int presort_window_width= 1940;  // 19.4us needed for all crosstalk corrections. 5us needed for pileup corrections.
int sort_window_width   = 200; //  2us - MAXIMUM (indiv. gates can be smaller)
int ct_index[4][4] = {{-1,0,1,2},{0,-1,1,2},{0,1,-1,2},{0,1,2,-1}}; // crosstalk index used in presort functions

long ppg_cycles_binning_factor = 10000000;      // Default of 10,000,000 converts 10ns to 100 millisecond binning
int ppg_cycles_gamma_gate_min = 1800;           // 26Na, S1140
int ppg_cycles_gamma_gate_max = 1820;           // 26Na, S1140

// Definitions for the Current cycle of this run
// These are derived from the ODB settings at BOR
long ppg_cycle_duration;             // Length of one cycle in timestamp units
char ppg_cycle_name[128];            // Name of this current cycle
long ppg_cycle_length;               // Number of patterns/durations for this current cycle
long ppg_cycle_pattern_duration[16]; // Length of each pattern in timestamp units
int  ppg_cycle_pattern_code[16];     // Index of each pattern for use with the ppg_patterns array
int  ppg_cycles_active;              // Cycles active or made inactive if set to Source/constant beam-on etc

// These variables are updated in pre_sort_enter at each PPG pattern change
long ppg_last_ptr_ts;     // Previous event timestamp. Avoids rare bug where single events are out of order.
int ppg_current_pattern;  // Index of the current PPG cycle pattern for use with the ppg_patterns array
int ppg_cycle_number;     // Current cycle number. Cycles counted from zero at beginning of run
long ppg_cycle_start;     // Timestamp of the start of the current cycle
long ppg_cycle_end;       // Timestamp of the end of the current cycle
int ppg_cycle_step;       // Current pattern number within this cycle. Patterns counted from zero at beginning of cycle
long ppg_pattern_start;   // Timestamp of the start of the current pattern
long ppg_pattern_end;     // Timestamp of the end of the current pattern
long ppg_bin_end;         // Timestamp of the end of the current bin (used for deadtime histogram)

//#######################################################################
//########          HISTOGRAM TITLES AND HANDLES               ##########
//#######################################################################

// Cycles & PPG handles
int ppg_patterns[N_PPG_PATTERNS]={ 0xC008,0xC002,0xC001,0xC004,0xC0F0,0xC009,0xC00A };
char ppg_handles[N_PPG_PATTERNS][32]={ "0xC008","0xC002","0xC001","0xC004","0xC0F0","0xC009","0xC00A" };
char ppg_names[N_PPG_PATTERNS][32]={
  "Move_Tape", "Background", "Beam_on_Implant", "Beam_off_Decay", "Source_data", "Continuous_Tape_Beam_on", "Continuous_Tape_Background"
};

char ge_cycle_code_titles[N_PPG_PATTERNS][HANDLE_LENGTH]={ "0xC008_Move_Tape","0xC002_Background","0xC001_Beam_on_Implant","0xC004_Beam_off_Decay","0xC0F0_Source_data","0xC009_Continuous_Tape_Beam_on", "0xC00A_Continuous_Tape_Background" };
char gg_cycle_code_titles[N_PPG_PATTERNS][HANDLE_LENGTH]={ "0xC008_GG_Move_Tape","0xC002_GG_Background","0xC001_GG_Beam_on_Implant","0xC004_GG_Beam_off_Decay","0xC0F0_GG_Source_data","0xC009_GG_Continuous_Tape_Beam_on", "0xC00A_GG_Continuous_Tape_Background" };
// Hitpatterns
char hit_handles[N_HITPAT][32]={ "q_hit","e_hit","t_hit","w_hit","r_hit", "s_hit", "d_hit" };
char   hit_names[N_HITPAT][32]={
  "Pulse_Height", "Energy", "Time", "Waveform", "Rate", "Subsys", "DetType",
};

static char ge_pu_class_handles[N_PU_CLASSES][HANDLE_LENGTH]={
  "PU0",                          //  0    = Pu=0 error
  "PU1_NHIT1", "PU1_NHIT1_error", //  1, 2 = Single Hit, single hit error
  "PU1_NHIT2", "PU2_NHIT1",       //  3, 4 = 2Hit pile-up, corrected Hit1, corrected Hit2
  "NHIT2_lateA1", "NHIT2_lateA2",   //  5, 6 = 2Hit pile-up,  separate Hit1,  separate Hit2
  "NHIT2_lateB1", "NHIT2_lateB2",   //  7, 8 = 2Hit pile-up,  separate Hit1,  separate Hit2
  "NHIT2_error",                  //  9,   = 2Hit pile-up event error, most likely q1 or q2 is zero
  "PU1_NHIT3", "PU2_NHIT2", "PU3_NHIT1", "NHIT3_error", // 10-13, Three pile-up events, Hit1, Hit2, Hit3, error
  "Other_PU" // 14
};
static char ge_pu_class_sum_titles[N_PU_CLASSES][HANDLE_LENGTH]={
  "PU_zero",
  "single_hit", "error_single_hit", // Single Hit, single hit error
  "2Hit_PU_Type_A_1stHit", "2Hit_PU_Type_A_2ndHit", "2Hit_PU_Type_B_1stHit", "2Hit_PU_Type_B_2ndHit", "2Hit_PU_Type_C_1stHit", "2Hit_PU_Type_C_2ndHit", "2Hit_PU_error", // Two pile-up events, corrected Hit1, corrected Hit2, separate Hit1, separate Hit2, error
  "3Hit_PU_1stHit", "3Hit_PU_2ndHit", "3Hit_PU_3rdHit", "3Hit_PU_error", // 3 pile-up events, Hit1, Hit2, Hit3, error
  "Other_PU"
};
static char ge_pu_class_2d_titles[N_PU_CLASSES][HANDLE_LENGTH]={
  "E_vs_k_PU_zero",
  "E_vs_k_single_hit", "E_vs_k_error_single_hit", // Single Hit, single hit error
  "E_vs_k_2Hit_PU_Type_A_1stHit", "E_vs_k_2Hit_PU_Type_A_2ndHit", "E_vs_k_2Hit_PU_Type_B_1stHit", "E_vs_k_2Hit_PU_Type_B_2ndHit", "E_vs_k_2Hit_PU_Type_C_1stHit", "E_vs_k_2Hit_PU_Type_C_2ndHit", "E_vs_k_2Hit_PU_error", // Two pile-up events, corrected Hit1, corrected Hit2, separate Hit1, separate Hit2, error
  "E_vs_k_3Hit_PU_1stHit", "E_vs_k_3Hit_PU_2ndHit", "E_vs_k_3Hit_PU_3rdHit", "E_vs_k_3Hit_PU_error", // 3 pile-up events, Hit1, Hit2, Hit3, error
  "E_vs_k_Other_PU"
};
static char ge_hem_handles[3][HANDLE_LENGTH] = {"ge_downstream","ge_upstream","ge_nope"};

char qed_psd_handles[N_QED_POS][HANDLE_LENGTH] = {"QED01_E_vs_psd","QED02_E_vs_psd","QED03_E_vs_psd","QED04_E_vs_psd","QED05_E_vs_psd","QED06_E_vs_psd"};
char qed_strips_handles[N_QED_POS][HANDLE_LENGTH]={"QED01_E_strips", "QED02_E_strips", "QED03_E_strips", "QED04_E_strips", "QED05_E_strips", "QED06_E_strips"};
char qed_hit_handles[N_QED_POS][HANDLE_LENGTH]={"QED01_PN_hit", "QED02_PN_hit", "QED03_PN_hit", "QED04_PN_hit", "QED05_PN_hit", "QED06_PN_hit"};
char qed_fb_handles[N_QED_POS][HANDLE_LENGTH]={"QED01_Front_Back", "QED02_Front_Back", "QED03_Front_Back", "QED04_Front_Back", "QED05_Front_Back", "QED06_Front_Back"};
char qed_totE_theta_handles[N_QED_POS][HANDLE_LENGTH]={"QED01_totalE_vs_theta", "QED02_totalE_vs_theta", "QED03_totalE_vs_theta", "QED04_totalE_vs_theta", "QED05_totalE_vs_theta", "QED06_totalE_vs_theta"};
char qed_E_theta_handles[N_QED_POS][HANDLE_LENGTH]={"QED01_E_vs_theta", "QED02_E_vs_theta", "QED03_E_vs_theta", "QED04_E_vs_theta", "QED05_E_vs_theta", "QED06_E_vs_theta"};
char qed_geE_theta_handles[N_QED_POS][HANDLE_LENGTH]={"QED01_geE_vs_theta", "QED02_geE_vs_theta", "QED03_geE_vs_theta", "QED04_geE_vs_theta", "QED05_geE_vs_theta", "QED06_geE_vs_theta"};
char qed_E_cycle_handles[N_QED_POS][HANDLE_LENGTH]={"cycle_vs_QED01_Energy", "cycle_vs_QED02_Energy", "cycle_vs_QED03_Energy", "cycle_vs_QED04_Energy", "cycle_vs_QED05_Energy", "cycle_vs_QED06_Energy", };
char qedx_dcs_omega_dt_handles[N_QED_POS][HANDLE_LENGTH] = {"QED01_DCS_omega_vs_dt","QED02_DCS_omega_vs_dt","QED03_DCS_omega_vs_dt","QED04_DCS_omega_vs_dt","QED05_DCS_omega_vs_dt","QED06_DCS_omega_vs_dt"};

char qed_geE_theta_clov_handles[N_CLOVER][HANDLE_LENGTH] = {
  "QED_Clover01E_vs_theta","QED_Clover02E_vs_theta","QED_Clover03E_vs_theta","QED_Clover04E_vs_theta",
  "QED_Clover05E_vs_theta","QED_Clover06E_vs_theta","QED_Clover07E_vs_theta","QED_Clover08E_vs_theta",
  "QED_Clover09E_vs_theta","QED_Clover10E_vs_theta","QED_Clover11E_vs_theta","QED_Clover12E_vs_theta",
  "QED_Clover13E_vs_theta","QED_Clover14E_vs_theta","QED_Clover15E_vs_theta","QED_Clover16E_vs_theta"
};

char qed_geE_theta_clov_t_handles[N_CLOVER][HANDLE_LENGTH] = {
  "QED_Clover01E_vs_theta_totEgated","QED_Clover02E_vs_theta_totEgated","QED_Clover03E_vs_theta_totEgated","QED_Clover04E_vs_theta_totEgated",
  "QED_Clover05E_vs_theta_totEgated","QED_Clover06E_vs_theta_totEgated","QED_Clover07E_vs_theta_totEgated","QED_Clover08E_vs_theta_totEgated",
  "QED_Clover09E_vs_theta_totEgated","QED_Clover10E_vs_theta_totEgated","QED_Clover11E_vs_theta_totEgated","QED_Clover12E_vs_theta_totEgated",
  "QED_Clover13E_vs_theta_totEgated","QED_Clover14E_vs_theta_totEgated","QED_Clover15E_vs_theta_totEgated","QED_Clover16E_vs_theta_totEgated"
};

char qedp_ge_theta_handles[N_QED_POS*N_QED_STRIPS][HANDLE_LENGTH]={
  "QED1P00_E_vs_theta", "QED1P01_E_vs_theta", "QED1P02_E_vs_theta", "QED1P03_E_vs_theta", "QED1P04_E_vs_theta", "QED1P05_E_vs_theta",
  "QED1P06_E_vs_theta", "QED1P07_E_vs_theta", "QED1P08_E_vs_theta", "QED1P09_E_vs_theta", "QED1P10_E_vs_theta", "QED1P11_E_vs_theta",
  "QED1P12_E_vs_theta", "QED1P13_E_vs_theta", "QED1P14_E_vs_theta", "QED1P15_E_vs_theta", "QED1P16_E_vs_theta", "QED1P17_E_vs_theta",
  "QED1P18_E_vs_theta", "QED1P19_E_vs_theta", "QED1P20_E_vs_theta", "QED1P21_E_vs_theta", "QED1P22_E_vs_theta", "QED1P23_E_vs_theta",
  "QED1P24_E_vs_theta", "QED1P25_E_vs_theta", "QED1P26_E_vs_theta", "QED1P27_E_vs_theta", "QED1P28_E_vs_theta", "QED1P29_E_vs_theta",
  "QED1P30_E_vs_theta", "QED1P31_E_vs_theta",

  "QED2P00_E_vs_theta", "QED2P01_E_vs_theta", "QED2P02_E_vs_theta", "QED2P03_E_vs_theta", "QED2P04_E_vs_theta", "QED2P05_E_vs_theta",
  "QED2P06_E_vs_theta", "QED2P07_E_vs_theta", "QED2P08_E_vs_theta", "QED2P09_E_vs_theta", "QED2P10_E_vs_theta", "QED2P11_E_vs_theta",
  "QED2P12_E_vs_theta", "QED2P13_E_vs_theta", "QED2P14_E_vs_theta", "QED2P15_E_vs_theta", "QED2P16_E_vs_theta", "QED2P17_E_vs_theta",
  "QED2P18_E_vs_theta", "QED2P19_E_vs_theta", "QED2P20_E_vs_theta", "QED2P21_E_vs_theta", "QED2P22_E_vs_theta", "QED2P23_E_vs_theta",
  "QED2P24_E_vs_theta", "QED2P25_E_vs_theta", "QED2P26_E_vs_theta", "QED2P27_E_vs_theta", "QED2P28_E_vs_theta", "QED2P29_E_vs_theta",
  "QED2P30_E_vs_theta", "QED2P31_E_vs_theta",

  "QED3P00_E_vs_theta", "QED3P01_E_vs_theta", "QED3P02_E_vs_theta", "QED3P03_E_vs_theta", "QED3P04_E_vs_theta", "QED3P05_E_vs_theta",
  "QED3P06_E_vs_theta", "QED3P07_E_vs_theta", "QED3P08_E_vs_theta", "QED3P09_E_vs_theta", "QED3P10_E_vs_theta", "QED3P11_E_vs_theta",
  "QED3P12_E_vs_theta", "QED3P13_E_vs_theta", "QED3P14_E_vs_theta", "QED3P15_E_vs_theta", "QED3P16_E_vs_theta", "QED3P17_E_vs_theta",
  "QED3P18_E_vs_theta", "QED3P19_E_vs_theta", "QED3P20_E_vs_theta", "QED3P21_E_vs_theta", "QED3P22_E_vs_theta", "QED3P23_E_vs_theta",
  "QED3P24_E_vs_theta", "QED3P25_E_vs_theta", "QED3P26_E_vs_theta", "QED3P27_E_vs_theta", "QED3P28_E_vs_theta", "QED3P29_E_vs_theta",
  "QED3P30_E_vs_theta", "QED3P31_E_vs_theta",

  "QED4P00_E_vs_theta", "QED4P01_E_vs_theta", "QED4P02_E_vs_theta", "QED4P03_E_vs_theta", "QED4P04_E_vs_theta", "QED4P05_E_vs_theta",
  "QED4P06_E_vs_theta", "QED4P07_E_vs_theta", "QED4P08_E_vs_theta", "QED4P09_E_vs_theta", "QED4P10_E_vs_theta", "QED4P11_E_vs_theta",
  "QED4P12_E_vs_theta", "QED4P13_E_vs_theta", "QED4P14_E_vs_theta", "QED4P15_E_vs_theta", "QED4P16_E_vs_theta", "QED4P17_E_vs_theta",
  "QED4P18_E_vs_theta", "QED4P19_E_vs_theta", "QED4P20_E_vs_theta", "QED4P21_E_vs_theta", "QED4P22_E_vs_theta", "QED4P23_E_vs_theta",
  "QED4P24_E_vs_theta", "QED4P25_E_vs_theta", "QED4P26_E_vs_theta", "QED4P27_E_vs_theta", "QED4P28_E_vs_theta", "QED4P29_E_vs_theta",
  "QED4P30_E_vs_theta", "QED4P31_E_vs_theta",

  "QED5P00_E_vs_theta", "QED5P01_E_vs_theta", "QED5P02_E_vs_theta", "QED5P03_E_vs_theta", "QED5P04_E_vs_theta", "QED5P05_E_vs_theta",
  "QED5P06_E_vs_theta", "QED5P07_E_vs_theta", "QED5P08_E_vs_theta", "QED5P09_E_vs_theta", "QED5P10_E_vs_theta", "QED5P11_E_vs_theta",
  "QED5P12_E_vs_theta", "QED5P13_E_vs_theta", "QED5P14_E_vs_theta", "QED5P15_E_vs_theta", "QED5P16_E_vs_theta", "QED5P17_E_vs_theta",
  "QED5P18_E_vs_theta", "QED5P19_E_vs_theta", "QED5P20_E_vs_theta", "QED5P21_E_vs_theta", "QED5P22_E_vs_theta", "QED5P23_E_vs_theta",
  "QED5P24_E_vs_theta", "QED5P25_E_vs_theta", "QED5P26_E_vs_theta", "QED5P27_E_vs_theta", "QED5P28_E_vs_theta", "QED5P29_E_vs_theta",
  "QED5P30_E_vs_theta", "QED5P31_E_vs_theta",

  "QED6P00_E_vs_theta", "QED6P01_E_vs_theta", "QED6P02_E_vs_theta", "QED6P03_E_vs_theta", "QED6P04_E_vs_theta", "QED6P05_E_vs_theta",
  "QED6P06_E_vs_theta", "QED6P07_E_vs_theta", "QED6P08_E_vs_theta", "QED6P09_E_vs_theta", "QED6P10_E_vs_theta", "QED6P11_E_vs_theta",
  "QED6P12_E_vs_theta", "QED6P13_E_vs_theta", "QED6P14_E_vs_theta", "QED6P15_E_vs_theta", "QED6P16_E_vs_theta", "QED6P17_E_vs_theta",
  "QED6P18_E_vs_theta", "QED6P19_E_vs_theta", "QED6P20_E_vs_theta", "QED6P21_E_vs_theta", "QED6P22_E_vs_theta", "QED6P23_E_vs_theta",
  "QED6P24_E_vs_theta", "QED6P25_E_vs_theta", "QED6P26_E_vs_theta", "QED6P27_E_vs_theta", "QED6P28_E_vs_theta", "QED6P29_E_vs_theta",
  "QED6P30_E_vs_theta", "QED6P31_E_vs_theta",
};

char qedn_ge_theta_handles[N_QED_POS*N_QED_STRIPS][HANDLE_LENGTH]={
  "QED1N00_E_vs_theta", "QED1N01_E_vs_theta", "QED1N02_E_vs_theta", "QED1N03_E_vs_theta", "QED1N04_E_vs_theta", "QED1N05_E_vs_theta",
  "QED1N06_E_vs_theta", "QED1N07_E_vs_theta", "QED1N08_E_vs_theta", "QED1N09_E_vs_theta", "QED1N10_E_vs_theta", "QED1N11_E_vs_theta",
  "QED1N12_E_vs_theta", "QED1N13_E_vs_theta", "QED1N14_E_vs_theta", "QED1N15_E_vs_theta", "QED1N16_E_vs_theta", "QED1N17_E_vs_theta",
  "QED1N18_E_vs_theta", "QED1N19_E_vs_theta", "QED1N20_E_vs_theta", "QED1N21_E_vs_theta", "QED1N22_E_vs_theta", "QED1N23_E_vs_theta",
  "QED1N24_E_vs_theta", "QED1N25_E_vs_theta", "QED1N26_E_vs_theta", "QED1N27_E_vs_theta", "QED1N28_E_vs_theta", "QED1N29_E_vs_theta",
  "QED1N30_E_vs_theta", "QED1N31_E_vs_theta",

  "QED2N00_E_vs_theta", "QED2N01_E_vs_theta", "QED2N02_E_vs_theta", "QED2N03_E_vs_theta", "QED2N04_E_vs_theta", "QED2N05_E_vs_theta",
  "QED2N06_E_vs_theta", "QED2N07_E_vs_theta", "QED2N08_E_vs_theta", "QED2N09_E_vs_theta", "QED2N10_E_vs_theta", "QED2N11_E_vs_theta",
  "QED2N12_E_vs_theta", "QED2N13_E_vs_theta", "QED2N14_E_vs_theta", "QED2N15_E_vs_theta", "QED2N16_E_vs_theta", "QED2N17_E_vs_theta",
  "QED2N18_E_vs_theta", "QED2N19_E_vs_theta", "QED2N20_E_vs_theta", "QED2N21_E_vs_theta", "QED2N22_E_vs_theta", "QED2N23_E_vs_theta",
  "QED2N24_E_vs_theta", "QED2N25_E_vs_theta", "QED2N26_E_vs_theta", "QED2N27_E_vs_theta", "QED2N28_E_vs_theta", "QED2N29_E_vs_theta",
  "QED2N30_E_vs_theta", "QED2N31_E_vs_theta",

  "QED3N00_E_vs_theta", "QED3N01_E_vs_theta", "QED3N02_E_vs_theta", "QED3N03_E_vs_theta", "QED3N04_E_vs_theta", "QED3N05_E_vs_theta",
  "QED3N06_E_vs_theta", "QED3N07_E_vs_theta", "QED3N08_E_vs_theta", "QED3N09_E_vs_theta", "QED3N10_E_vs_theta", "QED3N11_E_vs_theta",
  "QED3N12_E_vs_theta", "QED3N13_E_vs_theta", "QED3N14_E_vs_theta", "QED3N15_E_vs_theta", "QED3N16_E_vs_theta", "QED3N17_E_vs_theta",
  "QED3N18_E_vs_theta", "QED3N19_E_vs_theta", "QED3N20_E_vs_theta", "QED3N21_E_vs_theta", "QED3N22_E_vs_theta", "QED3N23_E_vs_theta",
  "QED3N24_E_vs_theta", "QED3N25_E_vs_theta", "QED3N26_E_vs_theta", "QED3N27_E_vs_theta", "QED3N28_E_vs_theta", "QED3N29_E_vs_theta",
  "QED3N30_E_vs_theta", "QED3N31_E_vs_theta",

  "QED4N00_E_vs_theta", "QED4N01_E_vs_theta", "QED4N02_E_vs_theta", "QED4N03_E_vs_theta", "QED4N04_E_vs_theta", "QED4N05_E_vs_theta",
  "QED4N06_E_vs_theta", "QED4N07_E_vs_theta", "QED4N08_E_vs_theta", "QED4N09_E_vs_theta", "QED4N10_E_vs_theta", "QED4N11_E_vs_theta",
  "QED4N12_E_vs_theta", "QED4N13_E_vs_theta", "QED4N14_E_vs_theta", "QED4N15_E_vs_theta", "QED4N16_E_vs_theta", "QED4N17_E_vs_theta",
  "QED4N18_E_vs_theta", "QED4N19_E_vs_theta", "QED4N20_E_vs_theta", "QED4N21_E_vs_theta", "QED4N22_E_vs_theta", "QED4N23_E_vs_theta",
  "QED4N24_E_vs_theta", "QED4N25_E_vs_theta", "QED4N26_E_vs_theta", "QED4N27_E_vs_theta", "QED4N28_E_vs_theta", "QED4N29_E_vs_theta",
  "QED4N30_E_vs_theta", "QED4N31_E_vs_theta",

  "QED5N00_E_vs_theta", "QED5N01_E_vs_theta", "QED5N02_E_vs_theta", "QED5N03_E_vs_theta", "QED5N04_E_vs_theta", "QED5N05_E_vs_theta",
  "QED5N06_E_vs_theta", "QED5N07_E_vs_theta", "QED5N08_E_vs_theta", "QED5N09_E_vs_theta", "QED5N10_E_vs_theta", "QED5N11_E_vs_theta",
  "QED5N12_E_vs_theta", "QED5N13_E_vs_theta", "QED5N14_E_vs_theta", "QED5N15_E_vs_theta", "QED5N16_E_vs_theta", "QED5N17_E_vs_theta",
  "QED5N18_E_vs_theta", "QED5N19_E_vs_theta", "QED5N20_E_vs_theta", "QED5N21_E_vs_theta", "QED5N22_E_vs_theta", "QED5N23_E_vs_theta",
  "QED5N24_E_vs_theta", "QED5N25_E_vs_theta", "QED5N26_E_vs_theta", "QED5N27_E_vs_theta", "QED5N28_E_vs_theta", "QED5N29_E_vs_theta",
  "QED5N30_E_vs_theta", "QED5N31_E_vs_theta",

  "QED6N00_E_vs_theta", "QED6N01_E_vs_theta", "QED6N02_E_vs_theta", "QED6N03_E_vs_theta", "QED6N04_E_vs_theta", "QED6N05_E_vs_theta",
  "QED6N06_E_vs_theta", "QED6N07_E_vs_theta", "QED6N08_E_vs_theta", "QED6N09_E_vs_theta", "QED6N10_E_vs_theta", "QED6N11_E_vs_theta",
  "QED6N12_E_vs_theta", "QED6N13_E_vs_theta", "QED6N14_E_vs_theta", "QED6N15_E_vs_theta", "QED6N16_E_vs_theta", "QED6N17_E_vs_theta",
  "QED6N18_E_vs_theta", "QED6N19_E_vs_theta", "QED6N20_E_vs_theta", "QED6N21_E_vs_theta", "QED6N22_E_vs_theta", "QED6N23_E_vs_theta",
  "QED6N24_E_vs_theta", "QED6N25_E_vs_theta", "QED6N26_E_vs_theta", "QED6N27_E_vs_theta", "QED6N28_E_vs_theta", "QED6N29_E_vs_theta",
  "QED6N30_E_vs_theta", "QED6N31_E_vs_theta",
};

char dt_handles[N_DT][HANDLE_LENGTH]={
  "dt_ge_ge",      "dt_ge_bgo",     "dt_ge_sep",             "dt_ge_zds",     // 0-3
  "dt_ge_pac",     "dt_ge_labr",    "dt_ge_rcmp",            "dt_pac_zds",    // 4-7
  "dt_pac_labr",   "dt_rcmp_rcmp",  "dt_ge_art",             "dt_labr_art",   // 8-11
  "dt_paces_art",  "dt_art_art",    "dt_art_tac",            "dt_zds_tac",    // 12-15
  "dt_labr_tac",   "dt_labr_zds",   "dt_dsw_dsw",            "dt_dsw_ge",     // 16-19
  "dt_dsw_art",    "dt_dsw_zds",    "dt_zds_GRIF_CAEN_10ns", "dt_zds_GRIF_CAEN_2ns", // 20-23
  "dt_dsw_dsw_2ns","dt_dsw_zds_2ns","dt_labr_labr",          "dt_ge_qed",   "dt_qed_qed",   // 24-28
  "dt_comp_comp", "dt_comp_ge" };

  char dcfd_handles[N_DT][HANDLE_LENGTH]={
    "dcfd_ge_ge",      "dcfd_ge_bgo",     "dcfd_ge_sep",             "dcfd_ge_zds",     // 0-3
    "dcfd_ge_pac",     "dcfd_ge_labr",    "dcfd_ge_rcmp",            "dcfd_pac_zds",    // 4-7
    "dcfd_pac_labr",   "dcfd_rcmp_rcmp",  "dcfd_ge_art",             "dcfd_labr_art",   // 8-11
    "dcfd_paces_art",  "dcfd_art_art",    "dcfd_art_tac",            "dcfd_zds_tac",    // 12-15
    "dcfd_labr_tac",   "dcfd_labr_zds",   "dcfd_dsw_dsw",            "dcfd_dsw_ge",     // 16-19
    "dcfd_dsw_art",    "dcfd_dsw_zds",    "dcfd_zds_GRIF_CAEN_10ns", "dcfd_zds_GRIF_CAEN_2ns", // 20-23
    "dcfd_dsw_dsw_2ns","dcfd_dsw_zds_2ns","dcfd_labr_labr",          "dcfd_ge_qed",   "dcfd_qed_qed",   // 24-28
    "dt_comp_comp", "dt_comp_ge" };

    char gg_energy_handles[N_HPGE][HANDLE_LENGTH]={
      "GRG01BN00A_GGEnergy","GRG01GN00A_GGEnergy","GRG01RN00A_GGEnergy","GRG01WN00A_GGEnergy", "GRG02BN00A_GGEnergy","GRG02GN00A_GGEnergy","GRG02RN00A_GGEnergy","GRG02WN00A_GGEnergy",
      "GRG03BN00A_GGEnergy","GRG03GN00A_GGEnergy","GRG03RN00A_GGEnergy","GRG03WN00A_GGEnergy", "GRG04BN00A_GGEnergy","GRG04GN00A_GGEnergy","GRG04RN00A_GGEnergy","GRG04WN00A_GGEnergy",
      "GRG05BN00A_GGEnergy","GRG05GN00A_GGEnergy","GRG05RN00A_GGEnergy","GRG05WN00A_GGEnergy", "GRG06BN00A_GGEnergy","GRG06GN00A_GGEnergy","GRG06RN00A_GGEnergy","GRG06WN00A_GGEnergy",
      "GRG07BN00A_GGEnergy","GRG07GN00A_GGEnergy","GRG07RN00A_GGEnergy","GRG07WN00A_GGEnergy", "GRG08BN00A_GGEnergy","GRG08GN00A_GGEnergy","GRG08RN00A_GGEnergy","GRG08WN00A_GGEnergy",
      "GRG09BN00A_GGEnergy","GRG09GN00A_GGEnergy","GRG09RN00A_GGEnergy","GRG09WN00A_GGEnergy", "GRG10BN00A_GGEnergy","GRG10GN00A_GGEnergy","GRG10RN00A_GGEnergy","GRG10WN00A_GGEnergy",
      "GRG11BN00A_GGEnergy","GRG11GN00A_GGEnergy","GRG11RN00A_GGEnergy","GRG11WN00A_GGEnergy", "GRG12BN00A_GGEnergy","GRG12GN00A_GGEnergy","GRG12RN00A_GGEnergy","GRG12WN00A_GGEnergy",
      "GRG13BN00A_GGEnergy","GRG13GN00A_GGEnergy","GRG13RN00A_GGEnergy","GRG13WN00A_GGEnergy", "GRG14BN00A_GGEnergy","GRG14GN00A_GGEnergy","GRG14RN00A_GGEnergy","GRG14WN00A_GGEnergy",
      "GRG15BN00A_GGEnergy","GRG15GN00A_GGEnergy","GRG15RN00A_GGEnergy","GRG15WN00A_GGEnergy", "GRG16BN00A_GGEnergy","GRG16GN00A_GGEnergy","GRG16RN00A_GGEnergy","GRG16WN00A_GGEnergy"
    };

    char ge_bgo_handles[N_BGO][HANDLE_LENGTH]={
      "Ge01BGO01","Ge01BGO02","Ge01BGO03","Ge01BGO04","Ge01BGO05", "Ge02BGO01","Ge02BGO02","Ge02BGO03","Ge02BGO04","Ge02BGO05",
      "Ge03BGO01","Ge03BGO02","Ge03BGO03","Ge03BGO04","Ge03BGO05", "Ge04BGO01","Ge04BGO02","Ge04BGO03","Ge04BGO04","Ge04BGO05",
      "Ge05BGO01","Ge05BGO02","Ge05BGO03","Ge05BGO04","Ge05BGO05", "Ge06BGO01","Ge06BGO02","Ge06BGO03","Ge06BGO04","Ge06BGO05",
      "Ge07BGO01","Ge07BGO02","Ge07BGO03","Ge07BGO04","Ge07BGO05", "Ge08BGO01","Ge08BGO02","Ge08BGO03","Ge08BGO04","Ge08BGO05",
      "Ge09BGO01","Ge09BGO02","Ge09BGO03","Ge09BGO04","Ge09BGO05",
      "Ge10BGO01","Ge10BGO02","Ge10BGO03","Ge10BGO04","Ge10BGO05", "Ge11BGO01","Ge11BGO02","Ge11BGO03","Ge11BGO04","Ge11BGO05",
      "Ge12BGO01","Ge12BGO02","Ge12BGO03","Ge12BGO04","Ge12BGO05", "Ge13BGO01","Ge13BGO02","Ge13BGO03","Ge13BGO04","Ge13BGO05",
      "Ge14BGO01","Ge14BGO02","Ge14BGO03","Ge14BGO04","Ge14BGO05", "Ge15BGO01","Ge15BGO02","Ge15BGO03","Ge15BGO04","Ge15BGO05",
      "Ge16BGO01","Ge16BGO02","Ge16BGO03","Ge16BGO04","Ge16BGO05", "Ge17BGO01","Ge17BGO02","Ge17BGO03","Ge17BGO04","Ge17BGO05",
      "Ge18BGO01","Ge18BGO02","Ge18BGO03","Ge18BGO04","Ge18BGO05", "Ge19BGO01","Ge19BGO02","Ge19BGO03","Ge19BGO04","Ge19BGO05",
      "Ge20BGO01","Ge20BGO02","Ge20BGO03","Ge20BGO04","Ge20BGO05", "Ge21BGO01","Ge21BGO02","Ge21BGO03","Ge21BGO04","Ge21BGO05",
      "Ge22BGO01","Ge22BGO02","Ge22BGO03","Ge22BGO04","Ge22BGO05", "Ge23BGO01","Ge23BGO02","Ge23BGO03","Ge23BGO04","Ge23BGO05",
      "Ge24BGO01","Ge24BGO02","Ge24BGO03","Ge24BGO04","Ge24BGO05", "Ge25BGO01","Ge25BGO02","Ge25BGO03","Ge25BGO04","Ge25BGO05",
      "Ge26BGO01","Ge26BGO02","Ge26BGO03","Ge26BGO04","Ge26BGO05", "Ge27BGO01","Ge27BGO02","Ge27BGO03","Ge27BGO04","Ge27BGO05",
      "Ge28BGO01","Ge28BGO02","Ge28BGO03","Ge28BGO04","Ge28BGO05", "Ge29BGO01","Ge29BGO02","Ge29BGO03","Ge29BGO04","Ge29BGO05",
      "Ge30BGO01","Ge30BGO02","Ge30BGO03","Ge30BGO04","Ge30BGO05", "Ge31BGO01","Ge31BGO02","Ge31BGO03","Ge31BGO04","Ge31BGO05",
      "Ge32BGO01","Ge32BGO02","Ge32BGO03","Ge32BGO04","Ge32BGO05", "Ge33BGO01","Ge33BGO02","Ge33BGO03","Ge33BGO04","Ge33BGO05",
      "Ge34BGO01","Ge34BGO02","Ge34BGO03","Ge34BGO04","Ge34BGO05", "Ge35BGO01","Ge35BGO02","Ge35BGO03","Ge35BGO04","Ge35BGO05",
      "Ge36BGO01","Ge36BGO02","Ge36BGO03","Ge36BGO04","Ge36BGO05", "Ge37BGO01","Ge37BGO02","Ge37BGO03","Ge37BGO04","Ge37BGO05",
      "Ge38BGO01","Ge38BGO02","Ge38BGO03","Ge38BGO04","Ge38BGO05", "Ge39BGO01","Ge39BGO02","Ge39BGO03","Ge39BGO04","Ge39BGO05",
      "Ge40BGO01","Ge40BGO02","Ge40BGO03","Ge40BGO04","Ge40BGO05", "Ge41BGO01","Ge41BGO02","Ge41BGO03","Ge41BGO04","Ge41BGO05",
      "Ge42BGO01","Ge42BGO02","Ge42BGO03","Ge42BGO04","Ge42BGO05", "Ge43BGO01","Ge43BGO02","Ge43BGO03","Ge43BGO04","Ge43BGO05",
      "Ge44BGO01","Ge44BGO02","Ge44BGO03","Ge44BGO04","Ge44BGO05", "Ge45BGO01","Ge45BGO02","Ge45BGO03","Ge45BGO04","Ge45BGO05",
      "Ge46BGO01","Ge46BGO02","Ge46BGO03","Ge46BGO04","Ge46BGO05", "Ge47BGO01","Ge47BGO02","Ge47BGO03","Ge47BGO04","Ge47BGO05",
      "Ge48BGO01","Ge48BGO02","Ge48BGO03","Ge48BGO04","Ge48BGO05", "Ge49BGO01","Ge49BGO02","Ge49BGO03","Ge49BGO04","Ge49BGO05",
      "Ge50BGO01","Ge50BGO02","Ge50BGO03","Ge50BGO04","Ge50BGO05", "Ge51BGO01","Ge51BGO02","Ge51BGO03","Ge51BGO04","Ge51BGO05",
      "Ge52BGO01","Ge52BGO02","Ge52BGO03","Ge52BGO04","Ge52BGO05", "Ge53BGO01","Ge53BGO02","Ge53BGO03","Ge53BGO04","Ge53BGO05",
      "Ge54BGO01","Ge54BGO02","Ge54BGO03","Ge54BGO04","Ge54BGO05", "Ge55BGO01","Ge55BGO02","Ge55BGO03","Ge55BGO04","Ge55BGO05",
      "Ge56BGO01","Ge56BGO02","Ge56BGO03","Ge56BGO04","Ge56BGO05", "Ge57BGO01","Ge57BGO02","Ge57BGO03","Ge57BGO04","Ge57BGO05",
      "Ge58BGO01","Ge58BGO02","Ge58BGO03","Ge58BGO04","Ge58BGO05", "Ge59BGO01","Ge59BGO02","Ge59BGO03","Ge59BGO04","Ge59BGO05",
      "Ge60BGO01","Ge60BGO02","Ge60BGO03","Ge60BGO04","Ge60BGO05", "Ge61BGO01","Ge61BGO02","Ge61BGO03","Ge61BGO04","Ge61BGO05",
      "Ge62BGO01","Ge62BGO02","Ge62BGO03","Ge62BGO04","Ge62BGO05", "Ge63BGO01","Ge63BGO02","Ge63BGO03","Ge63BGO04","Ge63BGO05",
      "Ge64BGO01","Ge64BGO02","Ge64BGO03","Ge64BGO04","Ge64BGO05"
    };

    //#######################################################################
    //########                Histogram pointers                   ##########
    //#######################################################################

    TH1I   *ge_cycle_activity, *zds_cycle_activity; // Activity over cycle time, sum of all cycles
    TH2I   *ge_e_vs_cycle_time;                     // Energy vs time within the cycle
    TH1I   *ge_cycle_code[N_PPG_PATTERNS];          // Energy spectrum for each PPG pattern
    TH2I   *gg_cycle_code[N_PPG_PATTERNS];          // Ge-Ge 2D histogram for each PPG pattern
    TH1I   *gea_cycle_num[MAX_CYCLES];               // Activity over cycle time for each indivdual cycle, GRGA
    TH1I   *gea_cycle_num_sh[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGA, single_hit only
    TH1I   *gea_cycle_num_pu[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGA, pileup only
    TH1I   *gea_cycle_num_dt[MAX_CYCLES];            // Deadtime over cycle time for each indivdual cycle, GRGA
    TH1I   *gea_cycle_num_g[MAX_CYCLES];             // Activity over cycle time for each indivdual cycle, GRGA, gamma-gated
    TH1I   *gea_cycle_num_sh_g[MAX_CYCLES];          // Activity over cycle time for each indivdual cycle, GRGA, gamma-gated, single_hit only
    TH1I   *geb_cycle_num[MAX_CYCLES];               // Activity over cycle time for each indivdual cycle, GRGB
    TH1I   *geb_cycle_num_sh[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGB, single_hit only
    TH1I   *geb_cycle_num_pu[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGB, pileup only
    TH1I   *geb_cycle_num_dt[MAX_CYCLES];            // Deadtime over cycle time for each indivdual cycle, GRGB
    TH1I   *geb_cycle_num_g[MAX_CYCLES];             // Activity over cycle time for each indivdual cycle, GRGB, gamma-gated
    TH1I   *geb_cycle_num_sh_g[MAX_CYCLES];          // Activity over cycle time for each indivdual cycle, GRGB, gamma-gated, single_hit only
    TH2I   *cycle_num_vs_ge;                        // 2D histogram of cycle # vs the cycle time.
    TH2I   *cycle_num_vs_ge_sh_g;                   // 2D histogram of cycle # vs the cycle time for 1809 gated/NP.
    TH2I   *cycle_num_vs_ge_dt;                     // 2D histogram of cycle # vs the cycle time for deadtime .
    TH2I   *cycle_num_vs_sh;                        // 2D histogram of cycle # vs the cycle time for NP.
    TH2I   *cycle_num_vs_pu;                        // 2D histogram of cycle # vs the cycle time for PU.
    TH2I   *cycle_num_vs_ge_b;                      // 2D histogram of cycle # vs the cycle time. GRGB
    TH2I   *cycle_num_vs_ge_b_sh_g;                 // 2D histogram of cycle # vs the cycle time for 1809 gated/NP. GRGB
    TH2I   *cycle_num_vs_ge_b_dt;                   // 2D histogram of cycle # vs the cycle time for deadtime . GRGB
    TH2I   *cycle_num_vs_sh_b;                      // 2D histogram of cycle # vs the cycle time for NP. GRGB
    TH2I   *cycle_num_vs_pu_b;                      // 2D histogram of cycle # vs the cycle time for PU. GRGB
    TH2I   *cycle_num_vs_geEnergy[N_HPGE];          // 2D histogram of cycle # vs the Ge energy spectrum for that cycle (Use for monitoring gain drifts).
    TH2I   *cycle_num_vs_qedEnergy[N_QED_POS];      // 2D histogram of cycle # vs the Ge energy spectrum for that cycle (Use for monitoring gain drifts).
    // Pulse height and energy
    TH1I   *ts_hist; // timestamp
    TH1I   *gc_hist; // GRIF-CAEN hitpattern for checking coincidences
    TH1I   *ph_hist[MAX_DAQSIZE];
    TH1I   *e_hist[MAX_DAQSIZE];
    //TH1I *wave_hist[MAX_DAQSIZE];

    TH1I  *hit_hist[N_HITPAT], *mult_hist[MAX_SUBSYS];

    // DESCANT Wall
    TH1I  *desw_tof[N_DES_WALL];                // Time-Of-Flight
    TH1I  *desw_tof_corr[N_DES_WALL];           // corrected Time-Of-Flight
    TH1I  *desw_tof_psd[N_DES_WALL];            // corrected Time-Of-Flight, PSD gated
    TH1I  *desw_psd[N_DES_WALL];                // Pulse Shape Discrimination
    TH1I  *ge_pu_type; // The value of pile-up type
    TH1I  *ge_nhits_type; // The value of nhits type
    TH1I  *ge_pu_class; // The value of pile-up class
    TH1I  *ge_sum_class[N_PU_CLASSES]; // Ge energy for each pile-up value
    TH2I  *ge_e_vs_k_class[N_PU_CLASSES]; // Ge energy vs k for each pile-up value
    TH2I  *ge_xtal_1hit, *ge_xtal_2hit, *ge_xtal_3hit; // Ge energy vs crystal number for 1, 2, 3 hit PU.
    TH1I  *ge_1hit[N_HPGE]; // Ge single hit events
    TH1I  *ge_2hit[N_HPGE]; // Ge 2-hit pileup events
    TH1I  *ge_3hit[N_HPGE]; // Ge 3-hit pileup events
    TH1I  *ge_pu_dt12; // Time difference between first and second Hit
    TH1I  *ge_pu_dt13; // Time difference between first and third Hit
    TH2I  *ge_e_vs_k_2hit_first[N_HPGE]; // Ge Hit1 energy vs k1 for 2-hit pileup events
    TH2I  *ge_e_vs_k_2hit_second[N_HPGE]; // Ge Hit2 energy vs k2 for 2-hit pileup events
    TH2I  *ge_PU2_e2_v_k_gatedxrays[N_HPGE]; // Ge e2 vs k2 for fixed e1 energy gates
    TH2I  *ge_PU2_e2_v_k_gated1408[N_HPGE]; // Ge e2 vs k2 for fixed e1 energy gates

    // for most subsystem-pair-combinations, there is a
    // a 1d time-difference and a 2d ecal-vs-ecal matrix
    TH2I *subsys_e_vs_e[MAX_SUBSYS][MAX_SUBSYS];
    TH1I *subsys_dt[MAX_SUBSYS][MAX_SUBSYS], *subsys_dcfd[MAX_SUBSYS][MAX_SUBSYS];
    TH1I *tac_lbl_ts_diff[N_TACS];

    // HPGe (ge_sum is sum of crystal energies, ge_sum_b is beta-gated)
    TH1I  *ge_ab_e[N_CLOVER], *ge_ab_sup_e[N_CLOVER], *ge_sum_ab, *ge_sum_ab_sup, *ge_sum_ab_sup_rej;
    TH1I  *ge_sum, *ge_sum_us, *ge_sum_ds, *ge_sum_hem[2], *ge_sum_ab_us, *ge_sum_ab_ds;
    TH1I  *ge_sum_b, *ge_sum_b_ab, *ge_sum_b_sep, *ge_sum_b_sep_brems, *ge_sum_b_ab_sep_brems, *ge_sum_b_zds;
    TH1I  *ge_sum_b_art, *ge_sum_b_art_brems, *ge_sum_b_artT, *ge_sum_b_artR, *ge_sum_b_artS;

    // ARIES, PACES and LaBr3
    TH1I  *aries_sum;  // aries_sum is sum of tile energies
    TH1I  *paces_sum, *paces_sum_b;  // paces_sum is sum of crystal energies, *paces_sum_b has a beta-coincidence
    TH1I  *labr_sum;  // labr_sum is sum of crystal energies

    // RCMP
    TH1I  *rcmp_sum, *rcmp_fb_sum;  // rcmp_sum is sum of strip energies, fb is with front-back coincidence
    TH2I  *rcmp_strips[N_RCMP_POS];
    TH2I  *rcmp_hit[N_RCMP_POS];
    TH2I  *rcmp_fb[N_RCMP_POS];
    TH2I  *rcmp_x_ge_hit, *rcmp_y_ge_hit; // rcmp strips vs Ge hitpatterns

    // QED
    TH1I  *qed_sum, *qed_fb_sum;  // qed_sum is sum of strip energies, fb is with front-back coincidence
    TH2I  *qed_strips[N_QED_POS];
    TH2I  *qed_hit[N_QED_POS];
    TH2I  *qed_fb[N_QED_POS];
    TH2I  *qed_p_ge_hit, *qed_n_ge_hit; // qed strips vs Ge hitpatterns
    TH2I  *qedp_ge_theta[N_QED_POS*N_QED_STRIPS], *qedn_ge_theta[N_QED_POS*N_QED_STRIPS]; // qed strip energy vs theta of a qed-Ge hit
    TH2I  *qedE_ge_theta_sum, *qed_geE_theta_sum, *qed_E_totE_sum_t, *qed_geE_totE_sum_t, *qedE_ge_theta_sum_t, *qed_geE_theta_sum_t, *qedE_ge_thetaI_sum_t, *qed_geE_thetaDiff_sum_t, *qed_geE_thetaI_sum_t;
    TH2I  *qedE_ge_theta_sum_c, *qed_geE_theta_sum_c, *qedE_ge_theta_sum_c_g, *qed_geE_theta_sum_c_g, *qedE_ge_theta_sum_c_s, *qed_geE_theta_sum_c_s, *ge_qed_c;
    TH2I  *qed_geE_theta_clov[N_CLOVER], *qed_geE_theta_clov_t[N_CLOVER], *qed_E_theta_dssd[N_QED_POS], *qed_geE_theta_dssd[N_QED_POS];
    TH2I  *qedE_ge_dt, *qed_geE_dt, *qedE_ge_dt_c, *qed_geE_dt_c, *qed_theta_dt_c, *qed_theta_dt, *qed_theta_dt_cfd, *qed_dcs_omega_dt, *qed_dcs_omega_dtx, *qedx_dcs_omega_dt[N_QED_POS], *qed_theta1_vs_theta2, *qed_theta1_azi, *qed_theta2_azi, *qed2_theta1_vs_theta2, *qed2_theta1_azi, *qed2_theta2_azi;
    TH1I  *qed_dcs_omega, *qed_dcs_omega_t, *qed_dcs_azi, *qed_dcs_azi_t, *qed_dcs_azi_tg, *qed_dcs_azi_TRWF, *qed_dcs_azi_TRWF_t, *qed_dcs_azi_TRWF_tg, *qed_delta_theta1_theta2, *qed_sum_theta1_theta2;
    TH1I  *qed_wf_dcs_azi, *qed_wf_omega;
    TH2I  *qed_angle_theta_g, *qed_angle_theta_s, *qed_angle_phi_s, *qed_angle_phi_a, *qed_angle_phi_b;
    TH1I  *qed_theta, *qed_phi_s, *qed_phi_a, *qed_phi_b;
    TH1I  *qed_dcs_azi_bins1, *qed_dcs_azi_bins2, *qed_dcs_azi_bins3, *qed_dcs_azi_bins4, *qed_dcs_azi_bins5, *qed_dcs_azi_bins6, *qed_dcs_azi_bins7, *qed_dcs_azi_bins8, *qed_dcs_azi_bins8a, *qed_dcs_azi_bins9, *qed_dcs_azi_bins10;
    TH1I  *qed_dcs_azi_TRWF_bins1, *qed_dcs_azi_TRWF_bins2, *qed_dcs_azi_TRWF_bins3, *qed_dcs_azi_TRWF_bins4, *qed_dcs_azi_TRWF_bins5, *qed_dcs_azi_TRWF_bins6, *qed_dcs_azi_TRWF_bins7, *qed_dcs_azi_TRWF_bins8, *qed_dcs_azi_TRWF_bins8a, *qed_dcs_azi_TRWF_bins9, *qed_dcs_azi_TRWF_bins10;
    TH1I  *qed_phi_bins1, *qed_phi_bins2, *qed_phi_bins3, *qed_phi_bins4, *qed_phi_bins5, *qed_phi_bins6, *qed_phi_bins7, *qed_phi_bins8, *qed_phi_bins9, *qed_phi_bins10;
    TH2I  *dcsaE_ge_theta, *dcsa_geE_theta, *dcsa_theta_azi, *qed_dcs_omega_dt_TRWF, *dcsa_theta_azi_ge;
    TH1I  *dcsa_cs_omega, *dcsa_cs_omega_ge, *dcsa_theta;
    TH2I  *dcsbE_ge_theta, *dcsb_geE_theta, *dcsb_theta_azi, *dcsb_theta_azi_ge;
    TH1I  *dcsb_cs_omega, *dcsb_cs_omega_ge, *dcsb_theta;
    TH2I  *qed_qed_23, *qed_qed_23_theta2, *qed_qed_23_theta3, *qed_qed_23_totv2, *qed_qed_23_totv3, *qed_qed_12, *qed_qed_12_theta1, *qed_qed_12_theta2, *qed_qed_12_totv1, *qed_qed_12_totv2, *qed_qed_14, *qed_qed_14_theta1, *qed_qed_14_theta4, *qed_qed_14_totv1, *qed_qed_14_totv4;
    TH1I  *qed_qed_23dt, *qed_qed_12dt, *qed_qed_14dt;
    TH2I  *qed_ge_weight;

    // DESCANT WALL
    TH1I  *desw_sum_e, *desw_sum_tof, *desw_sum_psd;  // Sums of energies and corTOF and PSD
    TH1I  *desw_sum_e_b, *desw_sum_tof_b;       // Beta-tagged Sums of energies and corTOF
    TH1I  *desw_sum_e_nn, *desw_sum_tof_nn;     // fold>2 Sums of energies and corTOF
    TH1I  *desw_sum_e_nn_a, *desw_sum_tof_nn_a; // fold>2, angle>60 Sums of energies and corTOF
    TH2I  *desw_psd_e, *desw_psd_tof;           // PSD vs energies or corrected-TOF
    TH2I  *desw_psd_q,*desw_psd_cc,*desw_q_cc,*desw_q_tof,*desw_cc_tof,*desw_psd_zdse; //

    // TAC spectra
    TH1I *tac_labr_hist[(int)((N_LABR)*(N_LABR-1)/2)+2]; // this index numbers are the LaBr-LaBr position numbers
    TH1I *tac_labr_hist_uncal[(int)((N_LABR)*(N_LABR-1)/2)+2]; // this index numbers are the LaBr-LaBr position numbers
    // One additional histogram (2_1) needed for Compton Walk corrections
    TH2I *tac_labr_CompWalk[N_LABR];         // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
    TH2I *tac_labr_CompWalk0;                // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
    TH1I *tac_gated_lbl[N_LABR];             // TAC-gated LBL energy spectrum to check anode threshold in analogue CFD
    TH1I *final_tac[N_TACS], *final_tac_sum; // Final TAC spectra after all calibration and alignment
    TH1I *tac_aries_lbl[N_LABR];             // this index number is the LaBr position number
    TH1I *tac_aries_art[N_ARIES];            // this index number is the Aries position number
    TH1I *tac_aries_lbl_sum;                 // ARIES TAC sum spectrum of all LBLs
    TH1I *tac_aries_art_sum;                 // ARIES TAC sum spectrum of all ARTs
    TH1I *aries_tac;                         // aries_tac gated on 1275keV peak
    TH1I *aries_tac_Egate;                   // aries_tac gated on 1275keV peak
    TH1I *aries_tac_artEn;                   // aries energy in coincidence with TAC
    TH2I *lblE_tac, *zdsE_tac, *ariesE_tac;  // lbl or zds or aries energy vs TAC
    TH2I *lbl_lbl_tac;                       // A special 3d histogram disguised as a 2d histogram

    // 2D Energy vs detector number
    TH2I *ge_xtal, *geb_xtal, *bgo_xtal, *bgof_xtal, *bgos_xtal, *bgob_xtal, *bgoa_xtal, *labr_xtal;
    TH2I *labr_tac_xtal, *paces_xtal, *sceptar_xtal, *aries_xtal, *art_tac_xtal, *desw_e_xtal, *desw_tof_xtal;

    TH1I  *dt_hist[N_DT], *dcfd_hist[N_DT], *dt_tacs_hist[N_TACS];

    // 2D hitpatterns
    TH2I *gg_hit, *bgobgo_hit, *aa_hit, *gea_hit, *lba_hit, *dsw_hit;

    // 2D Energy vs Energy Coincidence matrices
    TH2I *gea_self_dt,*geb_self_dt;
    TH2I *gg, *gg_ab, *gg_opp, *gg_ab_opp, *ge_bgo, *ge_paces, *ge_labr, *ge_rcmp, *labr_labr, *labr_zds, *labr_rcmp;
    TH2I *ge_art, *ge_zds, *paces_art, *labr_art, *art_art, *dsw_dsw, *ge_dsw, *art_dsw, *ge_qed, *qed_qed, *comp_comp, *ge_comp, *geadd_comp, *ge_dcs, *geadd_dcs, *comp_dcs;
    TH1I *gg_energy[N_HPGE];

    // Angular Correlation histograms
    TH2I  *gg_angcor_110[N_GE_ANG_CORR], *gg_angcor_145[N_GE_ANG_CORR], *ge_art_angcor[N_GRG_ART_ANG_CORR], *dsw_angcor[N_DSW_DSW_ANG_CORR];

    // Compton Polarimetry histograms
    TH2I  *comp_pol_angles_110, *comp_pol_angles_145, *gg_comp_pol_110[N_GE_COMP_POL], *gg_comp_pol_145[N_GE_COMP_POL];

    // Isomer Spectroscopy
    TH2I  *gg_dt, *gb_dt;
    TH1I  *ge_isomer_popu, *ge_isomer_depop;

    // Crosstalk Analysis
    TH2I  *ct_e_vs_dt_B[N_HPGE], *ct_e_vs_dt_G[N_HPGE], *ct_e_vs_dt_R[N_HPGE], *ct_e_vs_dt_W[N_HPGE];

    ////////////////////////////////////
    ////////////////////////////////////

    extern int tac_ts_offset[12]; // LBT (TAC) timestamp offset values.
    extern int tac_lbl_combo_offset[(int)((N_LABR)*(N_LABR-1)/2)+2]; // TAC coincidence combination offsets.

    // BGO HV alignment histograms
    TH1I *ge_bgo_gated[N_BGO];



    //#######################################################################
    //########                PRESORT Time Gates                   ##########
    //#######################################################################

    // The definition of the time difference gate in 10 nanosecond units.
    // The value is the maximum time difference in 10 nanosecond units.
    // The default values set here are replaced by the Global value at start of sorting.
    int bgo_window_min       =  0;
    int addback_window_min   =  0;
    int rcmp_fb_window_min   =  0;
    int qed_fb_window_min    =  0;
    int lbl_tac_window_min   =  0;
    int art_tac_window_min   =  0;
    int zds_tac_window_min   =  0;
    int desw_beta_window_min =  0;
    int bgo_window_max       = 20;
    int addback_window_max   = 20;
    int rcmp_fb_window_max   = 10;
    int qed_fb_window_max    = 10;
    int lbl_tac_window_max   = 25;
    int art_tac_window_max   = 25;
    int zds_tac_window_max   = 25;
    int desw_beta_window_max = 80;

    // The definition of the time difference gate in 10 nanosecond units.
    // First and second index are the subsystem index numbers
    // The value is the maximum time difference in 10 nanosecond units.
    // Default is 250 nanoseconds, replaced by the Global value at start of sorting
    int time_diff_gate_min[MAX_SUBSYS][MAX_SUBSYS];
    int time_diff_gate_max[MAX_SUBSYS][MAX_SUBSYS];


    int tac_labr_hist_index[N_LABR][N_LABR]; // index for filling tac_labr_hist from LBL id numbers
    // LBT (TAC) timestamp offset values.
    // These are subtracted in the apply_gains function
    //int tac_ts_offset[12] = {60,60,60,60,60,60,60,60,60,60,60,60}; // From Dec 2024
    //  int tac_ts_offset[12] = { 54, 64,406,137, 77,404,114,158, 0, 0, 0, 0}; // S2231_S2196_Nov2024
    // int tac_ts_offset[12] = {134, 48, 74, 59, 48,400,395,  0, 0, 0, 0, 0}; // S1723, Aug 2021
    int tac_ts_offset[12] = {60,60,60,60,60,60,60,60,60,60,60,60}; // Default 60 for Dec 2024 onwards. Set as Globals otherwise


    // TAC coincidence combination offsets.
    // These should not be here, they are a calibration and should be settable by the user.

    // These are the offsets for aligning the TAC calibrated energy based on the two LaBr3 energies in coincidence
    // 30 values
    // They are set from Globals
    int tac_lbl_combo_offset[(int)((N_LABR)*(N_LABR-1)/2)+2] = {
      0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
      0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0
    };

    //#######################################################################
    //########                FUNCTIONS                            ##########
    //#######################################################################

    extern void init_acos_table(void);
    // odb tables need to be transferred into config, which is saved with histos
    int init_default_histos(Config *cfg, Sort_status *arg)
    {
      Cal_coeff *cal;
      int i, j;

      // Initialize all pileup and crosstalk parameters to unset values
      for(i=0; i<odb_daqsize; i++){
        for(j=0; j<7; j++){
          pileupk1[i][j] = pileupk2[i][j] = pileupE1[i][j] = -1;
        }
        for(j=0; j<16; j++){
          crosstalk[i][0][j] = crosstalk[i][1][j] = crosstalk[i][2][j] = -1;
        }
      }

      cfg->odb_daqsize = odb_daqsize;
      for(i=0; i<odb_daqsize; i++){ // update config with odb info
        edit_calibration(cfg, chan_name[i], offsets[i], gains[i], quads[i], pileupk1[i], pileupk2[i], pileupE1[i],
          crosstalk[i][0], crosstalk[i][1], crosstalk[i][2], chan_address[i],  dtype_table[i], arg->cal_overwrite );
        }
        // ALSO need to transfer config info to the arrays that are used in sort
        for(i=0; i<odb_daqsize; i++){

          cal = cfg->calib[i];
          if( strcmp(chan_name[i], cal->name) != 0 ){ // conf not in odb order
            for(j=0; j<cfg->ncal; j++){ cal = cfg->calib[j];
              if( strcmp(chan_name[i], cal->name) == 0 ){ break; }
            }
            if( j == cfg->ncal ){ continue; } // not found in config
          }

          // overwrite = 0 => USE CONFIG NOT ODB for offset, gain, quads
          if( arg->cal_overwrite == 0 ){
            offsets[i]=cal->offset; gains[i]=cal->gain;  quads[i]=cal->quad;
          }

          // Pileup or crosstalk parameters do not exist in the MIDAS ODB so must always be copied from the config (file or config)
          if(strncmp(chan_name[i],"GRG",3)==0){
            for(j=0; j<7; j++){
              pileupk1[i][j] = (isnan(cal->pileupk1[j])) ? 0.0 : cal->pileupk1[j];
              pileupk2[i][j] = (isnan(cal->pileupk2[j])) ? 0.0 : cal->pileupk2[j];
              pileupE1[i][j] = (isnan(cal->pileupE1[j])) ? 0.0 : cal->pileupE1[j];
            }
            for(j=0; j<16; j++){
              crosstalk[i][0][j] = (isnan(cal->crosstalk0[j])) ? 0.0 : cal->crosstalk0[j];
              crosstalk[i][1][j] = (isnan(cal->crosstalk1[j])) ? 0.0 : cal->crosstalk1[j];
              crosstalk[i][2][j] = (isnan(cal->crosstalk2[j])) ? 0.0 : cal->crosstalk2[j];
            }
          }
        }

        init_acos_table();
        init_parameters_from_globals(cfg);
        init_chan_histos(cfg);
        init_histos(cfg, SUBSYS_HPGE_A); // always create Ge histos

        // Reset the deadtime counters and previous_trig_acc at BOR
        memset(subsys_deadtime_count,0,MAX_SUBSYS*sizeof(int));
        memset(previous_trig_acc,0,MAX_DAQSIZE*sizeof(int));
        ppg_bin_end = ppg_cycles_binning_factor;

        return(0);
      }

      int init_chan_histos(Config *cfg)
      {                      // 1d histograms for Q,E,T,Wf for each channel in odb
        char title[STRING_LEN], handle[STRING_LEN];
        int i, j, k, pos;

        open_folder(cfg, "Hits_and_Sums");
        open_folder(cfg, "Hits");
        for(i=0; i<N_HITPAT; i++){ // Create Hitpattern spectra
          sprintf(title,  "Hitpattern_%s",    hit_names[i] );
          hit_hist[i] = H1_BOOK(cfg, hit_handles[i], title, MAX_DAQSIZE, 0, MAX_DAQSIZE);
        }
        ts_hist = H1_BOOK(cfg, "ts", "Timestamp", 163840, 0, 163840);
        gc_hist = H1_BOOK(cfg, "gc", "ZDS GRIF-CAEN", 16, 0, 16);
        close_folder(cfg);
        open_folder(cfg, "Multiplicities");
        for(i=0; i<MAX_SUBSYS; i++){ mult_hist[i] = NULL;
          if( strncmp(subsys_handle[i],"XXX",3) == 0 ||
          strlen(subsys_handle[i]) == 0 ){ continue; }
          sprintf(title,  "%s_Multiplicity", subsys_handle[i] );
          sprintf(handle, "%s_Mult",         subsys_handle[i] );
          mult_hist[i] = H1_BOOK(cfg, handle, title, MULT_SPEC_LENGTH, 0, MULT_SPEC_LENGTH);
        }
        close_folder(cfg);
        close_folder(cfg);
        for(i=0; i<MAX_DAQSIZE; i++){
          if( i >= odb_daqsize ){ break; }
          for(j=0; j<MAX_SUBSYS-1; j++){ // get subsystem name from chan name
            if( strncmp(subsys_handle[j],"XXX",3) == 0 ||
            strlen(subsys_handle[j]) == 0 ){ continue; }
            if( memcmp(subsys_handle[j], chan_name[i], 3) == 0 ){ break; }
          }                              // stop at final entry "unknown"
          if(j==(MAX_SUBSYS-1)){ if(strncmp(chan_name[i],"GRS",3) == 0){ j=8; } } // GRS != BGO
          open_folder(cfg, subsys_name[j]);
          open_folder(cfg, "Energy");
          sprintf(title,  "%s_Energy",         chan_name[i] );
          sprintf(handle, "%s_E",              chan_name[i] );
          if( strcmp(subsys_handle[j],"LBT") == 0){
            e_hist[i] = H1_BOOK(cfg, handle, title, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
          }else{
            e_hist[i] = H1_BOOK(cfg, handle, title, E_SPECLEN, 0, E_SPECLEN);
          }
          close_folder(cfg);
          //      open_folder(cfg, "Waveform");
          //      sprintf(title,  "%s_Waveform",       chan_name[i] );
          //      sprintf(handle, "%s_w",              chan_name[i] );
          //      wave_hist[i] = H1_BOOK(cfg, handle, title, WV_SPEC_LENGTH, 0, WV_SPEC_LENGTH);
          //      close_folder(cfg);
          open_folder(cfg, "PulseHeight");
          sprintf(title,  "%s_Pulse_Height",   chan_name[i] );
          sprintf(handle, "%s_Q",              chan_name[i] );
          if( strcmp(subsys_handle[j],"LBT") == 0){
            ph_hist[i] = H1_BOOK(cfg, handle, title, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
          }else{
            ph_hist[i] = H1_BOOK(cfg, handle, title, E_SPECLEN, 0, E_SPECLEN);
          }
          close_folder(cfg);
          if( strcmp(subsys_handle[j],"DSW") == 0){
            pos  = crystal_table[i];
            if(pos>0 && pos<=N_DES_WALL){
              open_folder(cfg, "PSD");
              sprintf(title,  "%s_PSD",         chan_name[i] );
              sprintf(handle, "%s_PSD",         chan_name[i] );
              desw_psd[pos] = H1_BOOK(cfg, handle, title, E_PSD_SPEC_LENGTH, 0, E_PSD_SPEC_LENGTH);
              close_folder(cfg);
              open_folder(cfg, "Time_Of_Flight");
              sprintf(title,  "%s_TOF_PSD-gated", chan_name[i] );
              sprintf(handle, "%s_CTOF_PSDn",         chan_name[i] );
              desw_tof_psd[pos] = H1_BOOK(cfg, handle, title, E_TOF_SPEC_LENGTH, 0, E_TOF_SPEC_LENGTH);
              sprintf(title,  "%s_Corrected_TOF", chan_name[i] );
              sprintf(handle, "%s_CTOF",         chan_name[i] );
              desw_tof_corr[pos] = H1_BOOK(cfg, handle, title, E_TOF_SPEC_LENGTH, 0, E_TOF_SPEC_LENGTH);
              sprintf(title,  "%s_TOF",         chan_name[i] );
              sprintf(handle, "%s_TOF",         chan_name[i] );
              desw_tof[pos] = H1_BOOK(cfg, handle, title, E_TOF_SPEC_LENGTH, 0, E_TOF_SPEC_LENGTH);
              close_folder(cfg);
            }
          }
          close_folder(cfg);
          if( strncmp(subsys_handle[j],"XXX",3) == 0 ){     // suppress XXX channels
            ph_hist[i]->suppress = e_hist[i]->suppress =
            /*  cfd_hist[i]->suppress = wave_hist[i]->suppress = */ 1;
          }
        }
        return(0);
      }


      // Only create histograms when that subsystem is seen in the data
      // Do not create unused histos -> need to store subsystem types along with definitions ...
      // NOTE: only store one subsystem (even for 2d) - choose the most likely to save space/time
      typedef struct histo_def_struct {
        void **ptr; char title[80]; void *handle;
        int subsys; int xchan; int ychan; int count;
      } Histogram_definition;
      // histo folders have ptr=NULL and xchan=0
      // count != 0 is an array of histos
      //    either - array title will contain %d (padding such as %02d is allowed). The index will run from 0 to (count-1)
      //    or     - title is "" (a blank string) and handle is pointer to array of handles.
      // NOTE that if %d or an array of handles are used then the **ptr should be a pointer and not an address (ie. do not include the leading &)
      // --------------------
      // Set ychan to "SYMMETERIZE" to create a fully symmeterized 2d histogram
      // --------------------
      // default handles derived from title with the following substitutions
      //    "Energy"->"E", "CrystalNum"->"Xtal", "Des_Wall"->"DESW"
      //    "Upstream"->"US", "Downstream"->"DS", ""->"", ""->"",
      #define HISTO_DEF_SIZE 512
      Histogram_definition histodef_array[HISTO_DEF_SIZE] = {
        // Singles
        {NULL,                   "Hits_and_Sums/Sums",                                   },
        {(void **)&ge_sum_ab,    "Addback_Sum_Energy",      "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_ab,  "AB_Sum_En_betaTagged",    "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ab_sup,"AB_Sup_Sum_Energy",      "",                           SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ab_sup_rej,"AB_Sup_Sum_En_Rejected", "",                       SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum,       "Ge_Sum_Energy",           "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b,     "Ge_Sum_En_betaTagged",    "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_sep, "Ge_Sum_En_SceptarTagged", "Ge_Sum_E_B_SEP",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_zds, "Ge_Sum_En_ZdsTagged",     "Ge_Sum_E_B_ZDS",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_art, "Ge_Sum_En_AriesTagged",   "Ge_Sum_E_B_ART",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_artT, "Ge_Sum_En_Aries_Tri",   "Ge_Sum_E_B_ART",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_artR, "Ge_Sum_En_Aries_Rec",   "Ge_Sum_E_B_ART",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_artS, "Ge_Sum_En_Aries_Sqr",   "Ge_Sum_E_B_ART",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_art_brems, "Ge_Sum_En_AriesTagged_Brems",   "Ge_Sum_E_B_ART_Brems",SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_sep_brems, "Ge_Sum_En_SceptarTagged_Brems", "Ge_Sum_E_B_SEP_Brems",SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_ab_sep_brems, "AB_Sum_En_SceptarTagged_Brems", "AB_Sum_E_B_SEP_Brems",SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&paces_sum,    "PACES_Sum_Energy",        "",                          SUBSYS_PACES,   E_SPECLEN},
        {(void **)&paces_sum_b,  "PACES_Sum_En_betaTagged", "",                          SUBSYS_PACES,   E_SPECLEN},
        {(void **)&labr_sum,     "LaBr3_Sum_Energy",        "",                          SUBSYS_LABR_L,  E_SPECLEN},
        {(void **)&aries_sum,    "ARIES_Sum_Energy",        "",                          SUBSYS_ARIES_A, E_SPECLEN},
        {(void **)&rcmp_sum,     "RCMP_Sum_Energy",         "",                          SUBSYS_RCMP,    E_SPECLEN},
        {(void **)&rcmp_fb_sum,  "RCMP_Sum_FB_Energy",      "",                          SUBSYS_RCMP,    E_SPECLEN},
        {(void **)&qed_sum,      "QED_Sum_Energy",          "",                          SUBSYS_QED_STRIP, E_SPECLEN},
        {(void **)&qed_fb_sum,   "QED_Sum_FB_Energy",       "",                          SUBSYS_QED_STRIP, E_SPECLEN},
        {(void **)&desw_sum_e,   "DES_Wall_Sum_Energy",     "",                          SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof, "DES_Wall_Sum_TOF",        "",                          SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_psd, "DES_Wall_Sum_PSD",        "",                          SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_e_b,      "DES_Wall_Sum_En_betaTagged",  "DESW_Sum_E_B",     SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof_b,    "DES_Wall_Sum_TOF_betaTagged", "DESW_Sum_TOF_B",   SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_e_nn,     "DES_Wall_Sum_En_fold2",       "DESW_Sum_E_nn",    SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof_nn,   "DES_Wall_Sum_TOF_fold2",      "DESW_Sum_TOF_nn",  SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_e_nn_a,   "DES_Wall_Sum_En_fold2_ang60", "DESW_Sum_E_nn_a",  SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof_nn_a, "DES_Wall_Sum_TOF_fold2_ang60","DESW_Sum_TOF_nn_a",SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&ge_sum_us,    "Upstream_Ge_Sum_Energy",  "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ds,    "Downstream_Ge_Sum_Energy","",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ab_us, "Upstream_AB_Sum_Energy",  "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ab_ds, "Downstream_AB_Sum_Energy","",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **) ge_ab_e,      "Addback_%d",               "",                         SUBSYS_HPGE_A,  E_SPECLEN, 0, N_CLOVER},
        {(void **) ge_ab_sup_e,  "Addback_Suppressed_%d",    "",                         SUBSYS_HPGE_A,  E_SPECLEN, 0, N_CLOVER},
        {NULL,                   "Hits_and_Sums/Energy",     "",                         },
        {(void **)&ge_xtal,      "GeEnergy_CrystalNum",      "",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&ge_xtal_1hit, "GeEnergy_CrystalNum_1_hit","",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&ge_xtal_2hit, "GeEnergy_CrystalNum_2_hit","",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&ge_xtal_3hit, "GeEnergy_CrystalNum_3_hit","",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&bgo_xtal,     "BgoEnergy_CrystalNum",     "",                         SUBSYS_BGO,     320, E_2D_SPECLEN},
        {(void **)&bgof_xtal,    "BgoFrontEnergy_CrystalNum","",                         SUBSYS_BGO,     128, E_2D_SPECLEN},
        {(void **)&bgos_xtal,    "BgoSideEnergy_CrystalNum", "",                         SUBSYS_BGO,     128, E_2D_SPECLEN},
        {(void **)&bgob_xtal,    "BgoBackEnergy_CrystalNum", "",                         SUBSYS_BGO,      64, E_2D_SPECLEN},
        {(void **)&bgoa_xtal,    "BgoAncilEnergy_CrystalNum","",                         SUBSYS_LABR_BGO, 32, E_2D_SPECLEN},
        {(void **)&labr_xtal,    "Labr3Energy_CrystalNum",   "LabrE_Xtal",               SUBSYS_LABR_L,   16, E_2D_SPECLEN},
        {(void **)&paces_xtal,   "PacesEnergy_CrystalNum",   "PacesE_Xtal",              SUBSYS_PACES,    16, E_2D_SPECLEN},
        {(void **)&sceptar_xtal, "SceptarEnergy_CrystalNum", "SceptarE_Xtal",            SUBSYS_SCEPTAR,  32, E_2D_SPECLEN},
        {(void **)&aries_xtal,   "AriesEnergy_CrystalNum",   "AriesE_Xtal",              SUBSYS_ARIES_A,  80, E_2D_SPECLEN},
        // {(void **)&labr_tac_xtal,"TAC_LBL_ART_vs_LBL_Num", "TAC_ART_LBL_LBL_Xtal",       SUBSYS_TAC_ART,   16,  E_2D_SPECLEN},
        {(void **)&art_tac_xtal, "TAC_LBL_ART_vs_ART_Num", "TAC_ART_LBL_ART_Xtal",       SUBSYS_ARIES_A,  80,  E_2D_SPECLEN},
        {(void **)&geb_xtal,      "GeBEnergy_CrystalNum",      "",                       SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&desw_e_xtal,  "DESWall_En_DetNum",    "DSW_En_Xtal",                 SUBSYS_DESWALL, 64,  E_2D_SPECLEN},
        {(void **)&desw_tof_xtal,"DESWall_TOF_DetNum",   "DSW_TOF_Xtal",                SUBSYS_DESWALL, 64,  E_2D_SPECLEN},
        {(void **)&desw_psd_e,   "DESWall_PSD_En",       "DES_Wall_PSD_En",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_tof, "DES_Wall_PSD_TOF",      "DES_Wall_PSD_TOF",            SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_q,   "DESWall_PSD_q",       "DES_Wall_PSD_q",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_cc,  "DESWall_PSD_cc",       "DES_Wall_PSD_cc",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_q_cc,  "DESWall_q_cc",       "DES_Wall_q_cc",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_q_tof,  "DESWall_q_tof",       "DES_Wall_q_tof",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_cc_tof,  "DESWall_cc_tof",       "DES_Wall_cc_tof",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_zdse, "DES_Wall_PSD_ZDSEn",      "DES_Wall_PSD_ZDSEn",       SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **) rcmp_strips,  "RCS%02d_E_strips",         "",                         SUBSYS_RCMP,   2*N_RCMP_STRIPS, E_2D_RCMP_SPECLEN, N_RCMP_POS},
        {(void **)  qed_strips,  "",                qed_strips_handles[0],               SUBSYS_QED_STRIP,   2*N_QED_STRIPS, E_2D_QED_SPECLEN, N_QED_POS},
        {NULL,                   "Hits_and_Sums/Pileup",     "",                         },
        {(void **)&gea_self_dt,  "GeA_crystal_vs_self_dt",    "",                         SUBSYS_HPGE_A,  E_2D_SPECLEN, 64},
        {(void **)&geb_self_dt,  "GeB_crystal_vs_self_dt",    "",                         SUBSYS_HPGE_A,  E_2D_SPECLEN, 64},
        {(void **)&ge_pu_class,  "Pile_up_class",           "",                          SUBSYS_HPGE_A,         64},
        {(void **) ge_sum_class,    "", ge_pu_class_sum_titles[0],                       SUBSYS_HPGE_A,  E_SPECLEN,   0, N_PU_CLASSES},
        {(void **) ge_e_vs_k_class, "", ge_pu_class_2d_titles[0],                        SUBSYS_HPGE_A,       2048, 512, N_PU_CLASSES},
        {(void **)&ge_pu_type,   "Pile_up_type",            "",                          SUBSYS_HPGE_A,         64},
        {(void **)&ge_nhits_type,"nhits_type",              "",                          SUBSYS_HPGE_A,         64},
        {(void **) ge_1hit,      "Ge%02d_Single_hit",       "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {(void **) ge_2hit,      "Ge%02d_2_hit_pileup",     "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {(void **) ge_3hit,      "Ge%02d_3_hit_pileup",     "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {NULL,                   "Hits_and_Sums/Pileup_corrections",     "",                         },
        {(void **) ge_e_vs_k_2hit_first,      "Ge%02d_E_vs_k_1st_of_2hit",              "", SUBSYS_HPGE_A,  2048, 720, 64},
        {(void **) ge_e_vs_k_2hit_second,     "Ge%02d_E_vs_k_2nd_of_2hit",              "", SUBSYS_HPGE_A,  2048, 720, 64},
        {(void **) ge_PU2_e2_v_k_gatedxrays,  "Ge%02d_PU2_E2_vs_k2_E1gated_on_Xrays",   "", SUBSYS_HPGE_A,   256, 720, 64},
        {(void **) ge_PU2_e2_v_k_gated1408,   "Ge%02d_PU2_E2_vs_k2_E1gated_on_1408keV", "", SUBSYS_HPGE_A,   256, 720, 64},
        {NULL,                   "Hits_and_Sums/BGO_HV_Alignment",     "",                         },
        {(void **) ge_bgo_gated, "",   ge_bgo_handles[0],  SUBSYS_HPGE_A,  E_2D_SPECLEN, 0, N_BGO},
        // Coinc
        {NULL,                  "Hits_and_Sums/Delta_t"," "},
        {(void **) dt_hist,     "",        dt_handles[0],  SUBSYS_HPGE_A,  DT_SPEC_LENGTH, 0, N_DT }, // leave subsys as GE -> all always defined
        {(void **) dt_tacs_hist,"dt_labr_tac%d",      "",  SUBSYS_TAC_LABR,  DT_SPEC_LENGTH, 0, N_TACS },
        {(void **) tac_lbl_ts_diff,"TAC%02d timestamp offset", "",  SUBSYS_LABR_L,  DT_SPEC_LENGTH, 0, N_TACS },
        {(void **) dcfd_hist,     "",      dcfd_handles[0],  SUBSYS_HPGE_A,  DT_SPEC_LENGTH, 0, N_DT }, // leave subsys as GE -> all always defined
        {NULL,                  "Coinc/Coinc",        ""},
        {(void **)&gg_ab,       "Addback_GG",         "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&gg,          "GG",                 "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&ge_bgo,      "GeBGO",              "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_paces,    "GePaces",            "",  SUBSYS_PACES,   E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_labr,     "GeLabr",             "",  SUBSYS_LABR_L,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_zds,      "GeZds",              "",  SUBSYS_ZDS_A,   E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_art,      "GeAries",            "",  SUBSYS_ARIES_A, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_qed,      "GeQED",            "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&qed_qed,      "QEDQED",            "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_comp,      "GeCOMP",            "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&geadd_comp,   "GeAddCOMP",         "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&comp_comp,    "COMPCOMP",          "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_dcs,      "GeDCS",              "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&geadd_dcs,   "GeAddDCS",           "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&comp_dcs,    "COMP_DCS",           "",  SUBSYS_QED_STRIP, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&labr_art,    "LaBrAries",          "",  SUBSYS_LABR_L,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&paces_art,   "PacesAries",         "",  SUBSYS_PACES,   E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&art_art,     "AriesAries",         "",  SUBSYS_ARIES_A, E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&labr_labr,   "LaBrLabr",           "",  SUBSYS_LABR_L,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&labr_zds,    "LaBrZds",            "",  SUBSYS_LABR_L,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&dsw_dsw,     "DSWDSW",             "",  SUBSYS_DESWALL,E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&ge_dsw,      "GeDSW",              "",  SUBSYS_DESWALL,E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&art_dsw,     "ARTDSW",             "",  SUBSYS_DESWALL,E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_rcmp,     "GeRCMP",             "",  SUBSYS_RCMP,    E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&labr_rcmp,   "LaBrRCMP",           "",  SUBSYS_RCMP,    E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&gg_opp,      "GGoppo",             "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&gg_ab_opp,   "GG_Addback_oppo",    "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {NULL,                  "Coinc/Hits",         ""},
        {(void **)&gg_hit,      "GeGeHit",            "",  SUBSYS_HPGE_A,  64,  64},
        {(void **)&bgobgo_hit,  "BgoBgoHit",          "",  SUBSYS_BGO,    512, 512},
        {(void **)&gea_hit,     "GeAriesHit",         "",  SUBSYS_HPGE_A,  80,  64},
        {(void **)&lba_hit,     "LaBrAriesHit",       "",  SUBSYS_LABR_L,  16,  80},
        {(void **)&aa_hit,      "AriesAriesHit",      "",  SUBSYS_ARIES_A, 80,  80},
        {(void **)&dsw_hit,     "DSWDSWHit",          "",  SUBSYS_DESWALL,64,  64},
        {(void **) rcmp_hit,    "RCS%d_PN_hit",       "",  SUBSYS_RCMP, N_RCMP_STRIPS,     N_RCMP_STRIPS,     N_RCMP_POS},
        {(void **) rcmp_fb,     "RCS%d_Front_Back",  "",   SUBSYS_RCMP, E_2D_RCMP_SPECLEN, E_2D_RCMP_SPECLEN, N_RCMP_POS},
        {(void **)&rcmp_x_ge_hit,"RCS_Xstrips_vs_GeHit","",SUBSYS_RCMP, 192,  64},
        {(void **)&rcmp_y_ge_hit,"RCS_Ystrips_vs_GeHit","",SUBSYS_RCMP, 192,  64},
        {(void **) qed_hit,    "",    qed_hit_handles[0],  SUBSYS_QED_STRIP, N_QED_STRIPS*2,     N_QED_STRIPS*2,     N_QED_POS},
        {(void **) qed_fb,     "",    qed_fb_handles[0],   SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN, N_QED_POS},
        {(void **)&qed_p_ge_hit,"QED_Pstrips_vs_GeHit","",SUBSYS_QED_STRIP, 192,  64},
        {(void **)&qed_n_ge_hit,"QED_Nstrips_vs_GeHit","",SUBSYS_QED_STRIP, 192,  64},
        //  {(void **) qed_hit_trials,    "QED%d_TRIALS",    "",  SUBSYS_QED_STRIP, 2048,     2048, N_QED_POS},
        {NULL,                  "Ang_Corr/GG_Ang_Corr",         ""},
        {(void **) gg_angcor_110,"Ge_Ge_110mm_angular_bin%02d", "", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  SYMMETERIZE, N_GE_ANG_CORR},
        {(void **) gg_angcor_145,"Ge_Ge_145mm_angular_bin%02d", "", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  SYMMETERIZE, N_GE_ANG_CORR},
        {(void **) gg_energy,    "",          gg_energy_handles[0], SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {NULL,                   "Ang_Corr/DSW_DSW_Ang_Corr",   ""},
        {(void **) dsw_angcor,   "DSW_DSW_angular_bin%03d",     "", SUBSYS_DESWALL,DSW_ANGCOR_SPECLEN, SYMMETERIZE, N_DSW_DSW_ANG_CORR },
        {NULL,                   "Ang_Corr/GG_ART_Ang_Corr",    ""},
        {(void **) ge_art_angcor,"Ge_ART_angular_bin%03d",      "", SUBSYS_ARIES_A,  GE_ANGCOR_SPECLEN, GE_ANGCOR_SPECLEN, N_GRG_ART_ANG_CORR },
        {NULL,                   "Fast_Timing/LBL_Walk",        ""},
        {(void **)&tac_labr_CompWalk0,"TAC01_LBL01_00_CompWalk",       "", SUBSYS_TAC_LABR, ECAL_TAC_SPECLEN, 1440},
        {(void **) tac_labr_CompWalk,"TAC00_LBL00_%02d_CompWalk",       "", SUBSYS_TAC_LABR, ECAL_TAC_SPECLEN, 1440, N_LABR },
        {NULL,                   "Fast_Timing/TAC_Gated_LBL_Energy", ""},
        {(void **) tac_gated_lbl,"TAC_gated_LBL%02d",           "", SUBSYS_TAC_LABR, E_SPECLEN, 0, N_LABR },
        {NULL,                   "Fast_Timing/Calibrated_TACs", ""},
        {(void **)&final_tac_sum,"Calibrated_TAC_Sum",          "", SUBSYS_TAC_LABR, E_SPECLEN },
        {(void **) final_tac,    "Calibrated_TAC%02d",          "", SUBSYS_TAC_LABR, E_SPECLEN, 0, N_LABR },
        {(void **)&lbl_lbl_tac,  "LBL_LBL_vs_TAC",              "", SUBSYS_TAC_LABR, E_3D_TAC_SPECLEN, E_3D_LBL_SPECLEN},
        {NULL,                   "Fast_Timing/ART_TACs",        ""},
        {(void **)&tac_aries_lbl_sum, "TAC_ART_LBL_LBLSUM",     "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&tac_aries_art_sum, "TAC_ART_LBL_ARTSUM",     "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&aries_tac,         "TAC_ARIES_LaBr3_1275keV","", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&aries_tac_Egate,   "TAC_ARTE_LaBr3_1275keV", "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&aries_tac_artEn,   "TAC_ARIES_Energy",       "", SUBSYS_ARIES_A, E_SPECLEN      },
        {(void **) tac_aries_lbl,    "TAC_ART_LBL%d",           "", SUBSYS_ARIES_A, E_TAC_SPECLEN, E_TAC_SPECLEN, N_LABR  },
        {(void **) tac_aries_art,    "TAC_LBL_ART%d",           "", SUBSYS_ARIES_A, E_TAC_SPECLEN, E_TAC_SPECLEN, N_ARIES },
        {NULL,                   "Analysis/Isomer_Spec",        ""},
        {(void **)&gb_dt,      "betaG_time_diff_vs_gamma_energy","", SUBSYS_HPGE_A, DT_SPEC_LENGTH, E_2D_SPECLEN},
        {(void **)&gg_dt,      "GG_time_diff_vs_gamma_energy",   "", SUBSYS_HPGE_A, DT_SPEC_LENGTH, E_2D_SPECLEN},
        {(void **)&ge_isomer_popu, "GammaE_populating",          "", SUBSYS_HPGE_A,     E_SPECLEN},
        {(void **)&ge_isomer_depop,"GammaE_depopulating",        "", SUBSYS_HPGE_A,     E_SPECLEN},
        {NULL,                   "Analysis/Cycles",        ""},
        {(void **)&cycle_num_vs_ge,       "cycle_vs_Ge",         "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_sh_g,  "cycle_vs_GeE_NPU",    "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_dt,    "cycle_vs_Ge_DT",      "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_sh,       "cycle_vs_Ge_NPU",     "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_pu,       "cycle_vs_Ge_PU",      "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_b,     "cycle_vs_GeB",        "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_b_sh_g,"cycle_vs_GeBE_NPU",   "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_b_dt,  "cycle_vs_GeB_DT",     "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_sh_b,     "cycle_vs_GeB_NPU",    "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_pu_b,     "cycle_vs_GeB_PU",     "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&ge_e_vs_cycle_time,    "Time_within_cycle_vs_Ge","",SUBSYS_HPGE_A, CYCLE_SPEC_LENGTH, E_2D_SPECLEN},
        {(void **)&ge_cycle_activity, "HPGe_cycle_activity",     "", SUBSYS_HPGE_A,     CYCLE_SPEC_LENGTH},
        {(void **)&zds_cycle_activity, "ZDS_cycle_activity",     "", SUBSYS_ZDS_A,      CYCLE_SPEC_LENGTH},
        {(void **) ge_cycle_code,    "",    ge_cycle_code_titles[0], SUBSYS_HPGE_A,     E_SPECLEN, 0, N_PPG_PATTERNS},
        {(void **) gg_cycle_code,    "",    gg_cycle_code_titles[0], SUBSYS_HPGE_A,     E_SPECLEN, E_SPECLEN, N_PPG_PATTERNS},
        {(void **) gea_cycle_num,     "HPGeA_cycle%03d",         "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_sh,  "HPGeA_NPU_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_pu,  "HPGeA_PU_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_dt,  "HPGeA_DT_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_g,   "HPGeA_GeE_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_sh_g,"HPGeA_GeE_NPU_cycle%03d", "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num,     "HPGeB_cycle%03d",         "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_sh,  "HPGeB_NPU_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_pu,  "HPGeB_PU_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_dt,  "HPGeB_DT_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_g,   "HPGeB_GeE_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_sh_g,"HPGeB_GeE_NPU_cycle%03d", "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) cycle_num_vs_geEnergy,"cycle_vs_Ge%02d_Energy",      "", SUBSYS_HPGE_A, MAX_CYCLES, E_2D_SPECLEN, N_HPGE},
        {(void **) cycle_num_vs_qedEnergy,"",    qed_E_cycle_handles[0], SUBSYS_QED_STRIP, MAX_CYCLES, E_2D_QED_SPECLEN, N_QED_POS},
        {NULL,                   "Analysis/Crosstalk",        ""},
        {(void **) ct_e_vs_dt_B,       "Crosstalk_Blue_E_vs_dt_Ge%02d",  "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {(void **) ct_e_vs_dt_G,       "Crosstalk_Green_E_vs_dt_Ge%02d", "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {(void **) ct_e_vs_dt_R,       "Crosstalk_Red_E_vs_dt_Ge%02d",   "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {(void **) ct_e_vs_dt_W,       "Crosstalk_White_E_vs_dt_Ge%02d", "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {NULL,                   "Analysis/Comp_Pol",        ""},
        {(void **) gg_comp_pol_110,"GeGe_110mm_CompPol_bin%02d", "", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  GE_ANGCOR_SPECLEN, N_GE_COMP_POL},
        {(void **) gg_comp_pol_145,"GeGe_145mm_CompPol_bin%02d", "", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  GE_ANGCOR_SPECLEN, N_GE_COMP_POL},
        {NULL,                   "QED/Calibrations",        ""},
        {(void **) qedp_ge_theta,  "",     qedp_ge_theta_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_QED_POS*N_QED_STRIPS},
        {(void **) qedn_ge_theta,  "",     qedn_ge_theta_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_QED_POS*N_QED_STRIPS},
        {(void **) qed_geE_theta_clov,    "",qed_geE_theta_clov_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_CLOVER},
        {(void **) qed_geE_theta_clov_t,  "",qed_geE_theta_clov_t_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_CLOVER},
        {NULL,                   "QED/DSSD-Ge",        ""},
        {(void **)&qedE_ge_theta_sum, "QED_E_vs_theta",         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_theta_sum, "QED_GeE_vs_theta",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        //  {(void **)&qed_totE_theta_sum,"QED_totalE_vs_theta",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **) qed_E_theta_dssd,  "",qed_E_theta_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_QED_POS},
        {(void **) qed_geE_theta_dssd,"",qed_geE_theta_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_QED_POS},
        //  {(void **) qed_totE_theta,    "",qed_totE_theta_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192, N_QED_POS},
        {(void **)&qed_E_totE_sum_t,       "QED_E_vs_totEgated",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_geE_totE_sum_t,     "QED_GeE_vs_totEgated",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qedE_ge_theta_sum_t,    "QED_E_vs_theta_totEgated",         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_theta_sum_t,    "QED_GeE_vs_theta_totEgated",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_thetaI_sum_t,   "QED_GeE_vs_thetaI_totEgated",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qedE_ge_thetaI_sum_t,   "QED_E_vs_thetaI_totEgated",         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_thetaDiff_sum_t,"QED_GeE_vs_thetaDiff_totEgated",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {NULL,                   "QED/Compton-SiGe",        ""},
        {(void **)&ge_qed_c,             "GeQED_COMP",        "",                 SUBSYS_QED_STRIP,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&qedE_ge_theta_sum_c,  "COMP_QED_E_vs_theta",         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_theta_sum_c,  "COMP_QED_GeE_vs_theta",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qedE_ge_theta_sum_c_s,"COMP_QED_E_vs_theta_Si_first",         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_theta_sum_c_s,"COMP_QED_GeE_vs_theta_Si_first",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qedE_ge_theta_sum_c_g,"COMP_QED_E_vs_theta_Ge_first",         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_geE_theta_sum_c_g,"COMP_QED_GeE_vs_theta_Ge_first",       "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&qed_dcs_omega_dt,     "QED_DCS_omega_vs_dt",     "", SUBSYS_QED_STRIP, 1024,   192},
        {(void **)&qed_dcs_omega_dtx,    "QED_DCS_omega_vs_dt_diff",    "", SUBSYS_QED_STRIP, 1024,   192},
        {(void **) qedx_dcs_omega_dt,    "",qedx_dcs_omega_dt_handles[0], SUBSYS_QED_STRIP, 1024,   192, N_QED_POS},
        {(void **)&qed_dcs_omega,        "QED_DCS_omega",           "", SUBSYS_QED_STRIP, 192},
        {(void **)&qed_dcs_omega_t,      "QED_DCS_omega_t",           "", SUBSYS_QED_STRIP, 192},
        {(void **)&qed_dcs_azi_t,        "QED_DCS_azimuth_0_180",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi,          "QED_DCS_azimuth_70_110",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_tg,       "QED_DCS_azimuth_93_103",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins1,     "QED_DCS_azimuth2_0_180",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins2,     "QED_DCS_azimuth2_10_170",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins3,     "QED_DCS_azimuth2_20_160",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins4,     "QED_DCS_azimuth2_30_150",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins5,     "QED_DCS_azimuth2_40_140",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins6,     "QED_DCS_azimuth2_50_130",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins7,     "QED_DCS_azimuth2_60_120",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins8,     "QED_DCS_azimuth2_70_110",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins9,     "QED_DCS_azimuth2_80_100",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins10,    "QED_DCS_azimuth2_85_95",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_bins8a,    "QED_DCS_azimuth2_93_103",          "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_omega_dt_TRWF,      "QED_DCS_omega_vs_dt_TRWF",           "", SUBSYS_QED_STRIP, 4096, 192},
        {(void **)&qed_dcs_azi_TRWF_t,        "QED_DCS_azimuth_TRWF_0_180",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF,          "QED_DCS_azimuth_TRWF_70_110",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_tg,       "QED_DCS_azimuth_TRWF_93_103",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins1,     "QED_DCS_azimuth2_TRWF_0_180",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins2,     "QED_DCS_azimuth2_TRWF_10_170",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins3,     "QED_DCS_azimuth2_TRWF_20_160",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins4,     "QED_DCS_azimuth2_TRWF_30_150",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins5,     "QED_DCS_azimuth2_TRWF_40_140",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins6,     "QED_DCS_azimuth2_TRWF_50_130",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins7,     "QED_DCS_azimuth2_TRWF_60_120",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins8,     "QED_DCS_azimuth2_TRWF_70_110",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins9,     "QED_DCS_azimuth2_TRWF_80_100",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins10,    "QED_DCS_azimuth2_TRWF_85_95",         "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_dcs_azi_TRWF_bins8a,    "QED_DCS_azimuth2_TRWF_93_103",        "", SUBSYS_QED_STRIP, 384},
        {(void **)&qed_theta,            "COMP_QED_theta",           "",SUBSYS_QED_STRIP, 192},
        {(void **)&qed_wf_omega,         "COMP_QED_weight_omega",        "",SUBSYS_QED_STRIP, 192},
        {(void **)&qed_wf_dcs_azi,       "COMP_QED_weight_dcs_azi",     "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_theta1_vs_theta2, "COMP_QED_theta1_vs_theta2",  "",SUBSYS_QED_STRIP, 192, 192},
        {(void **)&qed_delta_theta1_theta2,"COMP_QED_delta_theta1_theta2",  "",SUBSYS_QED_STRIP, 192},
        {(void **)&qed_sum_theta1_theta2,"COMP_QED_sum_theta1_theta2",  "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_theta1_azi,       "COMP_QED_theta1_vs_azi",  "",SUBSYS_QED_STRIP, 192, 384},
        {(void **)&qed_theta2_azi,       "COMP_QED_theta2_vs_azi",  "",SUBSYS_QED_STRIP, 192, 384},
        {(void **)&qed2_theta1_vs_theta2,"COMP_QED2_theta1_vs_theta2",  "",SUBSYS_QED_STRIP, 192, 192},
        {(void **)&qed2_theta1_azi,      "COMP_QED2_theta1_vs_azi",  "",SUBSYS_QED_STRIP, 192, 384},
        {(void **)&qed2_theta2_azi,      "COMP_QED2_theta2_vs_azi",  "",SUBSYS_QED_STRIP, 192, 384},
        {(void **)&qed_ge_weight,        "COMP_QED_GE_weights",      "",SUBSYS_QED_STRIP, 6144, 64},
        {NULL,                   "QED/Dbl-COMPTON-SiGeGe",        ""},
        {(void **)&dcsa_theta,         "DCSA_theta",                                   "",SUBSYS_QED_STRIP, 192},
        {(void **)&dcsaE_ge_theta,     "DCSA_E_vs_ICStheta",                           "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&dcsa_geE_theta,     "DCSA_GeE_vs_ICStheta",                         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&dcsa_cs_omega_ge,   "DCSA_omega_GeGe_SiGeGe",                       "",SUBSYS_QED_STRIP, 192},
        {(void **)&dcsa_cs_omega,      "DCSA_omega_SiGe-SiGeGe",                       "",SUBSYS_QED_STRIP, 192},
        {(void **)&dcsa_theta_azi_ge,  "DCSA_ICStheta_vs_azimuth_GeGe_SiGeGe",         "",SUBSYS_QED_STRIP, 384, 384},
        {(void **)&dcsa_theta_azi,     "DCSA_ICStheta_vs_azimuth_SiGe_SiGeGe",         "",SUBSYS_QED_STRIP, 384, 384},
        {NULL,                   "QED/Dbl-COMPTON-SiSiGe",        ""},
        {(void **)&dcsbE_ge_theta,     "DCSB_E_vs_ICStheta",                           "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&dcsb_geE_theta,     "DCSB_GeE_vs_ICStheta",                         "",SUBSYS_QED_STRIP, E_2D_QED_SPECLEN,   192},
        {(void **)&dcsb_cs_omega_ge,   "DCSB_omega_GeGe-SiSiGe",                       "",SUBSYS_QED_STRIP, 192},
        {(void **)&dcsb_cs_omega,      "DCSB_omega_SiGe-SiSiGe",                       "",SUBSYS_QED_STRIP, 192},
        {(void **)&dcsb_theta_azi_ge,  "DCSB_ICStheta_vs_azimuth_GeGe_SiSiGe",         "",SUBSYS_QED_STRIP, 384, 384},
        {(void **)&dcsb_theta_azi,     "DCSB_ICStheta_vs_azimuth_SiGe_SiSiGe",         "",SUBSYS_QED_STRIP, 384, 384},
        //  {NULL,                   "QED/PSD",        ""},
        //  {(void **) qed_psd_e,      "",           qed_psd_handles[0],SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_SPECLEN, N_QED_POS},
        {NULL,                   "QED/Triples",        ""},
        {(void **)&qed_qed_23,  "QED2_QED3",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_23dt,  "QED2_QED3_dt",       "",     SUBSYS_QED_STRIP, DT_SPEC_LENGTH},
        {(void **)&qed_qed_23_theta2,  "QED2_QED3_E2_theta",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, 192},
        {(void **)&qed_qed_23_theta3,  "QED2_QED3_E3_theta",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, 192},
        {(void **)&qed_qed_23_totv2,  "QED2_QED3_E2_totE",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_23_totv3,  "QED2_QED3_E3_totE",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_12,  "QED1_QED2",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_12dt,  "QED1_QED2_dt",       "",     SUBSYS_QED_STRIP, DT_SPEC_LENGTH},
        {(void **)&qed_qed_12_theta1,  "QED1_QED2_E1_theta",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, 192},
        {(void **)&qed_qed_12_theta2,  "QED1_QED2_E2_theta",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, 192},
        {(void **)&qed_qed_12_totv1,  "QED1_QED2_E1_totE",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_12_totv2,  "QED1_QED2_E2_totE",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_14,  "QED1_QED4",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_14dt,  "QED1_QED4_dt",       "",     SUBSYS_QED_STRIP, DT_SPEC_LENGTH},
        {(void **)&qed_qed_14_theta1,  "QED1_QED4_E1_theta",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, 192},
        {(void **)&qed_qed_14_theta4,  "QED1_QED4_E4_theta",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, 192},
        {(void **)&qed_qed_14_totv1,  "QED1_QED4_E1_totE",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {(void **)&qed_qed_14_totv4,  "QED1_QED4_E4_totE",       "",     SUBSYS_QED_STRIP, E_2D_QED_SPECLEN, E_2D_QED_SPECLEN},
        {NULL,                   "QED/Misc.",        ""},
        {(void **)&qed_angle_phi_s,  "QED_angle_phi_Si_first",         "",SUBSYS_QED_STRIP, N_HPGE,   384},
        {(void **)&qed_angle_phi_a,  "QED_angle_phi_70_110",         "",SUBSYS_QED_STRIP, N_HPGE,   384},
        {(void **)&qed_angle_phi_b,  "QED_angle_phi_93_103",         "",SUBSYS_QED_STRIP, N_HPGE,   384},
        {(void **)&qed_phi_s,  "QED_phi_Si_first",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_a,  "QED_phi_70_110",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_b,  "QED_phi_93_103",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins1,  "QED_phi_weight_0_180",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins2,  "QED_phi_weight_10_170",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins3,  "QED_phi_weight_20_160",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins4,  "QED_phi_weight_30_150",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins5,  "QED_phi_weight_40_140",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins6,  "QED_phi_weight_50_130",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins7,  "QED_phi_weight_60_120",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins8,  "QED_phi_weight_70_110",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins9,  "QED_phi_weight_80_100",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qed_phi_bins10,  "QED_phi_weight_85_95",         "",SUBSYS_QED_STRIP, 384},
        {(void **)&qedE_ge_dt,  "QED_Ge_E_vs_dt",       "",SUBSYS_QED_STRIP, 1024, E_2D_QED_SPECLEN},
        {(void **)&qed_geE_dt,  "QED_GeE_vs_dt",        "",SUBSYS_QED_STRIP, 1024, E_2D_QED_SPECLEN},
        {(void **)&qed_theta_dt,  "QED_theta_vs_dt",        "",SUBSYS_QED_STRIP, 1024, 192},
        {(void **)&qed_theta_dt_cfd,  "QED_theta_vs_dt_cfd",        "",SUBSYS_QED_STRIP, 2048, 192},
        {(void **)&qedE_ge_dt_c,  "QED_COMP_E_vs_dt",        "",SUBSYS_QED_STRIP, 1024, E_2D_QED_SPECLEN},
        {(void **)&qed_geE_dt_c,  "QED_COMP_GeE_vs_dt",        "",SUBSYS_QED_STRIP, 1024, E_2D_QED_SPECLEN},
        {(void **)&qed_theta_dt_c,  "QED_COMP_theta_vs_dt",        "",SUBSYS_QED_STRIP, 1024, 192},
      }; // Note initialized array variable is CONST (not same as double-pointer)
      // TH1I *hist;  hist = (TH1I *) 0;   ptr = &hist = (TH1I **)addr;  *ptr =

      void get_histo_handle(char *result, char *title, char *handle )
      {
        if( strlen(handle) != 0 ){ memcpy(result, handle, strlen(handle)); return; }
        while( *title != 0 ){
          if( strcmp(title, "Energy")     == 0){ memcpy(result, "E",    1); title +=  6; result += 1; continue; }
          if( strcmp(title, "CrystalNum") == 0){ memcpy(result, "Xtal", 4); title += 10; result += 4; continue; }
          if( strcmp(title, "Upstream")   == 0){ memcpy(result, "US",   2); title +=  8; result += 2; continue; }
          if( strcmp(title, "Downstream") == 0){ memcpy(result, "DS",   2); title += 10; result += 2; continue; }
          if( strcmp(title, "Des_Wall")   == 0){ memcpy(result, "DESW", 4); title +=  8; result += 4; continue; }
          *result = *title; ++title; ++result;
        }
        *result = 0; return;
      }

      // only create histograms when that subsystem is seen in the data
      // this function is called everytime a new sybsystem is seen
      int init_histos(Config *cfg, int subsystem)
      {
        Histogram_definition *hptr;
        static Config *save_cfg;
        char tmp[STRING_LEN], **tptr;
        int i, j, k;

        if( cfg == NULL ){ cfg = save_cfg; } else { save_cfg = cfg; }
        subsys_initialized[subsystem] = 1;
        for(i=0; i<HISTO_DEF_SIZE; i++){ hptr = &histodef_array[i];
          if( hptr->ptr == NULL ){ // new folder or empty definition
            if( strlen(hptr->title) != 0 ){ memcpy(cfg->current_path, hptr->title, strlen(hptr->title)+1 ); }
            continue;
          }
          if( subsystem != hptr->subsys ){ continue; } // skip
          if( hptr->count == 0 ){ // single histogram
            get_histo_handle(tmp, hptr->title, (char *)hptr->handle );
            if( hptr->ychan == 0 ){ // 1d
              *(TH1I **)(hptr->ptr) = H1_BOOK(cfg, tmp, hptr->title, hptr->xchan, 0, hptr->xchan );
            } else {
              *(TH2I **)(hptr->ptr) = H2_BOOK(cfg, tmp, hptr->title, hptr->xchan, 0, hptr->xchan,
                hptr->ychan, 0, hptr->ychan );
              }
            } else { // array of histos (title and handle both the same)
              for(j=0; j<hptr->count; j++){
                sprintf(tmp, (strlen(hptr->title) == 0) ? hptr->handle+j*HANDLE_LENGTH : hptr->title, j);
                if( hptr->ychan == 0 ){ // 1d
                  *(TH1I **)(hptr->ptr+j) = H1_BOOK(cfg, tmp, tmp, hptr->xchan, 0, hptr->xchan );
                } else {
                  *(TH2I **)(hptr->ptr+j) = H2_BOOK(cfg, tmp, tmp, hptr->xchan, 0, hptr->xchan,
                    hptr->ychan, 0, hptr->ychan );
                  }
                }
              }
            } *cfg->current_path=0; // empty path

            // custom definitions that don't fit usual scheme
            if( subsystem == SUBSYS_TAC_LABR ){ // TAC coincidence pair spectra
              open_folder(cfg, "Fast_Timing");
              open_folder(cfg, "LBL_TAC_Combos");
              k=0; memset(tac_labr_hist_index, -1, N_LABR*N_LABR*sizeof(int));
              for(i=0; i<N_LABR; i++){
                for(j=(i+1); j<N_LABR; j++){
                  tac_labr_hist_index[i][j] = k;
                  sprintf(tmp,"TAC_%02d_%02d", i, j);
                  tac_labr_hist[k++] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
                }
              }
              sprintf(tmp,"TAC_%02d_%02d", 1, 0); // Add additional histogram (2_1) needed for Compton Walk corrections
              tac_labr_hist_index[1][0] = k;
              tac_labr_hist[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
              // and the uncalibrated versions
              k=0;
              for(i=0; i<N_LABR; i++){
                for(j=(i+1); j<N_LABR; j++){
                  sprintf(tmp,"uncalibrated_TAC_%02d_%02d", i, j);
                  tac_labr_hist_uncal[k++] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
                }
              }
              sprintf(tmp,"uncalibrated_TAC_%02d_%02d", 1, 0); // Add additional histogram (2_1) needed for Compton Walk corrections
              tac_labr_hist_uncal[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
              close_folder(cfg);
              close_folder(cfg);
            }
            // fill in subsys EvsE and Dt table pointers [** [X][<=Y]
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_HPGE_A ] = gg;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_BGO    ] = ge_bgo;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_PACES  ] = ge_paces;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_LABR_L ] = ge_labr;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_RCMP   ] = ge_rcmp;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_ZDS_A  ] = ge_zds;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_QED_PIXEL  ] = ge_qed;
            subsys_e_vs_e[SUBSYS_QED_PIXEL ][SUBSYS_QED_PIXEL  ] = qed_qed;
            subsys_e_vs_e[SUBSYS_PACES  ][SUBSYS_ARIES_A] = paces_art;
            subsys_e_vs_e[SUBSYS_LABR_L ][SUBSYS_LABR_L ] = labr_labr;
            subsys_e_vs_e[SUBSYS_LABR_L ][SUBSYS_ARIES_A] = labr_art;
            subsys_e_vs_e[SUBSYS_LABR_L ][SUBSYS_ZDS_A  ] = labr_zds;
            subsys_e_vs_e[SUBSYS_ARIES_A][SUBSYS_ARIES_A] = art_art;
            // Timestamp differences in 10ns units
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_HPGE_A  ] = dt_hist[ 0];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_PACES   ] = dt_hist[ 4];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_LABR_L  ] = dt_hist[ 5];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_RCMP    ] = dt_hist[ 6];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_ZDS_A   ] = dt_hist[ 3];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_ARIES_A ] = dt_hist[10];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_BGO     ] = dt_hist[ 1];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_SCEPTAR ] = dt_hist[ 2];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_DESWALL ] = dt_hist[19];
            subsys_dt[SUBSYS_PACES  ][SUBSYS_LABR_L  ] = dt_hist[ 8];
            subsys_dt[SUBSYS_PACES  ][SUBSYS_ARIES_A ] = dt_hist[12];
            subsys_dt[SUBSYS_PACES  ][SUBSYS_ZDS_A   ] = dt_hist[ 7];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_LABR_L  ] = dt_hist[26];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_ARIES_A ] = dt_hist[11];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_ZDS_A   ] = dt_hist[17];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_TAC_LABR] = dt_hist[16];
            subsys_dt[SUBSYS_RCMP   ][SUBSYS_RCMP    ] = dt_hist[ 9];
            subsys_dt[SUBSYS_ARIES_A][SUBSYS_ARIES_A ] = dt_hist[13];
            subsys_dt[SUBSYS_ARIES_A][SUBSYS_TAC_ART ] = dt_hist[14];
            subsys_dt[SUBSYS_ZDS_A  ][SUBSYS_TAC_ZDS ] = dt_hist[15];
            subsys_dt[SUBSYS_ZDS_A  ][SUBSYS_ZDS_B   ] = dt_hist[22];
            subsys_dt[SUBSYS_DESWALL][SUBSYS_DESWALL ] = dt_hist[18];
            subsys_dt[SUBSYS_DESWALL][SUBSYS_ARIES_A ] = dt_hist[20];
            subsys_dt[SUBSYS_DESWALL][SUBSYS_ZDS_A   ] = dt_hist[21];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_QED_PIXEL  ] = dt_hist[27];
            subsys_dt[SUBSYS_QED_PIXEL][SUBSYS_QED_PIXEL] = dt_hist[28];
            subsys_dt[SUBSYS_COMPTON][SUBSYS_COMPTON] = dt_hist[29];
            subsys_dt[SUBSYS_HPGE_A][SUBSYS_COMPTON] = dt_hist[30];
            // CFD differences in 10ns units
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_HPGE_A  ] = dcfd_hist[ 0];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_PACES   ] = dcfd_hist[ 4];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_LABR_L  ] = dcfd_hist[ 5];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_RCMP    ] = dcfd_hist[ 6];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_ZDS_A   ] = dcfd_hist[ 3];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_ARIES_A ] = dcfd_hist[10];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_BGO     ] = dcfd_hist[ 1];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_SCEPTAR ] = dcfd_hist[ 2];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_DESWALL ] = dcfd_hist[19];
            subsys_dcfd[SUBSYS_PACES  ][SUBSYS_LABR_L  ] = dcfd_hist[ 8];
            subsys_dcfd[SUBSYS_PACES  ][SUBSYS_ARIES_A ] = dcfd_hist[12];
            subsys_dcfd[SUBSYS_PACES  ][SUBSYS_ZDS_A   ] = dcfd_hist[ 7];
            subsys_dcfd[SUBSYS_LABR_L ][SUBSYS_LABR_L  ] = dcfd_hist[26];
            subsys_dcfd[SUBSYS_LABR_L ][SUBSYS_ARIES_A ] = dcfd_hist[11];
            subsys_dcfd[SUBSYS_LABR_L ][SUBSYS_ZDS_A   ] = dcfd_hist[17];
            subsys_dcfd[SUBSYS_LABR_L ][SUBSYS_TAC_LABR] = dcfd_hist[16];
            subsys_dcfd[SUBSYS_RCMP   ][SUBSYS_RCMP    ] = dcfd_hist[ 9];
            subsys_dcfd[SUBSYS_ARIES_A][SUBSYS_ARIES_A ] = dcfd_hist[13];
            subsys_dcfd[SUBSYS_ARIES_A][SUBSYS_TAC_ART ] = dcfd_hist[14];
            subsys_dcfd[SUBSYS_ZDS_A  ][SUBSYS_TAC_ZDS ] = dcfd_hist[15];
            subsys_dcfd[SUBSYS_ZDS_A  ][SUBSYS_ZDS_B   ] = dcfd_hist[22];
            subsys_dcfd[SUBSYS_DESWALL][SUBSYS_DESWALL ] = dcfd_hist[18];
            subsys_dcfd[SUBSYS_DESWALL][SUBSYS_ARIES_A ] = dcfd_hist[20];
            subsys_dcfd[SUBSYS_DESWALL][SUBSYS_ZDS_A   ] = dcfd_hist[21];
            subsys_dcfd[SUBSYS_HPGE_A ][SUBSYS_QED_PIXEL  ] = dcfd_hist[27];
            subsys_dcfd[SUBSYS_QED_PIXEL][SUBSYS_QED_PIXEL] = dcfd_hist[28];
            subsys_dcfd[SUBSYS_COMPTON][SUBSYS_COMPTON] = dcfd_hist[29];
            subsys_dcfd[SUBSYS_HPGE_A][SUBSYS_COMPTON] = dcfd_hist[30];

            return(0);
          }


          int init_parameters_from_globals(Config *cfg){
            // Here set parameters from the Globals
            // time-difference conditions for subsys-subsys coincidences
            // PRESORT time-difference Conditions
            // TAC energy offset parameters based on each LBL-LBL combination
            int i,j,k,c1,c2,index;
            Global *global;
            char tmp[32];
            const char *ptr;

            // Initialize all time differences between subsystems to be the default 250ns
            for(i=0; i<MAX_SUBSYS; i++){
              for(j=0; j<MAX_SUBSYS; j++){
                time_diff_gate_min[i][j] = 0;  // default is 0 nanoseconds
                time_diff_gate_max[i][j] = 25; // default is 250 nanoseconds
              }
            }

            // Intialize all PRESORT timing windows to their defaults
            bgo_window_min = addback_window_min = rcmp_fb_window_min = qed_fb_window_min = lbl_tac_window_min = zds_tac_window_min = art_tac_window_min = desw_beta_window_min = 0;
            bgo_window_max = 20;
            addback_window_max = 20;
            rcmp_fb_window_max = 10;
            qed_fb_window_max  = 10;
            lbl_tac_window_max = 25;
            art_tac_window_max = 25;
            zds_tac_window_max = 25;
            desw_beta_window_max = 50;

            // Initalize the cycles gamma-ray energy gate values
            ppg_cycles_gamma_gate_min = 1800;
            ppg_cycles_gamma_gate_max = 1820;
            ppg_cycles_binning_factor = 10000000; // Default of 10,000,000 converts 10ns to 100 millisecond binning

            // Search the globals for time difference settings and overwrite their values
            for(i=0; i<cfg->nglobal; i++){
              global = cfg->globals[i];
              sprintf(tmp,"%s",global->name);
              if(strncmp(tmp,"time_diff_",10) == 0){
                // This global is a time difference value
                // Identify the subsystem types and then save the value in the correct place
                for(j=0; j<MAX_SUBSYS; j++){
                  if(strlen(subsys_handle[j])<2){ continue; }
                  if((ptr = strstr(tmp,subsys_handle[j])) > 0){
                    ptr += strlen(subsys_handle[j]);
                    // Identify the second subsystem type
                    for(k=0; k<MAX_SUBSYS; k++){
                      if(strlen(subsys_handle[k])<2){ continue; }
                      if(strstr(ptr,subsys_handle[k]) > 0){
                        // save the value in the correct place
                        fprintf(stdout,"time_diff_%s_%s [%d,%d] set to %d,%d\n",subsys_handle[j],subsys_handle[k],j,k,global->min,global->max);
                        time_diff_gate_min[j][k] = global->min;
                        time_diff_gate_max[j][k] = global->max;
                        time_diff_gate_min[k][j] = global->min;
                        time_diff_gate_max[k][j] = global->max;
                        break;
                      }
                    } break;
                  }
                }
              }else if(strncmp(tmp,"TAC-Offset-",11) == 0){
                // TAC energy offset parameters based on each LBL-LBL combination
                // Only use the maximum value, ignore the minimum
                // Derive the index from the name
                sscanf(tmp,"TAC-Offset-%02d-%02d",&c1,&c2);
                c1--; c2--;
                if(c1>=0 && c1<N_LABR && c2>=0 && c2<N_LABR){
                  index = tac_labr_hist_index[c1][c2];
                  tac_lbl_combo_offset[index] = global->max;
                }else{
                  fprintf(stderr,"Problem decoding TAC Offset from Global, %s\n",tmp);
                }
              }else if(strncmp(tmp,"Timestamp-offset-LBT",20) == 0){
                // LBT timestamp offset value
                // Derive the index from the name
                sscanf(tmp,"Timestamp-offset-LBT%02d",&c1);
                c1--;
                if(c1>=0 && c1<N_TACS){
                  tac_ts_offset[c1] = global->max;
                }else{
                  fprintf(stderr,"Problem decoding LBT timestamp offset from Global, %s\n",tmp);
                }
              }else if(strncmp(tmp,"presort_time_diff_addback",25) == 0){
                addback_window_min = global->min; addback_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_zds_tac",25) == 0){
                zds_tac_window_min = global->min; zds_tac_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_art_tac",25) == 0){
                art_tac_window_min = global->min; art_tac_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_lbl_tac",25) == 0){
                lbl_tac_window_min = global->min; lbl_tac_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_desw_beta",27) == 0){
                desw_beta_window_min = global->min; desw_beta_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_suppression",29) == 0){
                bgo_window_min = global->min; bgo_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_qed_front-back",32) == 0){
                qed_fb_window_min = global->min; qed_fb_window_max = global->max;
              }else if(strncmp(tmp,"presort_time_diff_rcmp_front-back",33) == 0){
                rcmp_fb_window_min = global->min; rcmp_fb_window_max = global->max;
              }else if(strncmp(tmp,"cycles_bin_size_in_ms",21) == 0){
                ppg_cycles_binning_factor = (int)(global->max * 100000); // milliseconds to 10ns timestamp unit conversion
              }else if(strncmp(tmp,"cycles_gamma_gate",17) == 0){
                ppg_cycles_gamma_gate_min = global->min; ppg_cycles_gamma_gate_min = global->max;
              }else if(strncmp(tmp,"presort_window_width",20) == 0){
                presort_window_width = global->max;
              }else if(strncmp(tmp,"sort_window_width",17) == 0){
                sort_window_width = global->max;
              }// end of for(i=0; i<cfg->nglobal; i++){
              }
              return(0);
            }


            //#######################################################################
            //###########   READ XML ODB DUMP FROM START OF DATA FILE   #############
            //#######################################################################

            // Note - the odb names do not distinguish between subtypes of detectors
            // e.g Ge A and B channels
            // the subsystem names will be extended to include this information
            // (and the odb-specific names are only used below)

            #define MAX_ODB_SUBSYS 24
            #define ODBHANDLE_GRG   0
            #define ODBHANDLE_GRS   1
            #define ODBHANDLE_SEP   2
            #define ODBHANDLE_PAC   3
            #define ODBHANDLE_LBS   4
            #define ODBHANDLE_LBT   5
            #define ODBHANDLE_LBL   6
            #define ODBHANDLE_DSC   7
            #define ODBHANDLE_ART   8
            #define ODBHANDLE_ZDS   9
            #define ODBHANDLE_RCS  10
            #define ODBHANDLE_XXX  11
            #define ODBHANDLE_DSW  12
            #define ODBHANDLE_DSG  13
            #define ODBHANDLE_DAL  14
            #define ODBHANDLE_DAT  15
            #define ODBHANDLE_QED  16
            #define ODBHANDLE_UNK  23
            static char odb_handle[MAX_ODB_SUBSYS][8] = {
              "GRG", "GRS", "SEP", "PAC",  //  0- 3
              "LBS", "LBT", "LBL", "DSC",  //  4- 7
              "ART", "ZDS", "RCS", "XXX",  //  8-11
              "DSW", "DSG", "DAL", "DAT",  //  12-15
              "QED",    "",    "",    "",
              "",    "",    "",    "UNK"
            };

            static char   path[256];
            static char dirname[64],value[32],type[32];
            extern char midas_runtitle[SYS_PATH_LENGTH];

            static void *arrayptr;
            int read_odb_items(int len, int *bank_data)
            {
              char *path_ptr, *ptr, *str, *odb_data = (char *)bank_data, posn[2], odb_ppg_current[128];
              int i, j, c = '<', d = '>', dtype=0, active=0, index=0, ppg_index, odb_ppg_cycle_count=-1;

              // The currently set cycle is included in the ODB dump after all cycles are defined.
              // So we need to unpack them all and then select the values we need.
              ppg_cycles odb_ppg_cycle[MAX_ODB_PPG_CYCLES];

              ptr = odb_data;  path_ptr = path;
              while(1){
                if( (str = strchr(ptr,c)) == NULL ){ break; }
                ptr = str;
                if( (str = strchr(ptr,d)) == NULL ){ break; }

                if( strncmp(ptr,"<!--",4) == 0 || strncmp(ptr,"<odb", 4) == 0 ||
                strncmp(ptr,"</odb",5) == 0 ){ // comment - skip
                } else if( strncmp(ptr,"<dir ",5) == 0 ){
                  if( strncmp(ptr,"<dir name=\"",11) == 0 ){
                    i=11; while( *(ptr+i) != '"' && *(ptr+i) != d ){ ++i; }
                  }
                  memcpy(dirname, ptr+11, i-11); dirname[i-11] = '\0';
                  if( *(ptr+1+i) == '/' ){ ptr=str+1; continue; }
                  //if( sscanf(ptr,"<dir name=\"%s\">", dirname) < 1 ){
                  //   fprintf(stderr,"can't read dirname\n"); ptr=str+1; continue;
                  //}
                  //if( strncmp(dirname+strlen(dirname)-3,"\"/>",3) == 0 ){
                  //   ptr=str+1; continue;
                  //}
                  //if( dirname[strlen(dirname)-1]=='>'  ){
                  //   dirname[strlen(dirname)-1]='\0';
                  //}
                  //if( dirname[strlen(dirname)-1]=='\"' ){
                  //  dirname[strlen(dirname)-1]='\0';
                  //}
                  *path_ptr = '/'; strcpy(path_ptr+1, dirname);
                  path_ptr += strlen(dirname)+1;
                  *path_ptr = '\0';
                  if( strncmp(path,"/PPG/Cycles/",12) == 0 ){
                    odb_ppg_cycle_count++;
                    strcpy(odb_ppg_cycle[odb_ppg_cycle_count].name, dirname);
                  }
                } else if( strncmp(ptr,"</dir>",6) == 0 ){
                  while(1){
                    if( --path_ptr < path ){ path_ptr = path;  *path_ptr = '\0';  break; }
                    if( *path_ptr == '/' ){ *path_ptr = '\0';  break; }
                  }
                  index=0; // for debugger to stop here
                } else if( strncasecmp(ptr,"<key name=\"Run Title\" type=\"STRING\"", 35) == 0 ){
                  ptr = str+1;
                  if( (str = strchr(ptr,c)) == NULL ){ break; }
                  i = (str-ptr) > SYS_PATH_LENGTH-1 ? SYS_PATH_LENGTH-1 : (str-ptr);
                  memcpy( midas_runtitle, ptr, i ); midas_runtitle[i] = 0;
                  ptr += i+1;
                  if( (str = strchr(ptr,d)) == NULL ){ break; }
                } else if( strncasecmp(ptr,"<key name=\"Current\"", 19) == 0 &&
                strncmp(path,"/PPG",4) == 0 ){
                  *str = ' '; if( (str = strchr(str,c)) == NULL ){ break; }
                  *str = ' '; if( (str = strchr(str,d)) == NULL ){ break; }
                  if( sscanf(ptr,"<key name=\"Current\" type=\"STRING\" size=\"%d\" %s /key>", &dtype, value) < 2 ){
                    fprintf(stderr,"can't read key name for Current PPG cycle\n"); ptr=str+1; continue;
                  }
                  strcpy(odb_ppg_current,value);
                  ptr += i+1;
                } else if( strncasecmp(ptr,"<key name=\"prg_ddtm\"", 20) == 0 &&
                strncmp(path,"/DAQ/params/grif16/template/0",29) == 0 ){
                  if( sscanf(ptr,"<key name=\"prg_ddtm\" type=\"DWORD\">%d</key>", &subsys_prg_ddtm[SUBSYS_HPGE_A]) < 1 ){
                    fprintf(stderr,"can't read key value for /DAQ/params/grif16/template/0/prg_ddtm\n"); ptr=str+1; continue;
                  }
                  fprintf(stdout,"Read in Det type 0 (HPGE A) prg_ddtm as %d\n",subsys_prg_ddtm[SUBSYS_HPGE_A]);
                  while( *(ptr) != '/' ){ ++ptr; } while( *(ptr) != '<' ){ ++ptr; }
                } else if( strncasecmp(ptr,"<key name=\"prg_ddtm\"", 20) == 0 &&
                strncmp(path,"/DAQ/params/grif16/template/1",29) == 0 ){
                  if( sscanf(ptr,"<key name=\"prg_ddtm\" type=\"DWORD\">%d</key>", &subsys_prg_ddtm[SUBSYS_HPGE_B]) < 1 ){
                    fprintf(stderr,"can't read key value for /DAQ/params/grif16/template/1/prg_ddtm\n"); ptr=str+1; continue;
                  }
                  fprintf(stdout,"Read in Det type 1 (HPGE B) prg_ddtm as %d\n",subsys_prg_ddtm[SUBSYS_HPGE_B]);
                  while( *(ptr) != '/' ){ ++ptr; } while( *(ptr) != '<' ){ ++ptr; }
                } else if( strncmp(ptr,"</keyarray>",10) == 0 ){ active = 0; arrayptr = (void *)('\0');
                if( strncmp(path,"/PPG/Cycles/",12) == 0 ){ odb_ppg_cycle[odb_ppg_cycle_count].length = index+1; }
              } else if( strncmp(ptr,"<keyarray ",10) == 0 ){
                if( strncmp(path,"/PPG/Cycles/",12) == 0 ){
                  if( sscanf(ptr,"<keyarray name=\"%s", value) < 1 ){
                    fprintf(stderr,"can't read PPG keyarray entry\n"); ptr=str+1; continue;
                  }
                  if( value[strlen(value)-1]=='\"' ){ value[strlen(value)-1]='\0'; }
                  if( strcmp(value,"PPGcodes") == 0 ){
                    active = 1; arrayptr = (void *)odb_ppg_cycle[odb_ppg_cycle_count].codes; dtype=0;
                  }
                  if( strcmp(value,"durations") == 0 ){
                    active = 1; arrayptr = (void *)odb_ppg_cycle[odb_ppg_cycle_count].durations; dtype=0;
                  }
                  ptr=str+1;
                  continue;
                }
                if( strcmp(path,"/DAQ/params/MSC") != 0 &&
                strcmp(path,"/DAQ/MSC")        != 0 &&
                strcmp(path,"/DAQ/PSC")        != 0 ){  ptr=str+1; continue; }
                if( sscanf(ptr,"<keyarray name=\"%s", value) < 1 ){
                  fprintf(stderr,"can't read keyarray entry\n"); ptr=str+1; continue;
                }
                if( value[strlen(value)-1]=='\"' ){ value[strlen(value)-1]='\0'; }
                if( strcmp(value,"PSC") == 0 || strcmp(value,"MSC") == 0 ){
                  active = 1; arrayptr = (void *)addr_table; dtype=1;
                }
                if( strcmp(value,"chan") == 0 ){
                  active = 1; arrayptr = (void *)chan_name; dtype=3;
                }
                if( strcmp(value,"datatype") == 0 ){
                  //active = 1; arrayptr = (void *)dtype_table; dtype=1;
                  active = 1; arrayptr = (void *)dtype_table; dtype=0;
                }
                if( strcmp(value,"gain") == 0 ){
                  active = 1; arrayptr = (void *)gain_table; dtype=2;
                }
                if( strcmp(value,"offset") == 0 ){
                  active = 1; arrayptr = (void *)offs_table; dtype=2;
                }
                if( strcmp(value,"quadratic") == 0 ){
                  active = 1; arrayptr = (void *)quad_table; dtype=2;
                }
              } else if( strncmp(ptr,"<value index=",13) == 0 ){
                if( !active ){ ptr=str+1; continue; }
                // remove the >< surrounding the value, and move str to the end of the line
                *str = ' '; if( (str = strchr(str,c)) == NULL ){ break; }
                *str = ' '; if( (str = strchr(str,d)) == NULL ){ break; }
                if( sscanf(ptr,"<value index=\"%d\" %s /value>", &index, value) < 2 ){
                  fprintf(stderr,"can't read value entry\n");
                }
                if( index < 0 || index >= MAX_DAQSIZE ){
                  fprintf(stderr,"index %d out of range\n", index);
                }
                // index starts at zero, odb_daqsize is count
                if( index >= odb_daqsize ){ odb_daqsize = index+1; }
                if(        dtype == 0 ){  // int
                  if( sscanf(value,"%d", (((int *)arrayptr)+index)) < 1 ){
                    fprintf(stderr,"can't read value %s\n", value);
                  }
                } else if( dtype == 1 ){  // short int
                  if( sscanf(value,"%hd", (((short *)arrayptr)+index)) < 1 ){
                    fprintf(stderr,"can't read value %s\n", value);
                  }
                } else if( dtype == 2 ){  // float
                  if( sscanf(value,"%f", (((float *)arrayptr)+index)) < 1 ){
                    fprintf(stderr,"can't read value %s\n", value);
                  }
                } else {                 // string
                  strncpy(arrayptr+index*CHAN_NAMELEN, value, CHAN_NAMELEN);
                  *((char *)arrayptr+(index+1)*CHAN_NAMELEN - 1) = '\0';
                }
              }
              ptr=str+1;
            }
            fprintf(stdout,"odb record: %d bytes\n", len);

            // PPG: Now all cycles were unpacked and we identified Current
            // Copy the relevant ODB PPG pattern into the global variables for this run
            index=-1;
            for(i=0; i<MAX_ODB_PPG_CYCLES; i++){
              if( strncmp(odb_ppg_cycle[i].name,odb_ppg_current,strlen(odb_ppg_current)) == 0 ){
                index=i;
                break;
              }
            }
            if(index<0){
              fprintf(stderr,"Failed to locate Current PPG Cycle in ODB Cycles\n");
            }else{
              strcpy(ppg_cycle_name,odb_ppg_cycle[index].name);
              fprintf(stdout,"PPG cycle for this run is named %s:\n",ppg_cycle_name);
              ppg_cycle_duration=0;
              for(i=0; i<odb_ppg_cycle[index].length; i++){
                ppg_index=-1;
                for(j=0; j<N_PPG_PATTERNS; j++){ if( (odb_ppg_cycle[index].codes[i] & 0xFFFF) == ppg_patterns[j] ){ ppg_index = j; break; } }
                if(ppg_index<0){
                  fprintf(stderr,"unrecognized ppg pattern, 0x%04X\n", (odb_ppg_cycle[index].codes[i] & 0xFFFF));
                  gen_derived_odb_tables();
                  return(-1);
                }
                ppg_cycle_pattern_code[i] = ppg_index;

                if(odb_ppg_cycle[index].durations[i] == -1){
                  // Infinte duration
                  odb_ppg_cycle[index].length = i+1;
                }else{
                  ppg_cycles_active = 1;
                  ppg_cycle_length = odb_ppg_cycle[index].length;
                  ppg_cycle_pattern_duration[i] = (long)odb_ppg_cycle[index].durations[i]*100; // Convert from ODB microseconds to timestamp 10 nanosecond units
                  ppg_cycle_duration += ppg_cycle_pattern_duration[i];
                }
                fprintf(stdout,"PPG PATTERN %d: %s (%s) for %10.4f milliseconds (%015ld timestamps)\n", i,
                ppg_handles[ppg_cycle_pattern_code[i]], ppg_names[ppg_cycle_pattern_code[i]],(double)(ppg_cycle_pattern_duration[i]/100000),ppg_cycle_pattern_duration[i]);
              }
              if(ppg_cycle_duration == 0){
                fprintf(stdout,"PPG cycle duration is infinite, ie. no cycles\n");
                fprintf(stdout,"Setting PPG cycle duration to 15 seconds for diagnostics\n");
                ppg_cycle_duration = 1500000000; // 15 seconds
                ppg_cycle_pattern_duration[0] = 1500000000; // 15 seconds
                //  fprintf(stdout,"Setting PPG cycle duration to 5 minutes for diagnostics\n");
                //  ppg_cycle_duration = 30000000000; // 5 minutes
                //  ppg_cycle_pattern_duration[0] = 30000000000; // 5 minutes
                ppg_cycles_active = 1;
                ppg_cycle_length = odb_ppg_cycle[index].length;
              }else{
                fprintf(stdout,"PPG cycle duration is %10.4f seconds\n",(double)(ppg_cycle_duration/100000000));
              }
              // Set the initial cycle settings
              ppg_current_pattern = ppg_cycle_pattern_code[0];  // Index of the current PPG cycle pattern for use with the ppg_patterns array
              ppg_cycle_number = 0;                             // Current cycle number. Cycles counted from zero at beginning of run
              ppg_cycle_start = 0;                              // Timestamp of the start of the current cycle
              ppg_cycle_end = ppg_cycle_duration;               // Timestamp of the end of the current cycle
              ppg_cycle_step = 0;                               // Current pattern number within this cycle. Patterns counted from zero at beginning of cycle
              ppg_pattern_start = 0;                            // Timestamp of the start of the current pattern
              ppg_pattern_end = ppg_cycle_pattern_duration[0];  // Timestamp of the end of the current pattern
              fprintf(stdout,"Cycle %04d, start/finish [%ld/%ld]: step %d, %s, start/finish [%ld/%ld], ppg_current_pattern=%d\n",
              ppg_cycle_number, ppg_cycle_start, ppg_cycle_end, ppg_cycle_step, ppg_handles[ppg_current_pattern], ppg_pattern_start, ppg_pattern_end, ppg_current_pattern);
            }

            // arrays typically around 500 entries [one per "chan"] each entry with ...
            //   daq-address, name, type, gains etc.
            //
            gen_derived_odb_tables();

            return(0);
          }

          // This lookup table reorders strips that have already been reordered in the ODB...
          // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25966
          // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25968
          static int reorder_rcmp_strips[7][2][32] = {
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS1 X
            {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // RCS1 Y
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS2 X
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // RCS2 X
            {{1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30},  // RCS3 X
            {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // RCS3 Y
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS4 X
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // RCS4 Y
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS5 X
            {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // RCS5 Y
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS6 X
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}  // RCS6 Y
          };

          // This lookup table reorders strips that have already been reordered in the ODB...
          // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25966
          // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25968
          static int reorder_qed_strips[7][2][32] = {
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // QED1 P
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // QED1 N

            {{1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30},  // QED2 P
            {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // QED2 N

            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // QED3 P
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // QED3 N

            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // QED4 P
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // QED4

            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // QED5 P
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // QED5 N

            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // QED6 P
            {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}  // QED6 N
          };

          // original odb arrays were read into {addr_table,chan_name,dtype_table(+gains)}
          // extract extra details stored in channel names (and record for later)
          // (these details include crystal/element numbers and polarities)
          // [use above for subsystem (no longer use datatype to determine subsystems)]
          extern int read_caen_odb_addresses(int odb_daqsize, unsigned short *addr_table);
          int gen_derived_odb_tables()
          {
            int i, j, tmp, subsys, pos, element, output_type;
            char sys_name[64], crystal, polarity, type;

            read_caen_odb_addresses(odb_daqsize, (unsigned short *)addrs);

            // generate reverse mapping of address to channel number
            //  (most of this array is undefined and stays at -1)
            memset(address_chan, 0xFF, sizeof(address_chan)); // set to -1
            for(i=0; i<MAX_ADDRESS && i<odb_daqsize; i++){
              address_chan[ (unsigned short)chan_address[i] ] = i;
            }

            memset(crystal_table,  0xff, MAX_DAQSIZE*sizeof(int)); // initialise all to -1
            memset(element_table,  0xff, MAX_DAQSIZE*sizeof(int));
            memset(polarity_table, 0xff, MAX_DAQSIZE*sizeof(int));
            memset(subsys_table,   0xff, MAX_DAQSIZE*sizeof(int));
            for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
              if( (tmp=sscanf(chan_name[i], "%3c%d%c%c%d%c", sys_name, &pos, &crystal, &polarity, &element, &type)) != 6 ){
                fprintf(stderr,"can't decode name[%s] decoded %d of 6 items\n", chan_name[i], tmp );
                continue;
              }
              for(j=0; j<MAX_ODB_SUBSYS; j++){
                if( strncmp(sys_name, odb_handle[j], 3) == 0 ){ subsys = j; break; }
              }
              if( j == MAX_ODB_SUBSYS ){ subsys = j-1; // use final entry: "unknown"
              fprintf(stderr,"Unknown subsystem[%s] in %s\n", sys_name, chan_name[i]);
            }

            // Mention bad detector types (no longer relied on for subsystem id)
            if( dtype_table[i] < 0 || dtype_table[i] >= 16 ){
              fprintf(stderr,"bad datatype[%d] at table position %d\n", dtype_table[i], i);
            }

            // Some detector elements have more than one output (HPGe A and B)
            // 1 is A, 0 is B, -1 is X or unknown
            output_type = type=='A' ? 1 : (type=='B' ? 0 : -1);
            if(subsys == ODBHANDLE_ZDS){ if(type=='X'){ output_type = 1; } } // Older ZDS convention

            // Polarity: 1 is N, 0 is P or T or S, -1 is anything else
            if(        polarity == 'N' ){ polarity_table[i] = 1;
            } else if( polarity == 'P' ){ polarity_table[i] = 0;
            } else if( polarity == 'T' ){ polarity_table[i] = 0; // TAC signal
            } else if( polarity == 'S' ){ polarity_table[i] = 1; // ARIES Standard Ouput signal
            } else if( polarity == 'F' ){ polarity_table[i] = 0; // ARIES Fast Output signal
            } else if( polarity == 'X' ){ polarity_table[i] = 0; // XXX type
            } else { fprintf(stderr,"unknown polarity[=%c] in %s\n", polarity, chan_name[i]); }

            // Record crystal and element numbers [** Naming schemes are subsystem-dependant **]
            switch(subsys){
              case ODBHANDLE_LBL: case ODBHANDLE_LBS: // LaBr,Paces, Aries and Zds
              case ODBHANDLE_LBT: case ODBHANDLE_SEP: case ODBHANDLE_ART:
              case ODBHANDLE_DAL: case ODBHANDLE_DAT:
              case ODBHANDLE_PAC: case ODBHANDLE_ZDS: case ODBHANDLE_DSW:
              crystal_table[i] = pos;
              if(        crystal == 'A' ){ element_table[i] = 1;
              } else if( crystal == 'B' ){ element_table[i] = 2;
              } else if( crystal == 'C' ){ element_table[i] = 3;
              } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3, ZDS, ART, LBT, SEP
              } else {
                fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
              } break;
              case ODBHANDLE_RCS:
              crystal_table[i] = pos;
              element_table[i] = reorder_rcmp_strips[pos][polarity_table[i]][element];
              break;
              case ODBHANDLE_QED:
              crystal_table[i] = pos;
              element_table[i] = reorder_qed_strips[pos][polarity_table[i]][element];
              break;
              case ODBHANDLE_GRG: case ODBHANDLE_GRS:
              element_table[i] = element;
              pos -= 1; pos *=4;
              if(        crystal == 'B' ){ crystal_table[i] = pos;
              } else if( crystal == 'G' ){ crystal_table[i] = pos+1;
              } else if( crystal == 'R' ){ crystal_table[i] = pos+2;
              } else if( crystal == 'W' ){ crystal_table[i] = pos+3;
              } else if( crystal == 'X' ){ crystal_table[i] = -1; // crystal undefined
              } else {
                fprintf(stderr,"unknown crystal[=%c] in %s\n", crystal, chan_name[i]);
              } break;
              default: break;
            }

            // set full subsystem id (including polarity/output-type etc)
            switch(subsys){
              case ODBHANDLE_GRS: subsys_table[i] = SUBSYS_BGO;       break;
              case ODBHANDLE_SEP: subsys_table[i] = SUBSYS_SCEPTAR;   break;
              case ODBHANDLE_PAC: subsys_table[i] = SUBSYS_PACES;     break;
              case ODBHANDLE_LBS: subsys_table[i] = SUBSYS_LABR_BGO;  break;
              case ODBHANDLE_LBL: subsys_table[i] = SUBSYS_LABR_L;    break;
              case ODBHANDLE_DAL: subsys_table[i] = SUBSYS_LABR_L;    break;
              case ODBHANDLE_DSC: subsys_table[i] = SUBSYS_DESCANT;   break;
              case ODBHANDLE_RCS: subsys_table[i] = SUBSYS_RCMP;      break;
              case ODBHANDLE_QED: subsys_table[i] = SUBSYS_QED_STRIP; break;
              case ODBHANDLE_DSW: subsys_table[i] = SUBSYS_DESWALL;  break;
              case ODBHANDLE_DSG: subsys_table[i] = SUBSYS_DSG;  break;
              case ODBHANDLE_GRG: subsys_table[i] = (output_type == 1) ? SUBSYS_HPGE_A :SUBSYS_HPGE_B; break;
              case ODBHANDLE_ZDS: subsys_table[i] = (output_type == 1) ? SUBSYS_ZDS_A  :SUBSYS_ZDS_B;  break;
              case ODBHANDLE_ART: subsys_table[i] = (polarity_table[i] == 1) ? SUBSYS_ARIES_A:SUBSYS_ARIES_B;break;
              case ODBHANDLE_XXX: subsys_table[i] = SUBSYS_IGNORE;    break;
              case ODBHANDLE_UNK: subsys_table[i] = SUBSYS_UNKNOWN;   break;
              case ODBHANDLE_DAT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
              }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
              break;
              case ODBHANDLE_LBT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
              }else if(crystal_table[i]>8){ subsys_table[i] = SUBSYS_TAC_ART;
              }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
              break;
            }
          }
          memset(subsys_initialized, 0, sizeof(int)*MAX_SUBSYS );

          return(0);
        }
