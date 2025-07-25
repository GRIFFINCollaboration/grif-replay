// default_sort.h
// Header file for default_sort.c


//#######################################################################
//########         Subsystem and Detector definitions          ##########
//#######################################################################

// do not alter order without also changing subsys_e_vs_e, subsys_dt
#define MAX_SUBSYS      24
#define SUBSYS_HPGE_A    0
#define SUBSYS_PACES     1
#define SUBSYS_LABR_L    2
#define SUBSYS_RCMP      3
#define SUBSYS_ARIES_A   4
#define SUBSYS_ZDS_A     5 // GRIF16
#define SUBSYS_TAC_LABR  6
#define SUBSYS_LABR_BGO  7
#define SUBSYS_BGO       8
#define SUBSYS_SCEPTAR   9
#define SUBSYS_DESCANT  10
#define SUBSYS_DESWALL  11
#define SUBSYS_DSG      12
#define SUBSYS_IGNORE   13
#define SUBSYS_HPGE_B   16
#define SUBSYS_ARIES_B  17 // CAEN
#define SUBSYS_ZDS_B    18 // CAEN
#define SUBSYS_TAC_ZDS  19
#define SUBSYS_TAC_ART  20
#define SUBSYS_UNKNOWN  23
static char subsys_handle[MAX_SUBSYS][8] = {
  "GRGA", "PAC",  "LBL",  "RCS",
  "ARTA", "ZDSA", "LBT",  "LBS",
  "BGO",  "SEP",  "DSC",  "DSW",
  "DSG", "XXX1", "XXX2", "XXX3",
  "GRGB", "ARTB", "ZDSB", "", // secondary names start after #16
  "",     "",     "",     "UNK"
};
static char subsys_name[MAX_SUBSYS][STRING_LEN] = {
  "Griffin", "PACES",   "LaBrX",   "RCMP",     //  0- 3
  "ARIES",   "ZDSA",    "TAC_LBL",   "LaBrS",    //  4- 7
  "BGO",     "Sceptar", "Descant", "DES_WALL", //  8-11
  "Des_Ancil", "Ignore1", "Ignore2", "Ignore3",  // 12-15
  "Grif_B",  "ARS_B",   "ZDS_B",   "TAC_ZDS",      // 16-19
  "TAC_ART",      "",        "",        "Unknown"   // 20-23
}; // final entry will be used if not found - make sure it is not empty
// #####################################################################

#define N_CLOVER 16
#define N_HPGE 64
#define N_BGO 320
#define N_ARIES 76
#define N_LABR 8
#define N_TACS 12
#define N_RCMP_POS 6
#define N_RCMP_STRIPS 32
#define N_DES_WALL 60

//#######################################################################
//########                Histogram axis lengths               ##########
//#######################################################################

#define MULT_SPEC_LENGTH         128
#define E_SPECLEN               8192
#define E_TAC_SPECLEN          16384
#define ECAL_TAC_SPECLEN        1024
#define E_TOF_SPEC_LENGTH       8192
#define E_PSD_SPEC_LENGTH       1024
#define CYCLE_SPEC_LENGTH       1024  // At 100 millisecond binning this supports 17 minute cycles (1024 seconds)
#define E_2D_TOF_SPECLEN        1024
#define E_2D_SPECLEN            4096
#define E_2D_RCMP_SPECLEN       6400
#define E_3D_LBL_SPECLEN      160000  // 400*400
#define E_3D_TAC_SPECLEN         512
//#define T_SPEC_LENGTH     8192
//#define WV_SPEC_LENGTH    4096
#define DT_SPEC_LENGTH          1024
#define GE_ANGCOR_SPECLEN       4096
#define DSW_ANGCOR_SPECLEN      4096

//#######################################################################
//########             PPG variables and patterns              ##########
//#######################################################################

// PPG patterns and handles
#define N_PPG_PATTERNS  7
int ppg_patterns[N_PPG_PATTERNS]={ 0xC008,0xC002,0xC001,0xC004,0xC0F0,0xC009,0xC00A };
char ppg_handles[N_PPG_PATTERNS][32]={ "0xC008","0xC002","0xC001","0xC004","0xC0F0","0xC009","0xC00A" };
char ppg_names[N_PPG_PATTERNS][32]={
  "Move_Tape", "Background", "Beam_on_Implant", "Beam_off_Decay", "Source_data", "Continuous_Tape_Beam_on", "Continuous_Tape_Background"
};
#define MAX_ODB_PPG_CYCLES 50
typedef struct ppg_cycles_struct {
   char name[128];  int length; int codes[16]; int durations[16];
} ppg_cycles;

// Definitions for the Current cycle of this run
// These are derived from the ODB settings at BOR
long ppg_cycle_duration;             // Length of one cycle in timestamp units
char ppg_cycle_name[128];            // Name of this current cycle
long ppg_cycle_length;               // Number of patterns/durations for this current cycle
long ppg_cycle_pattern_duration[16]; // Length of each pattern in timestamp units
int  ppg_cycle_pattern_code[16];     // Index of each pattern for use with the ppg_patterns array
int  ppg_cycles_active;              // Cycles active or made inactive if set to Source/constant beam-on etc.

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

// Spectra for cycles
// The binning factor and gamma-energy gates will ultimately be set from a Global at BOR
#define MAX_CYCLES 512                          // Used as an axis length for some 2d histograms so must be a multiple of 16
long ppg_cycles_binning_factor = 10000000;      // Default of 10,000,000 converts 10ns to 100 millisecond binning
int ppg_cycles_gamma_gate_min = 1800;           // 26Na, S1140
int ppg_cycles_gamma_gate_max = 1820;           // 26Na, S1140
char ge_cycle_code_titles[N_PPG_PATTERNS][HANDLE_LENGTH]={ "0xC008_Move_Tape","0xC002_Background","0xC001_Beam_on_Implant","0xC004_Beam_off_Decay","0xC0F0_Source_data","0xC009_Continuous_Tape_Beam_on", "0xC00A_Continuous_Tape_Background" };
char gg_cycle_code_titles[N_PPG_PATTERNS][HANDLE_LENGTH]={ "0xC008_GG_Move_Tape","0xC002_GG_Background","0xC001_GG_Beam_on_Implant","0xC004_GG_Beam_off_Decay","0xC0F0_GG_Source_data","0xC009_GG_Continuous_Tape_Beam_on", "0xC00A_GG_Continuous_Tape_Background" };
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

//#######################################################################
//########        Individual channel singles HISTOGRAMS        ##########
//#######################################################################

// Pulse height and energy
TH1I   *ts_hist; // timestamp
TH1I   *gc_hist; // GRIF-CAEN hitpattern for checking coincidences
TH1I   *ph_hist[MAX_DAQSIZE];
TH1I   *e_hist[MAX_DAQSIZE];
//TH1I *wave_hist[MAX_DAQSIZE];

// Hitpatterns
#define N_HITPAT  7
char hit_handles[N_HITPAT][32]={ "q_hit","e_hit","t_hit","w_hit","r_hit", "s_hit", "d_hit" };
char   hit_names[N_HITPAT][32]={
  "Pulse_Height", "Energy", "Time", "Waveform", "Rate", "Subsys", "DetType",
};
TH1I  *hit_hist[N_HITPAT], *mult_hist[MAX_SUBSYS];

// DESCANT Wall
TH1I  *desw_tof[N_DES_WALL];                // Time-Of-Flight
TH1I  *desw_tof_corr[N_DES_WALL];           // corrected Time-Of-Flight
TH1I  *desw_tof_psd[N_DES_WALL];            // corrected Time-Of-Flight, PSD gated
TH1I  *desw_psd[N_DES_WALL];                // Pulse Shape Discrimination

//#######################################################################
//########                PRESORT Time Gates                   ##########
//#######################################################################

// The definition of the time difference gate in 10 nanosecond units.
// The value is the maximum time difference in 10 nanosecond units.
// The default values set here are replaced by the Global value at start of sorting.
static int bgo_window_min = 0;
static int addback_window_min = 0;
static int rcmp_fb_window_min = 0;
static int lbl_tac_window_min = 0;
static int art_tac_window_min = 0;
static int zds_tac_window_min = 0;
static int desw_beta_window_min = 0;
static int bgo_window_max = 20;
static int addback_window_max = 20;
static int rcmp_fb_window_max = 10;
static int lbl_tac_window_max = 25;
static int art_tac_window_max = 25;
static int zds_tac_window_max = 25;
static int desw_beta_window_max = 80;

//#######################################################################
//########                Coincidence Time Gates               ##########
//#######################################################################

// The definition of the time difference gate in 10 nanosecond units.
// First and second index are the subsystem index numbers
// The value is the maximum time difference in 10 nanosecond units.
// Default is 250 nanoseconds, replaced by the Global value at start of sorting
static int time_diff_gate_min[MAX_SUBSYS][MAX_SUBSYS];
static int time_diff_gate_max[MAX_SUBSYS][MAX_SUBSYS];

//#######################################################################
//########          Sums and coincidence  HISTOGRAMS           ##########
//#######################################################################

// HPGe pileup
#define N_PU_CLASSES 15
// Pileup Class definitions
#define PU_ERROR             0
#define PU_SINGLE_HIT        1
#define PU_SINGLE_HIT_ERROR  2
#define PU_2HIT_A1ST         3
#define PU_2HIT_A2ND         4
#define PU_2HIT_B1ST         5
#define PU_2HIT_B2ND         6
#define PU_2HIT_C1ST         7
#define PU_2HIT_C2ND         8
#define PU_2HIT_ERROR        9
#define PU_3HIT_1ST         10
#define PU_3HIT_2ND         11
#define PU_3HIT_3RD         12
#define PU_3HIT_ERROR       13
#define PU_OTHER            14
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
TH1I *subsys_dt[MAX_SUBSYS][MAX_SUBSYS];
TH1I *tac_lbl_ts_diff[N_TACS];

// HPGe (ge_sum is sum of crystal energies, ge_sum_b is beta-gated)
TH1I  *ge_ab_e[N_CLOVER], *ge_sum_ab;
TH1I  *ge_sum, *ge_sum_us, *ge_sum_ds, *ge_sum_ab_us, *ge_sum_ab_ds;
TH1I  *ge_sum_b, *ge_sum_b, *ge_sum_b_sep, *ge_sum_b_zds, *ge_sum_b_art, *ge_sum_b_art_brems;

// ARIES, PACES and LaBr3
TH1I  *aries_sum;  // aries_sum is sum of tile energies
TH1I  *paces_sum;  // paces_sum is sum of crystal energies
TH1I  *labr_sum;  // labr_sum is sum of crystal energies

// RCMP
TH1I  *rcmp_sum, *rcmp_fb_sum;  // rcmp_sum is sum of strip energies, fb is with front-back coincidence
TH2I  *rcmp_strips[N_RCMP_POS];
TH2I  *rcmp_hit[N_RCMP_POS];
TH2I  *rcmp_fb[N_RCMP_POS];
TH2I  *rcmp_x_ge_hit, *rcmp_y_ge_hit; // rcmp strips vs Ge hitpatterns

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
int tac_labr_hist_index[N_LABR][N_LABR]; // index for filling tac_labr_hist from LBL id numbers
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
TH2I *ge_xtal, *bgo_xtal, *bgof_xtal, *bgos_xtal, *bgob_xtal, *bgoa_xtal, *labr_xtal;
TH2I *labr_tac_xtal, *paces_xtal, *aries_xtal, *art_tac_xtal, *desw_e_xtal, *desw_tof_xtal;

#define N_DT 27   // Time difference
char dt_handles[N_DT][HANDLE_LENGTH]={
  "dt_ge_ge",      "dt_ge_bgo",     "dt_ge_sep",             "dt_ge_zds",     // 0-3
  "dt_ge_pac",     "dt_ge_labr",    "dt_ge_rcmp",            "dt_pac_zds",    // 4-7
  "dt_pac_labr",   "dt_rcmp_rcmp",  "dt_ge_art",             "dt_labr_art",   // 8-11
  "dt_paces_art",  "dt_art_art",    "dt_art_tac",            "dt_zds_tac",    // 12-15
  "dt_labr_tac",   "dt_labr_zds",   "dt_dsw_dsw",            "dt_dsw_ge",     // 16-19
  "dt_dsw_art",    "dt_dsw_zds",    "dt_zds_GRIF_CAEN_10ns", "dt_zds_GRIF_CAEN_2ns", // 20-23
  "dt_dsw_dsw_2ns","dt_dsw_zds_2ns", "dt_labr_labr"  };                                      // 24-26
  TH1I  *dt_hist[N_DT], *dt_tacs_hist[N_TACS];

  // 2D hitpatterns
  TH2I *gg_hit, *bgobgo_hit, *aa_hit, *gea_hit, *lba_hit, *dsw_hit;

  // 2d Energy vs Energy Coincidence matrices
  TH2I *gg, *gg_ab, *gg_opp, *gg_ab_opp, *ge_bgo, *ge_paces, *ge_labr, *ge_rcmp, *labr_labr, *labr_zds, *labr_rcmp;
  TH2I *ge_art, *ge_zds, *paces_art, *labr_art, *art_art, *dsw_dsw, *ge_dsw, *art_dsw;

  // Angular Correlation histograms
  #define N_GE_ANG_CORR       52
  #define N_GRG_ART_ANG_CORR 114
  #define N_DSW_DSW_ANG_CORR  42
  TH2I  *gg_angcor_110[N_GE_ANG_CORR];
  TH2I  *gg_angcor_145[N_GE_ANG_CORR];
  TH2I  *ge_art_angcor[N_GRG_ART_ANG_CORR];
  TH2I  *dsw_angcor[N_DSW_DSW_ANG_CORR];

  // Isomer Spectroscopy
  TH2I  *gg_dt, *gb_dt;
  TH1I  *ge_isomer_popu, *ge_isomer_depop;

    // Crosstalk Analysis
    TH2I  *ct_e_vs_dt_B[N_HPGE], *ct_e_vs_dt_G[N_HPGE], *ct_e_vs_dt_R[N_HPGE], *ct_e_vs_dt_W[N_HPGE];

  ////////////////////////////////////
  ////////////////////////////////////

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

  // BGO HV alignment histograms
  TH1I *ge_bgo_gated[N_BGO];
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
