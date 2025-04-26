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
#define E_2D_TOF_SPECLEN        1024
#define E_2D_SPECLEN            4096
#define E_2D_RCMP_SPECLEN       6400
#define E_3D_LBL_SPECLEN      160000  // 400*400
#define E_3D_TAC_SPECLEN        1024
//#define T_SPEC_LENGTH     8192
//#define WV_SPEC_LENGTH    4096
#define DT_SPEC_LENGTH          1024
#define GE_ANGCOR_SPECLEN       4096
#define DSW_ANGCOR_SPECLEN      4096

//#######################################################################
//########        Individual channel singles HISTOGRAMS        ##########
//#######################################################################

// Pulse height and energy
TH1I   *ts_hist; // timestamp
TH1I   *gc_hist; // GRIF-CAEN hitpattern for checking coincidences
TH1I   *ph_hist[MAX_DAQSIZE];
TH1I    *e_hist[MAX_DAQSIZE];
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

// DESCANT WALL
TH1I  *desw_sum_e, *desw_sum_tof, *desw_sum_psd;  // Sums of energies and corTOF and PSD
TH1I  *desw_sum_e_b, *desw_sum_tof_b;       // Beta-tagged Sums of energies and corTOF
TH1I  *desw_sum_e_nn, *desw_sum_tof_nn;     // fold>2 Sums of energies and corTOF
TH1I  *desw_sum_e_nn_a, *desw_sum_tof_nn_a; // fold>2, angle>60 Sums of energies and corTOF
TH2I  *desw_psd_e, *desw_psd_tof;           // PSD vs energies or corrected-TOF

// TAC spectra
TH1I *tac_labr_hist[(int)((N_LABR)*(N_LABR-1)/2)+2]; // this index numbers are the LaBr-LaBr position numbers
TH1I *tac_labr_hist_uncal[(int)((N_LABR)*(N_LABR-1)/2)+2]; // this index numbers are the LaBr-LaBr position numbers
// One additional histogram (2_1) needed for Compton Walk corrections
TH2I *tac_labr_CompWalk[N_LABR];         // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
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
  TH2I *gg, *gg_ab, *gg_opp, *ge_paces, *ge_labr, *ge_rcmp, *labr_labr, *labr_zds, *labr_rcmp;
  TH2I *ge_art, *ge_zds, *paces_art, *labr_art, *art_art, *dsw_dsw, *ge_dsw, *art_dsw;

  // Angular Correlation histograms
  #define N_GE_ANG_CORR       52
  #define N_GRG_ART_ANG_CORR 114
  #define N_DSW_DSW_ANG_CORR  42
  TH2I  *gg_angcor_110[N_GE_ANG_CORR];
  TH2I  *gg_angcor_145[N_GE_ANG_CORR];
  TH2I  *ge_art_angcor[N_GRG_ART_ANG_CORR];
  TH2I  *dsw_angcor[N_DSW_DSW_ANG_CORR];

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
