// default_sort.h
// Header file for default_sort.c


//#######################################################################
//########         Subsystem and Detector definitions          ##########
//#######################################################################

#define SUBSYS_HPGE      0
#define SUBSYS_BGO       1
#define SUBSYS_SCEPTAR   2
#define SUBSYS_PACES     3
#define SUBSYS_LABR_BGO  4
#define SUBSYS_LABR_T    5
#define SUBSYS_LABR_L    6
#define SUBSYS_DESCANT   7
#define SUBSYS_ARIES     8
#define SUBSYS_ZDS       9
#define SUBSYS_RCMP     10
#define SUBSYS_DES_WALL 12

#define MAX_SUBSYS 24
static char subsys_handle[MAX_SUBSYS][8] = {
  "GRG", "GRS", "SEP",  "PAC",
  "LBS", "LBT", "LBL",  "DSC",
  "ART", "ZDS", "RCS",  "XXX",
  "DSW",    "",    "",    "",
  "",    "",    "",    "",
  "",    "",    "",    ""
};
static char subsys_name[MAX_SUBSYS][STRING_LEN] = {
   "Griffin",  "BGO",   "SCEPTAR",   "PACES", //  0- 3
   "LaBrS",    "LaBrT", "LaBrX",   "Descant", //  4- 7
   "ARIES",    "ZDS",   "RCMP",        "XXX", //  8-11
   "DES_WALL",    "",      "",         "",    // 12-15
   "",         "",      "",         "",
   "",         "",       "",        "Unknown"
}; // final entry will be used if not found - make sure it is not empty

#define N_CLOVER 16
#define N_ARIES 76
#define N_LABR 8
#define N_RCMP_POS 6
#define N_RCMP_STRIPS 32
#define N_DES_WALL 60

//#######################################################################
//########                Histogram axis lengths               ##########
//#######################################################################

#define MULT_SPEC_LENGTH         128
#define E_SPEC_LENGTH           8192
#define E_TAC_SPEC_LENGTH      16384
#define E_TOF_SPEC_LENGTH       8192
#define E_PSD_SPEC_LENGTH       1024
#define E_2D_TOF_SPEC_LENGTH    1024
#define E_2D_SPEC_LENGTH        4096
#define E_2D_RCMP_SPEC_LENGTH   6400
//#define T_SPEC_LENGTH     8192
//#define WV_SPEC_LENGTH    4096
#define DT_SPEC_LENGTH          1024
#define GE_ANG_CORR_SPEC_LENGTH 4096
#define DSW_ANG_CORR_SPEC_LENGTH 4096
#define SYMMETERIZE 0 // When the ybins are set to zero for a 2D histogram it will be symmeterized

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


//#######################################################################
//########          Sums and coincidence  HISTOGRAMS           ##########
//#######################################################################

// HPGe
TH1I  *ge_ab_e[N_CLOVER];
TH1I  *ge_sum, *ge_sum_ab;  // ge_sum is sum of crystal energies
TH1I  *ge_sum_us, *ge_sum_ds, *ge_sum_ab_us, *ge_sum_ab_ds;  // ge_sum is sum of crystal energies
TH1I  *ge_sum_b; // beta-gated gamma sum spectrum
TH1I  *ge_sum_b, *ge_sum_b_sep, *ge_sum_b_zds, *ge_sum_b_art; // beta-gated gamma sum spectrum per ancillary

// ARIES, PACES and LaBr3
TH1I  *aries_sum;  // aries_sum is sum of tile energies
TH1I  *paces_sum;  // paces_sum is sum of crystal energies
TH1I  *labr_sum;  // labr_sum is sum of crystal energies

// RCMP
TH1I  *rcmp_sum, *rcmp_fb_sum;  // rcmp_sum is sum of strip energies, fb is with front-back coincidence
char rcmp_strips_handles[N_RCMP_POS][32]={ "RCS01_E_strips","RCS02_E_strips","RCS03_E_strips","RCS04_E_strips","RCS05_E_strips","RCS06_E_strips" };
TH2I  *rcmp_strips[N_RCMP_POS];
char rcmp_hit_handles[N_RCMP_POS][32]={ "RCS01_PN_hits","RCS02_PN_hits","RCS03_PN_hits","RCS04_PN_hits","RCS05_PN_hits","RCS06_PN_hits" };
char rcmp_fb_handles[N_RCMP_POS][32]={ "RCS01_Front_Back","RCS02_Front_Back","RCS03_Front_Back","RCS04_Front_Back","RCS05_Front_Back","RCS06_Front_Back" }; // front-back
TH2I  *rcmp_hit[N_RCMP_POS];
TH2I  *rcmp_fb[N_RCMP_POS];

// DESCANT WALL
TH1I  *desw_tof[N_DES_WALL];  // DESCANT Wall Time-Of-Flight
TH1I  *desw_tof_corr[N_DES_WALL];  // DESCANT Wall corrected Time-Of-Flight
TH1I  *desw_tof_psd[N_DES_WALL];  // DESCANT Wall corrected Time-Of-Flight, PSD gated
TH1I  *desw_psd[N_DES_WALL];  // DESCANT Wall Pulse Shape Discrimination
TH1I  *desw_sum_e, *desw_sum_tof, *desw_sum_psd;  // DESCANT Wall Sums of energies and corrected time-of-fligts and PSD
TH1I  *desw_sum_e_b, *desw_sum_tof_b;  // DESCANT Wall Beta-tagged Sums of energies and corrected time-of-fligts
TH1I  *desw_sum_e_nn, *desw_sum_tof_nn;  // DESCANT Wall fold>2 Sums of energies and corrected time-of-fligts
TH1I  *desw_sum_e_nn_a, *desw_sum_tof_nn_a;  // DESCANT Wall fold>2, angle>60 Sums of energies and corrected time-of-fligts
TH2I  *desw_psd_e, *desw_psd_tof; // DESCANT Wall PSD vs energies or corrected-TOF

// TAC spectra
TH1I *tac_labr_hist[(int)((N_LABR)*(N_LABR-1)/2)]; // this index numbers are the LaBr-LaBr position numbers
TH1I *tac_aries_lbl_hist[N_LABR];  // this index number is the LaBr position number
TH1I *tac_aries_art_hist[N_ARIES];  // this index number is the Aries position number
TH1I *tac_aries_lbl_sum;  // ARIES TAC sum spectrum of all LBLs
TH1I *tac_aries_art_sum;  // ARIES TAC sum spectrum of all ARTs
TH1I *aries_tac;  // aries_tac gated on 1275keV peak
TH1I *aries_tac_Egate;  // aries_tac gated on 1275keV peak
TH1I *aries_tac_artEn;  // aries energy in coincidence with TAC
TH2I *lblE_tac, *zdsE_tac, *ariesE_tac;  // lbl or zds or aries energy vs TAC

// Energy vs detector number 2D histograms
TH2I  *ge_xtal, *bgo_xtal, *bgof_xtal, *bgos_xtal, *bgob_xtal, *bgoa_xtal, *labr_xtal, *labr_tac_xtal, *paces_xtal, *aries_xtal, *art_tac_xtal, *desw_e_xtal, *desw_tof_xtal;

// Time difference spectra
#define N_DT 26
char dt_handles[N_DT][32]={ "dt_ge_ge", "dt_ge_bgo", "dt_ge_sep", "dt_ge_zds",  // 0-3
                            "dt_ge_pac", "dt_ge_labr", "dt_ge_rcmp", "dt_pac_zds", // 4-7
                            "dt_pac_labr", "dt_rcmp_rcmp", "dt_ge_art", "dt_labr_art", // 8-11
                            "dt_paces_art", "dt_art_art", "dt_art_tac", "dt_zds_tac", "dt_labr_tac", "dt_labr_zds", // 12-17
                            "dt_dsw_dsw", "dt_dsw_ge", "dt_dsw_art", "dt_dsw_zds", // 18-21
                            "dt_zds_GRIF_CAEN_10ns", "dt_zds_GRIF_CAEN_2ns", "dt_dsw_dsw_2ns", "dt_dsw_zds_2ns"  }; // 22-25
TH1I  *dt_hist[N_DT];
TH1I  *dt_tacs_hist[N_LABR];

// Two-dimensional hitpatterns
TH2I *gg_hit, *bgobgo_hit, *aa_hit, *gea_hit, *lba_hit, *dsw_hit;

// En-En Coincidence matrices
TH2I *gg, *gg_ab, *gg_opp, *ge_paces, *ge_labr, *ge_rcmp, *labr_labr, *labr_zds, *labr_rcmp, *ge_art, *ge_zds, *paces_art, *labr_art, *art_art, *dsw_dsw, *ge_dsw, *art_dsw;

// Angular Correlation histograms
#define N_GE_ANG_CORR 52
TH2I  *gg_ang_corr_hist[N_GE_ANG_CORR];
#define N_GRG_ART_ANG_CORR 114
TH2I  *grg_art_ang_corr_hist[N_GRG_ART_ANG_CORR];
#define N_DSW_DSW_ANG_CORR 42
TH2I  *dsw_dsw_ang_corr_hist[N_DSW_DSW_ANG_CORR];
