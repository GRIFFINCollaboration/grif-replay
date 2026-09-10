// default_sort_griffin.h
// Header file for default_init_griffin.c,  default_presort_griffin.c, default_sort_griffin.c

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
#define E_2D_QED_SPECLEN        2048
#define E_3D_LBL_SPECLEN      160000  // 400*400
#define E_3D_TAC_SPECLEN         512
//#define T_SPEC_LENGTH     8192
//#define WV_SPEC_LENGTH    4096
#define DT_SPEC_LENGTH          1024
#define GE_ANGCOR_SPECLEN       4096
#define DSW_ANGCOR_SPECLEN      4096
#define QED_STRIP_THRESHOLD       30  // 30keV in all strips is required
#define QED_PIXEL_THRESHOLD       50  // 50keV used in building coincidences
#define NUM_QED_REORDERS          10
#define QED_COMPTON              0xF
#define QED_GAMMA_ENERGY         511  //  511keV
#define QED_GAMMA2_ENERGY       1274  // 1274keV
#define QED_GAMMA_ENERGY_WINDOW   20  //   20keV
#define QED_ANGLE_WINDOW          15  // 15 degrees

#define N_PPG_PATTERNS  7
#define MAX_ODB_PPG_CYCLES 50
#define MAX_CYCLES 512                          // Used as an axis length for some 2d histograms so must be a multiple of 16
#define N_HITPAT  7

#define N_PU_CLASSES 15         // HPGe pileup
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

#define N_DT 31   // Time difference
#define N_GE_ANG_CORR       52
#define N_GRG_ART_ANG_CORR 114
#define N_DSW_DSW_ANG_CORR  42
#define N_GE_COMP_POL       12 // 12 bins is 180 degrees divided into a bin width of 15 degrees

//#######################################################################
//######## Variables used across default init, presort and sort #########
//#######################################################################
extern int DEBUG_OUTPUT;
extern char subsys_handle[MAX_SUBSYS][8], subysy_name[MAX_SUBSYS][STRING_LEN];

extern int          odb_daqsize;// number of daq channels currently defined in the odb
extern int         subsys_table[MAX_DAQSIZE];
extern int        crystal_table[MAX_DAQSIZE]; // Ge/BGO have 4 "crystals" per clover
extern int        element_table[MAX_DAQSIZE]; // BGO have 5 elements per crystal
extern int       polarity_table[MAX_DAQSIZE]; // 1 is negative, 0 is positive, -1 is unset
extern short       address_chan[MAX_ADDRESS], *chan_address;
extern short  addr_table[MAX_DAQSIZE], *addrs;
extern char    chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
extern int   dtype_table[MAX_DAQSIZE], *dtypes;
extern float  gain_table[MAX_DAQSIZE], *gains, *offsets, *quads;
extern float  offs_table[MAX_DAQSIZE];
extern float  quad_table[MAX_DAQSIZE];
extern float  pileupk1[MAX_DAQSIZE][7];
extern float  pileupk2[MAX_DAQSIZE][7];
extern float  pileupE1[MAX_DAQSIZE][7];
extern float  crosstalk[MAX_DAQSIZE][3][16];
extern int subsys_initialized[MAX_SUBSYS];
extern int subsys_deadtime_count[MAX_SUBSYS];
extern int subsys_prg_ddtm[MAX_SUBSYS];
extern int previous_trig_acc[MAX_DAQSIZE];
extern Grif_event grif_event[PTR_BUFSIZE];

extern int presort_window_width;
extern int sort_window_width;
extern int ct_index[4][4];

typedef struct ppg_cycles_struct {
  char name[128];  int length; int codes[16]; int durations[16];
} ppg_cycles;

//#######################################################################
//######## Function declarations used across default init, presort and sort #########
//#######################################################################

extern float spread(int val);
extern int init_parameters_from_globals(Config *cfg);
extern int perform_pileup_correction(Grif_event *ptr, Grif_event *alt, int dt, int chan, int chan2, int i, int end_idx);
extern int init_chan_histos(Config *cfg);
extern int init_histos(Config *cfg, int subsystem);
extern int fill_chan_histos(Grif_event *ptr);
extern int fill_singles_histos(Grif_event *ptr);
extern int fill_coinc_histos(int win_idx, int frag_idx);

//#######################################################################
//########             PPG variables and patterns              ##########
//#######################################################################

// PPG patterns and handles
extern int ppg_patterns[N_PPG_PATTERNS];
extern char ppg_handles[N_PPG_PATTERNS][32];
extern char ppg_names[N_PPG_PATTERNS][32];

// Definitions for the Current cycle of this run
// These are derived from the ODB settings at BOR
extern long ppg_cycle_duration;             // Length of one cycle in timestamp units
extern char ppg_cycle_name[128];            // Name of this current cycle
extern long ppg_cycle_length;               // Number of patterns/durations for this current cycle
extern long ppg_cycle_pattern_duration[16]; // Length of each pattern in timestamp units
extern int  ppg_cycle_pattern_code[16];     // Index of each pattern for use with the ppg_patterns array
extern int  ppg_cycles_active;              // Cycles active or made inactive if set to Source/constant beam-on etc.

// These variables are updated in pre_sort_enter at each PPG pattern change
extern long ppg_last_ptr_ts;     // Previous event timestamp. Avoids rare bug where single events are out of order.
extern int ppg_current_pattern;  // Index of the current PPG cycle pattern for use with the ppg_patterns array
extern int ppg_cycle_number;     // Current cycle number. Cycles counted from zero at beginning of run
extern long ppg_cycle_start;     // Timestamp of the start of the current cycle
extern long ppg_cycle_end;       // Timestamp of the end of the current cycle
extern int ppg_cycle_step;       // Current pattern number within this cycle. Patterns counted from zero at beginning of cycle
extern long ppg_pattern_start;   // Timestamp of the start of the current pattern
extern long ppg_pattern_end;     // Timestamp of the end of the current pattern
extern long ppg_bin_end;         // Timestamp of the end of the current bin (used for deadtime histogram)

// Spectra for cycles
// The binning factor and gamma-energy gates will ultimately be set from a Global at BOR
extern long ppg_cycles_binning_factor;
extern int ppg_cycles_gamma_gate_min, ppg_cycles_gamma_gate_max;
extern TH1I   *ge_cycle_activity, *zds_cycle_activity; // Activity over cycle time, sum of all cycles
extern TH2I   *ge_e_vs_cycle_time;                     // Energy vs time within the cycle
extern TH1I   *ge_cycle_code[N_PPG_PATTERNS];          // Energy spectrum for each PPG pattern
extern TH2I   *gg_cycle_code[N_PPG_PATTERNS];          // Ge-Ge 2D histogram for each PPG pattern
extern TH1I   *gea_cycle_num[MAX_CYCLES];               // Activity over cycle time for each indivdual cycle, GRGA
extern TH1I   *gea_cycle_num_sh[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGA, single_hit only
extern TH1I   *gea_cycle_num_pu[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGA, pileup only
extern TH1I   *gea_cycle_num_dt[MAX_CYCLES];            // Deadtime over cycle time for each indivdual cycle, GRGA
extern TH1I   *gea_cycle_num_g[MAX_CYCLES];             // Activity over cycle time for each indivdual cycle, GRGA, gamma-gated
extern TH1I   *gea_cycle_num_sh_g[MAX_CYCLES];          // Activity over cycle time for each indivdual cycle, GRGA, gamma-gated, single_hit only
extern TH1I   *geb_cycle_num[MAX_CYCLES];               // Activity over cycle time for each indivdual cycle, GRGB
extern TH1I   *geb_cycle_num_sh[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGB, single_hit only
extern TH1I   *geb_cycle_num_pu[MAX_CYCLES];            // Activity over cycle time for each indivdual cycle, GRGB, pileup only
extern TH1I   *geb_cycle_num_dt[MAX_CYCLES];            // Deadtime over cycle time for each indivdual cycle, GRGB
extern TH1I   *geb_cycle_num_g[MAX_CYCLES];             // Activity over cycle time for each indivdual cycle, GRGB, gamma-gated
extern TH1I   *geb_cycle_num_sh_g[MAX_CYCLES];          // Activity over cycle time for each indivdual cycle, GRGB, gamma-gated, single_hit only
extern TH2I   *cycle_num_vs_ge;                        // 2D histogram of cycle # vs the cycle time.
extern TH2I   *cycle_num_vs_ge_sh_g;                   // 2D histogram of cycle # vs the cycle time for 1809 gated/NP.
extern TH2I   *cycle_num_vs_ge_dt;                     // 2D histogram of cycle # vs the cycle time for deadtime .
extern TH2I   *cycle_num_vs_sh;                        // 2D histogram of cycle # vs the cycle time for NP.
extern TH2I   *cycle_num_vs_pu;                        // 2D histogram of cycle # vs the cycle time for PU.
extern TH2I   *cycle_num_vs_ge_b;                      // 2D histogram of cycle # vs the cycle time. GRGB
extern TH2I   *cycle_num_vs_ge_b_sh_g;                 // 2D histogram of cycle # vs the cycle time for 1809 gated/NP. GRGB
extern TH2I   *cycle_num_vs_ge_b_dt;                   // 2D histogram of cycle # vs the cycle time for deadtime . GRGB
extern TH2I   *cycle_num_vs_sh_b;                      // 2D histogram of cycle # vs the cycle time for NP. GRGB
extern TH2I   *cycle_num_vs_pu_b;                      // 2D histogram of cycle # vs the cycle time for PU. GRGB
extern TH2I   *cycle_num_vs_geEnergy[N_HPGE];          // 2D histogram of cycle # vs the Ge energy spectrum for that cycle (Use for monitoring gain drifts).
extern TH2I   *cycle_num_vs_qedEnergy[N_QED_POS];      // 2D histogram of cycle # vs the Ge energy spectrum for that cycle (Use for monitoring gain drifts).

//#######################################################################
//########        Individual channel singles HISTOGRAMS        ##########
//#######################################################################

// Pulse height and energy
extern TH1I   *ts_hist; // timestamp
extern TH1I   *gc_hist; // GRIF-CAEN hitpattern for checking coincidences
extern TH1I   *ph_hist[MAX_DAQSIZE];
extern TH1I   *e_hist[MAX_DAQSIZE];
//TH1I *wave_hist[MAX_DAQSIZE];

extern TH1I  *hit_hist[N_HITPAT], *mult_hist[MAX_SUBSYS];

// DESCANT Wall
extern TH1I  *desw_tof[N_DES_WALL];                // Time-Of-Flight
extern TH1I  *desw_tof_corr[N_DES_WALL];           // corrected Time-Of-Flight
extern TH1I  *desw_tof_psd[N_DES_WALL];            // corrected Time-Of-Flight, PSD gated
extern TH1I  *desw_psd[N_DES_WALL];                // Pulse Shape Discrimination

//#######################################################################
//########                PRESORT Time Gates                   ##########
//#######################################################################

// The definition of the time difference gate in 10 nanosecond units.
// The value is the maximum time difference in 10 nanosecond units.
// The default values set here are replaced by the Global value at start of sorting.
extern int bgo_window_min, addback_window_min, rcmp_fb_window_min, qed_fb_window_min, lbl_tac_window_min;
extern int art_tac_window_min, zds_tac_window_min, desw_beta_window_min, bgo_window_max;
extern int addback_window_max, rcmp_fb_window_max, qed_fb_window_max, lbl_tac_window_max;
extern int art_tac_window_max, zds_tac_window_max, desw_beta_window_max;

//#######################################################################
//########                Coincidence Time Gates               ##########
//#######################################################################

// The definition of the time difference gate in 10 nanosecond units.
// First and second index are the subsystem index numbers
// The value is the maximum time difference in 10 nanosecond units.
// Default is 250 nanoseconds, replaced by the Global value at start of sorting
extern int time_diff_gate_min[MAX_SUBSYS][MAX_SUBSYS];
extern int time_diff_gate_max[MAX_SUBSYS][MAX_SUBSYS];

//#######################################################################
//########          Sums and coincidence  HISTOGRAMS           ##########
//#######################################################################

extern TH1I  *ge_pu_type; // The value of pile-up type
extern TH1I  *ge_nhits_type; // The value of nhits type
extern TH1I  *ge_pu_class; // The value of pile-up class
extern TH1I  *ge_sum_class[N_PU_CLASSES]; // Ge energy for each pile-up value
extern TH2I  *ge_e_vs_k_class[N_PU_CLASSES]; // Ge energy vs k for each pile-up value
extern TH2I  *ge_xtal_1hit, *ge_xtal_2hit, *ge_xtal_3hit; // Ge energy vs crystal number for 1, 2, 3 hit PU.
extern TH1I  *ge_1hit[N_HPGE]; // Ge single hit events
extern TH1I  *ge_2hit[N_HPGE]; // Ge 2-hit pileup events
extern TH1I  *ge_3hit[N_HPGE]; // Ge 3-hit pileup events
extern TH1I  *ge_pu_dt12; // Time difference between first and second Hit
extern TH1I  *ge_pu_dt13; // Time difference between first and third Hit
extern TH2I  *ge_e_vs_k_2hit_first[N_HPGE]; // Ge Hit1 energy vs k1 for 2-hit pileup events
extern TH2I  *ge_e_vs_k_2hit_second[N_HPGE]; // Ge Hit2 energy vs k2 for 2-hit pileup events
extern TH2I  *ge_PU2_e2_v_k_gatedxrays[N_HPGE]; // Ge e2 vs k2 for fixed e1 energy gates
extern TH2I  *ge_PU2_e2_v_k_gated1408[N_HPGE]; // Ge e2 vs k2 for fixed e1 energy gates

// for most subsystem-pair-combinations, there is a
// a 1d time-difference and a 2d ecal-vs-ecal matrix
extern TH2I *subsys_e_vs_e[MAX_SUBSYS][MAX_SUBSYS];
extern TH1I *subsys_dt[MAX_SUBSYS][MAX_SUBSYS], *subsys_dcfd[MAX_SUBSYS][MAX_SUBSYS];
extern TH1I *tac_lbl_ts_diff[N_TACS];

// HPGe (ge_sum is sum of crystal energies, ge_sum_b is beta-gated)
extern TH1I  *ge_ab_e[N_CLOVER], *ge_ab_sup_e[N_CLOVER], *ge_sum_ab, *ge_sum_ab_sup, *ge_sum_ab_sup_rej;
extern TH1I  *ge_sum, *ge_sum_us, *ge_sum_ds, *ge_sum_hem[2], *ge_sum_ab_us, *ge_sum_ab_ds;
extern TH1I  *ge_sum_b, *ge_sum_b_ab, *ge_sum_b_sep, *ge_sum_b_sep_brems, *ge_sum_b_ab_sep_brems, *ge_sum_b_zds;
extern TH1I  *ge_sum_b_art, *ge_sum_b_art_brems, *ge_sum_b_artT, *ge_sum_b_artR, *ge_sum_b_artS;

// ARIES, PACES and LaBr3
extern TH1I  *aries_sum, *paces_sum, *paces_sum_b, *labr_sum;

// RCMP
extern TH1I  *rcmp_sum, *rcmp_fb_sum;
extern TH2I  *rcmp_strips[N_RCMP_POS], *rcmp_hit[N_RCMP_POS], *rcmp_fb[N_RCMP_POS], *rcmp_x_ge_hit, *rcmp_y_ge_hit;

// QED
extern TH1I  *qed_sum, *qed_fb_sum;  // qed_sum is sum of strip energies, fb is with front-back coincidence
extern TH2I  *qed_strips[N_QED_POS], *qed_hit[N_QED_POS], *qed_fb[N_QED_POS], *qed_p_ge_hit, *qed_n_ge_hit;
extern TH2I  *qedp_ge_theta[N_QED_POS*N_QED_STRIPS], *qedn_ge_theta[N_QED_POS*N_QED_STRIPS];
extern TH2I  *qedE_ge_theta_sum, *qed_geE_theta_sum, *qed_E_totE_sum_t, *qed_geE_totE_sum_t, *qedE_ge_theta_sum_t, *qed_geE_theta_sum_t, *qedE_ge_thetaI_sum_t, *qed_geE_thetaDiff_sum_t, *qed_geE_thetaI_sum_t;
extern TH2I  *qedE_ge_theta_sum_c, *qed_geE_theta_sum_c, *qedE_ge_theta_sum_c_g, *qed_geE_theta_sum_c_g, *qedE_ge_theta_sum_c_s, *qed_geE_theta_sum_c_s, *ge_qed_c;
extern TH2I  *qed_geE_theta_clov[N_CLOVER], *qed_geE_theta_clov_t[N_CLOVER], *qed_E_theta_dssd[N_QED_POS], *qed_geE_theta_dssd[N_QED_POS];
extern TH2I  *qed_angle_test_g, *qed_angle_test_s, *qedE_ge_dt, *qed_geE_dt, *qedE_ge_dt_c, *qed_geE_dt_c, *qed_theta_dt_c, *qed_theta_dt, *qed_theta_dt_cfd, *qed_dcs_omega_dt, *qed_dcs_omega_dtx, *qedx_dcs_omega_dt[N_QED_POS], *qed_theta1_vs_theta2, *qed_theta1_azi, *qed_theta2_azi, *qed2_theta1_vs_theta2, *qed2_theta1_azi, *qed2_theta2_azi;
extern TH1I  *qed_dcs_omega, *qed_dcs_omega_t, *qed_dcs_azi, *qed_dcs_azi_t, *qed_dcs_azi_tg, *qed_dcs_azi_TRWF, *qed_dcs_azi_TRWF_t, *qed_dcs_azi_TRWF_tg, *qed_delta_theta1_theta2, *qed_sum_theta1_theta2;
extern TH1I  *qed_wf_dcs_azi, *qed_wf_omega;
extern TH2I  *qed_angle_theta_g, *qed_angle_theta_s, *qed_angle_phi_s, *qed_angle_phi_a, *qed_angle_phi_b;
extern TH1I  *qed_theta, *qed_phi_s, *qed_phi_a, *qed_phi_b;
extern TH1I  *qed_dcs_azi_bins1, *qed_dcs_azi_bins2, *qed_dcs_azi_bins3, *qed_dcs_azi_bins4, *qed_dcs_azi_bins5, *qed_dcs_azi_bins6, *qed_dcs_azi_bins7, *qed_dcs_azi_bins8, *qed_dcs_azi_bins8a, *qed_dcs_azi_bins9, *qed_dcs_azi_bins10;
extern TH1I  *qed_dcs_azi_TRWF_bins1, *qed_dcs_azi_TRWF_bins2, *qed_dcs_azi_TRWF_bins3, *qed_dcs_azi_TRWF_bins4, *qed_dcs_azi_TRWF_bins5, *qed_dcs_azi_TRWF_bins6, *qed_dcs_azi_TRWF_bins7, *qed_dcs_azi_TRWF_bins8, *qed_dcs_azi_TRWF_bins8a, *qed_dcs_azi_TRWF_bins9, *qed_dcs_azi_TRWF_bins10;
extern TH1I  *qed_phi_bins1, *qed_phi_bins2, *qed_phi_bins3, *qed_phi_bins4, *qed_phi_bins5, *qed_phi_bins6, *qed_phi_bins7, *qed_phi_bins8, *qed_phi_bins9, *qed_phi_bins10;
extern TH2I  *dcsaE_ge_theta, *dcsa_geE_theta, *dcsa_theta_azi, *qed_dcs_omega_dt_TRWF, *dcsa_theta_azi_ge;
extern TH1I  *dcsa_cs_omega, *dcsa_cs_omega_ge, *dcsa_theta;
extern TH2I  *dcsbE_ge_theta, *dcsb_geE_theta, *dcsb_theta_azi, *dcsb_theta_azi_ge;
extern TH1I  *dcsb_cs_omega, *dcsb_cs_omega_ge, *dcsb_theta;
extern TH2I  *qed_qed_23, *qed_qed_23_theta2, *qed_qed_23_theta3, *qed_qed_23_totv2, *qed_qed_23_totv3, *qed_qed_12, *qed_qed_12_theta1, *qed_qed_12_theta2, *qed_qed_12_totv1, *qed_qed_12_totv2, *qed_qed_14, *qed_qed_14_theta1, *qed_qed_14_theta4, *qed_qed_14_totv1, *qed_qed_14_totv4;
extern TH1I  *qed_qed_23dt, *qed_qed_12dt, *qed_qed_14dt;
extern TH2I  *qed_ge_weight;

// DESCANT WALL
extern TH1I  *desw_sum_e, *desw_sum_tof, *desw_sum_psd;  // Sums of energies and corTOF and PSD
extern TH1I  *desw_sum_e_b, *desw_sum_tof_b;       // Beta-tagged Sums of energies and corTOF
extern TH1I  *desw_sum_e_nn, *desw_sum_tof_nn;     // fold>2 Sums of energies and corTOF
extern TH1I  *desw_sum_e_nn_a, *desw_sum_tof_nn_a; // fold>2, angle>60 Sums of energies and corTOF
extern TH2I  *desw_psd_e, *desw_psd_tof;           // PSD vs energies or corrected-TOF
extern TH2I  *desw_psd_q,*desw_psd_cc,*desw_q_cc,*desw_q_tof,*desw_cc_tof,*desw_psd_zdse; //

// TAC spectra
extern TH1I *tac_labr_hist[(int)((N_LABR)*(N_LABR-1)/2)+2]; // this index numbers are the LaBr-LaBr position numbers
extern TH1I *tac_labr_hist_uncal[(int)((N_LABR)*(N_LABR-1)/2)+2]; // this index numbers are the LaBr-LaBr position numbers
// One additional histogram (2_1) needed for Compton Walk corrections
extern TH2I *tac_labr_CompWalk[N_LABR];         // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
extern TH2I *tac_labr_CompWalk0;                // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
extern int tac_labr_hist_index[N_LABR][N_LABR]; // index for filling tac_labr_hist from LBL id numbers
extern TH1I *tac_gated_lbl[N_LABR];             // TAC-gated LBL energy spectrum to check anode threshold in analogue CFD
extern TH1I *final_tac[N_TACS], *final_tac_sum; // Final TAC spectra after all calibration and alignment
extern TH1I *tac_aries_lbl[N_LABR];             // this index number is the LaBr position number
extern TH1I *tac_aries_art[N_ARIES];            // this index number is the Aries position number
extern TH1I *tac_aries_lbl_sum;                 // ARIES TAC sum spectrum of all LBLs
extern TH1I *tac_aries_art_sum;                 // ARIES TAC sum spectrum of all ARTs
extern TH1I *aries_tac;                         // aries_tac gated on 1275keV peak
extern TH1I *aries_tac_Egate;                   // aries_tac gated on 1275keV peak
extern TH1I *aries_tac_artEn;                   // aries energy in coincidence with TAC
extern TH2I *lblE_tac, *zdsE_tac, *ariesE_tac;  // lbl or zds or aries energy vs TAC
extern TH2I *lbl_lbl_tac;                       // A special 3d histogram disguised as a 2d histogram

// 2D Energy vs detector number
extern TH2I *ge_xtal, *geb_xtal, *bgo_xtal, *bgof_xtal, *bgos_xtal, *bgob_xtal, *bgoa_xtal, *labr_xtal;
extern TH2I *labr_tac_xtal, *paces_xtal, *sceptar_xtal, *aries_xtal, *art_tac_xtal, *desw_e_xtal, *desw_tof_xtal;

extern TH1I  *dt_hist[N_DT], *dcfd_hist[N_DT], *dt_tacs_hist[N_TACS];

// 2D hitpatterns
extern TH2I *gg_hit, *bgobgo_hit, *aa_hit, *gea_hit, *lba_hit, *dsw_hit;

// 2D Energy vs Energy Coincidence matrices
extern TH2I *gea_self_dt,*geb_self_dt;
extern TH2I *gg, *gg_ab, *gg_opp, *gg_ab_opp, *ge_bgo, *ge_paces, *ge_labr, *ge_rcmp, *labr_labr, *labr_zds, *labr_rcmp;
extern TH2I *ge_art, *ge_zds, *paces_art, *labr_art, *art_art, *dsw_dsw, *ge_dsw, *art_dsw, *ge_qed, *qed_qed, *comp_comp, *ge_comp, *geadd_comp, *ge_dcs, *geadd_dcs, *comp_dcs;
extern TH1I *gg_energy[N_HPGE];

// Angular Correlation histograms
extern TH2I  *gg_angcor_110[N_GE_ANG_CORR], *gg_angcor_145[N_GE_ANG_CORR], *ge_art_angcor[N_GRG_ART_ANG_CORR], *dsw_angcor[N_DSW_DSW_ANG_CORR];

// Compton Polarimetry histograms
extern TH2I  *comp_pol_angles_110, *comp_pol_angles_145, *gg_comp_pol_110[N_GE_COMP_POL], *gg_comp_pol_145[N_GE_COMP_POL];

// Isomer Spectroscopy
extern TH1I  *ge_isomer_popu, *ge_isomer_depop;
extern TH2I  *gg_dt, *gb_dt;

// Crosstalk Analysis
extern TH2I  *ct_e_vs_dt_B[N_HPGE], *ct_e_vs_dt_G[N_HPGE], *ct_e_vs_dt_R[N_HPGE], *ct_e_vs_dt_W[N_HPGE];

////////////////////////////////////
////////////////////////////////////

extern int tac_ts_offset[12]; // LBT (TAC) timestamp offset values.
extern int tac_lbl_combo_offset[(int)((N_LABR)*(N_LABR-1)/2)+2]; // TAC coincidence combination offsets.

// BGO HV alignment histograms
extern TH1I *ge_bgo_gated[N_BGO];
