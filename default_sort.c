#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "config.h"
#include "grif-format.h"
#include "histogram.h"
#include "grif-angles.h"

#define SUBSYS_HPGE     0
#define SUBSYS_BGO      1
#define SUBSYS_SCEPTAR  2
#define SUBSYS_PACES    3
#define SUBSYS_LABR_BGO 4
#define SUBSYS_LABR_T   5
#define SUBSYS_LABR_L   6
#define SUBSYS_DESCANT  7
#define SUBSYS_ARIES    8
#define SUBSYS_ZDS      9
#define SUBSYS_RCMP    10

#define MAX_SUBSYS 24
static char subsys_handle[MAX_SUBSYS][8] = {
  "GRG", "GRS", "SEP",  "PAC",
  "LBS", "LBT", "LBL",  "DSC",
  "ART", "ZDS", "RCS",  "XXX",
  "",    "",    "",    "",
  "",    "",    "",    "",
  "",    "",    "",    ""
};
static char subsys_name[MAX_SUBSYS][STRING_LEN] = {
   "Griffin",  "BGO",   "SCEPTAR",   "PACES", //  0- 3
   "LaBrS",    "LaBrT", "LaBrX",   "Descant", //  4- 7
   "ARIES",    "ZDS",   "RCMP",        "XXX", //  8-11
   "",         "",      "",         "",       // 12-15
   "",         "",      "",         "",
   "",         "",       "",        "Unknown"
}; // final entry will be used if not found - make sure it is not empty


int          odb_daqsize;// number of daq channels currently defined in the odb
int     subsys_dtype_mat[MAX_SUBSYS][MAX_SUBSYS]; // map using names and dtypes
int         subsys_dtype[MAX_SUBSYS];             // subsys->dtype
int         dtype_subsys[MAX_SUBSYS];             // dtype->subsys
int        crystal_table[MAX_DAQSIZE]; // Ge/BGO have 4 "crystals" per clover
int        element_table[MAX_DAQSIZE]; // BGO have 5 elements per crystal
int       polarity_table[MAX_DAQSIZE]; // 1 is negative, 0 is positive, -1 is unset
int         output_table[MAX_DAQSIZE]; // 1 is A, 0 is B, -1 is X or unknown
short       address_chan[MAX_ADDRESS];
static short  addr_table[MAX_DAQSIZE]; short   *addrs = addr_table;
       char    chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
static int   dtype_table[MAX_DAQSIZE]; int    *dtypes = dtype_table;
static float  gain_table[MAX_DAQSIZE]; float   *gains = gain_table;
static float  offs_table[MAX_DAQSIZE]; float *offsets = offs_table;
static float  quad_table[MAX_DAQSIZE]; float   *quads = quad_table;
static short *chan_address = addr_table;
extern Grif_event grif_event[MAX_COINC_EVENTS];

// Default sort function declarations
extern int init_chan_histos(Config *cfg);
extern int init_singles_histos(Config *cfg);
extern int init_coinc_histos(Config *cfg);
extern int fill_chan_histos(Grif_event *ptr);
extern int fill_singles_histos(Grif_event *ptr);
extern int fill_coinc_histos(int win_idx, int frag_idx);

// odb tables need to be transferred into config, which is saved with histos
int init_default_histos(Config *cfg, Sort_status *arg)
{
   Cal_coeff *cal;
   int i, j;
   cfg->odb_daqsize = odb_daqsize;
   for(i=0; i<odb_daqsize; i++){ // update config with odb info
      edit_calibration(cfg, chan_name[i], offsets[i], gains[i], quads[i],
                       chan_address[i],  dtype_table[i], arg->cal_overwrite );
   }
   // ALSO need to transfer config info to the arrays that are used in sort
   if( arg->cal_overwrite == 0 ){ // overwrite = 0 => USE CONFIG NOT ODB
      for(i=0; i<odb_daqsize; i++){
         cal = cfg->calib[i];
         if( strcmp(chan_name[i], cal->name) != 0 ){ // conf not in odb order
            for(j=0; j<cfg->ncal; j++){ cal = cfg->calib[j];
               if( strcmp(chan_name[i], cal->name) == 0 ){ break; }
            }
            if( j == cfg->ncal ){ continue; } // not found in config
         }
         offsets[i]=cal->offset; gains[i]=cal->gain;  quads[i]=cal->quad;
      }
   }
   init_chan_histos(cfg);
   init_singles_histos(cfg);
   init_coinc_histos(cfg);

   return(0);
}

//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#######################################################################

#define NUM_CLOVER 16

float spread(int val){ return( val + rand()/(1.0*RAND_MAX) ); }
int GetIDfromAddress(int addr) { return( address_chan[addr] ); }

int apply_gains(Grif_event *ptr)
{
  int tac_ts_offset[8] = {48,60,100,49,119,59,133,14};
   float energy;
   int chan;
   if( (chan=ptr->chan) >= odb_daqsize ){
      fprintf(stderr,"unpack_event: ignored event in chan:%d [0x%04x]\n",
            	                                  chan, ptr->address );
      return(-1);
   }

   ptr->energy = energy = ( ptr->integ == 0 ) ? 0 : spread(ptr->q)/ptr->integ;
   ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);

   // Assign the Sub System index based on dtype
   // The dtype to subsys mapping was determined from the PSC table in the function gen_derived_odb_tables()
   if( ptr->dtype >= 0 && ptr->dtype < MAX_SUBSYS ){
      ptr->subsys = dtype_subsys[ptr->dtype];
   } else { ptr->subsys = MAX_SUBSYS-1; }


   // The TAC module produces its output signal around 2 microseconds later
   // than the start and stop detector signals are processed.
   if( ptr->subsys == SUBSYS_LABR_T){
     ptr->ts -= 150 + tac_ts_offset[crystal_table[ptr->chan]-1]; // Subtract 2 microseconds from TAC timestamps plus a more precise offset
   }

   // NOBODY CURRENTLY USES e2,e3,e4 ...
   //if( ptr->integ2 != 0 ){
   //   energy = ptr->energy2 = spread(ptr->q2)/ptr->integ2;
   //   ptr->e2cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
   //}
   //if( ptr->integ3 != 0 ){
   //   energy = ptr->energy3 = spread(ptr->q3)/ptr->integ3;
   //   ptr->e3cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
   //}
   //if( ptr->integ4 != 0 ){
   //   energy = ptr->energy4 = spread(ptr->q4)/ptr->integ4;
   //   ptr->e4cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
   //}
   return(0);
}

int default_sort(int win_idx, int frag_idx, int flag)
{
   Grif_event *ptr;
   int i;

   // sort first event, even if window only contains that event
   // (usually require at least two)
   for(i=win_idx; ; i++){ ptr = &grif_event[i];
      if( i >= MAX_COINC_EVENTS ){ i=0; } // WRAP
      if( i != win_idx && flag == SORT_ONE ){ break; }
      if( i != win_idx && i==frag_idx ){ break; }
      if( ptr->dtype == 15 ){ if( i==frag_idx ){ break; } continue; } // scalar
      if( ptr->chan == -1 ){
         printf("DefSort: UNKNOWN_CHAN type=%d\n", ptr->dtype);
         if( i==frag_idx ){ break; } continue;
      } //  ????
      fill_chan_histos(ptr);
      fill_singles_histos(ptr);
      if( i==frag_idx ){ break; }
   }
   fill_coinc_histos(win_idx, frag_idx);
   return(0);
}

// Presort - do Suppression and Addback here
//  - frag_idx is about to leave coinc window (which ends at end_idx)
//    check other frags in window for possible suppression and/or summing
//  also calculate multiplicities[store in frag_idx only]
int pre_sort(int frag_idx, int end_idx)
{
  Grif_event *alt, *ptr = &grif_event[frag_idx];
  int bgo_window = 20, addback_window = 20;
  int rcmp_fb_window = 10;
  int art_tac_window = 25;
  int i, dt;

  //printf("\n");

  //if( ptr->dtype !=  0 ){ return(0); } // not Ge event
  //if( ptr->dtype == 15 ){ return(0); } // scalar
  i = frag_idx; ptr->fold = 1;
  while( i != end_idx ){ // need at least two events in window
    if( ++i >=  MAX_COINC_EVENTS ){ i=0; } alt = &grif_event[i]; // WRAP

    // Determine fold
    if( alt->subsys == ptr->subsys ){ ++ptr->fold; }

    // Determine absolute time difference between timestamps
    dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }

    // SubSystem-specific pre-processing
    switch(ptr->subsys){
      case SUBSYS_HPGE:
      if( dt < bgo_window && alt->subsys == SUBSYS_BGO && !ptr->suppress ){
        // could alternatively use crystal numbers rather than clover#
        //    (don't currently have this for BGO)
        if( crystal_table[ptr->chan]/16 == crystal_table[alt->chan]/16 ){ ptr->suppress = 1; }
      }
      // Germanium addback -
      //    earliest fragment has the sum energy, others are marked -1
      // Ensure GRG events are both A or both B type using output_table
      // Remember the other crystal channel number in ab_alt_chan for use in Compton Polarimetry
      if( dt < addback_window && alt->subsys == SUBSYS_HPGE && (output_table[ptr->chan] == output_table[alt->chan])){
        if( alt->esum >= 0 && crystal_table[alt->chan]/16 == crystal_table[ptr->chan]/16 ){
          ptr->esum += alt->esum; alt->esum = -1; ptr->ab_alt_chan = alt->chan;
        }
      }
      break;
      case SUBSYS_RCMP:
      // RCMP Front-Back coincidence
      // Ensure its the same DSSD and that the two events are front and back
      // The charged particles enter the P side and this has superior energy resolution
      // Ensure the energy collected in the front and back is similar
      ptr->esum = -1; // Need to exclude any noise and random coincidences.
      if( dt < rcmp_fb_window && alt->subsys == SUBSYS_RCMP && (ptr->ecal>0 && ptr->ecal<32768)){
        if((crystal_table[ptr->chan] == crystal_table[alt->chan]) && (polarity_table[ptr->chan] != polarity_table[alt->chan]) && (alt->ecal > 0  && ptr->ecal<32768)){
          if( ((ptr->ecal / alt->ecal)<=1.1 && (ptr->ecal / alt->ecal)>=0.9)){
            // Ensure esum comes from P side, but use this timestamp
            ptr->esum = polarity_table[ptr->chan]==0 ? ptr->ecal : (polarity_table[ptr->chan]==1 ? alt->ecal : -1);
          }
        }
      }
      break;
      case SUBSYS_LABR_T:
      // TAC spectra
      // For TAC08 the start is ARIES and the stop is any of the LaBr3. So this is three detector types.
      // Here in the presort we will remember the ARIES tile that is in coincidence with the TAC.
      // In the TAC event we save the tile chan as ab_alt_chan, and the tile energy as e2cal.
      // So later in the main coincidence loop we only need to examine LBL and TAC.
      if(dt < art_tac_window && alt->subsys == SUBSYS_ARIES && alt->ecal > 5){
      ptr->ab_alt_chan = alt->chan; ptr->e2cal = alt->ecal;
      }
      break;
      default: // Unrecognized or unprocessed subsys type
      break;
    }// end of switch

  }// end of while
  return(0);
}

//#######################################################################
//########        Individual channel singles HISTOGRAMS        ##########
//#######################################################################

#define MULT_SPEC_LENGTH   128
#define E_SPEC_LENGTH     8192
#define E_TAC_SPEC_LENGTH 16384
#define E_2D_SPEC_LENGTH  4096
#define E_2D_RCMP_SPEC_LENGTH  6400
//#define T_SPEC_LENGTH     8192
//#define WV_SPEC_LENGTH    4096

#define N_HITPAT  7
char hit_handles[N_HITPAT][32]={ "q_hit","e_hit","t_hit","w_hit","r_hit", "s_hit", "d_hit" };
char   hit_names[N_HITPAT][32]={
   "Pulse_Height", "Energy", "Time", "Waveform", "Rate", "Subsys", "DetType",
};
TH1I *ts_hist;

TH1I  *hit_hist[N_HITPAT], *mult_hist[MAX_SUBSYS];
TH1I   *ph_hist[MAX_DAQSIZE];
TH1I    *e_hist[MAX_DAQSIZE];
//TH1I *wave_hist[MAX_DAQSIZE];

int init_chan_histos(Config *cfg)
{                      // 1d histograms for Q,E,T,Wf for each channel in odb
   char title[STRING_LEN], handle[STRING_LEN];
   int i, j;

   open_folder(cfg, "Hits_and_Sums");
   open_folder(cfg, "Hits");
   for(i=0; i<N_HITPAT; i++){ // Create Hitpattern spectra
      sprintf(title,  "Hitpattern_%s",    hit_names[i] );
      hit_hist[i] = H1_BOOK(cfg, hit_handles[i], title, MAX_DAQSIZE, 0, MAX_DAQSIZE);
   }
   ts_hist = H1_BOOK(cfg, "ts", "Timestamp", 163840, 0, 163840);
   close_folder(cfg);
   open_folder(cfg, "Multiplicities");
   for(i=0; i<MAX_SUBSYS; i++){ mult_hist[i] = NULL;
      if( strcmp(subsys_handle[i],"XXX") == 0 ||
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
         if( strcmp(subsys_handle[j],"XXX") == 0 ||
                   strlen(subsys_handle[j]) == 0 ){ continue; }
         if( memcmp(subsys_handle[j], chan_name[i], 3) == 0 ){ break; }
      }                              // stop at final entry "unknown"
      open_folder(cfg, subsys_name[j]);
      open_folder(cfg, "Energy");
      sprintf(title,  "%s_Energy",         chan_name[i] );
      sprintf(handle, "%s_E",              chan_name[i] );
      if( strcmp(subsys_handle[j],"LBT") == 0){
        e_hist[i] = H1_BOOK(cfg, handle, title, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
      }else{
        e_hist[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
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
        ph_hist[i] = H1_BOOK(cfg, handle, title, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
      }else{
        ph_hist[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
      }
      close_folder(cfg);
      close_folder(cfg);
      if( strcmp(subsys_handle[j],"XXX") == 0 ){     // suppress XXX channels
         ph_hist[i]->suppress = e_hist[i]->suppress =
            /*  cfd_hist[i]->suppress = wave_hist[i]->suppress = */ 1;
      }
   }
   return(0);
}

int fill_chan_histos(Grif_event *ptr)
{
   static int event;
   int chan, sys;

   if( (chan = ptr->chan) == -1 ){
      return(-1);
   }
   if( ++event < 16384 ){
      ts_hist -> Fill(ts_hist, event,  (int)(ptr->timestamp/100));
   } else if( (event % 1000) == 0 ){
      ts_hist -> Fill(ts_hist, 16367+(int)(event/1000),  (int)(ptr->timestamp/100));
   }
   ph_hist[chan] -> Fill(ph_hist[chan],  (int)ptr->energy,     1);
   e_hist[chan]  -> Fill(e_hist[chan],   (int)ptr->ecal,       1);

   hit_hist[0]   -> Fill(hit_hist[0],    chan,            1);
   if( ptr->ecal        >= 1 ){ hit_hist[1] -> Fill(hit_hist[1], chan, 1);
                          hit_hist[6] -> Fill(hit_hist[6], ptr->dtype, 1);
                            sys = ptr->subsys;
                            if( sys >=0 && sys < MAX_SUBSYS ){
                              hit_hist[5] -> Fill(hit_hist[5], sys, 1);
                            }
                         }
   if( ptr->cfd         != 0 ){ hit_hist[2] -> Fill(hit_hist[2], chan, 1); }
   if( ptr->wf_present  != 0 ){ hit_hist[3] -> Fill(hit_hist[3], chan, 1); }
   if( ptr->scl_present != 0 ){ hit_hist[4] -> Fill(hit_hist[4], chan, 1); }
   if( ptr->scl_present != 0 ){ hit_hist[4] -> Fill(hit_hist[4], chan, 1); }  // Why is scaler incremented twice? Copy/paste error?

   return(0);
}

//#######################################################################
//########               Sums and coinc  HISTOGRAMS            ##########
//#######################################################################

TH1I  *ge_ab_e[NUM_CLOVER];
TH1I  *ge_sum, *ge_sum_ab;  // ge_sum is sum of crystal energies
TH1I  *ge_sum_us, *ge_sum_ds, *ge_sum_ab_us, *ge_sum_ab_ds;  // ge_sum is sum of crystal energies
TH1I  *ge_sum_b; // beta-gated gamma sum spectrum
TH1I  *paces_sum;  // paces_sum is sum of crystal energies
TH1I  *labr_sum;  // labr_sum is sum of crystal energies
TH1I  *aries_sum;  // aries_sum is sum of tile energies
TH1I  *aries_tac;  // aries_tac gated on 1475keV peak
TH1I  *rcmp_sum, *rcmp_fb_sum;  // rcmp_sum is sum of strip energies, fb is with front-back coincidence
TH2I  *ge_xtal, *bgo_xtal, *bgof_xtal, *bgos_xtal, *bgob_xtal, *bgoa_xtal, *labr_xtal, *paces_xtal, *aries_xtal;

#define N_ARIES 76
#define N_LABR 8
#define N_RCMP_POS 6
#define N_RCMP_STRIPS 32
char rcmp_strips_handles[N_RCMP_POS][32]={ "RCS01_E_strips","RCS02_E_strips","RCS03_E_strips","RCS04_E_strips","RCS05_E_strips","RCS06_E_strips" };
TH2I  *rcmp_strips[N_RCMP_POS];

int init_singles_histos(Config *cfg)
{
  char title[STRING_LEN], handle[STRING_LEN];
  int i;
  open_folder(cfg, "Hits_and_Sums");
  open_folder(cfg, "Sums");
  sprintf(title,  "Addback_Sum_Energy"); sprintf(handle, "Addback_Sum_E");
  ge_sum_ab = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_Energy"); sprintf(handle, "Ge_Sum_E");
  ge_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_En_betaTagged"); sprintf(handle, "Ge_Sum_E_B");
  ge_sum_b = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "PACES_Sum_Energy"); sprintf(handle, "Paces_Sum_E");
  paces_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "LaBr3_Sum_Energy"); sprintf(handle, "Labr_Sum_E");
  labr_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "ARIES_Sum_Energy"); sprintf(handle, "Aries_Sum_E");
  aries_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "RCMP_Sum_Energy"); sprintf(handle, "RCMP_Sum_E");
  rcmp_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "RCMP_Sum_FB_Energy"); sprintf(handle, "RCMP_Sum_FB_E");
  rcmp_fb_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Upstream_Ge_Sum_Energy"); sprintf(handle, "US_Ge_Sum_E");
  ge_sum_us = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Downstream_Ge_Sum_Energy"); sprintf(handle, "DS_Ge_Sum_E");
  ge_sum_ds = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Upstream_AB_Sum_Energy"); sprintf(handle, "US_AB_Sum_E");
  ge_sum_ab_us = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Downstream_AB_Sum_Energy"); sprintf(handle, "DS_AB_Sum_E");
  ge_sum_ab_ds = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  for(i=0; i<NUM_CLOVER; i++){
    sprintf(title,  "Addback_%d", i );
    sprintf(handle, "Addback_%d", i );
    ge_ab_e[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  }
  close_folder(cfg);
  open_folder(cfg, "Energy");
  sprintf(title, "GeEnergy_CrystalNum"); sprintf(handle, "GeEnergy_Xtal");
  ge_xtal     = H2_BOOK(cfg, handle, title, 64, 0, 64,
			                                    E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "BgoEnergy_CrystalNum"); sprintf(handle, "BgoEnergy_Xtal");
  bgo_xtal     = H2_BOOK(cfg, handle, title, 320, 0, 320,
			                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "BgoFrontEnergy_CrystalNum"); sprintf(handle, "BgoFrontEnergy_Xtal");
  bgof_xtal     = H2_BOOK(cfg, handle, title, 128, 0, 128,
			                                       E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "BgoSideEnergy_CrystalNum"); sprintf(handle, "BgoSideEnergy_Xtal");
  bgos_xtal     = H2_BOOK(cfg, handle, title, 128, 0, 128,
			                                       E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "BgoBackEnergy_CrystalNum"); sprintf(handle, "BgoBackEnergy_Xtal");
  bgob_xtal     = H2_BOOK(cfg, handle, title, 64, 0, 64,
			                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "BgoAncilEnergy_CrystalNum"); sprintf(handle, "BgoAncilEnergy_Xtal");
  bgoa_xtal     = H2_BOOK(cfg, handle, title, 32, 0, 32,
			                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "LaBr3Energy_CrystalNum"); sprintf(handle, "LabrEnergy_Xtal");
  labr_xtal     = H2_BOOK(cfg, handle, title, 16, 0, 16,
			                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "PACESEnergy_CrystalNum"); sprintf(handle, "PacesEnergy_Xtal");
  paces_xtal     = H2_BOOK(cfg, handle, title, 16, 0, 16,
                                            E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "ARIESEnergy_CrystalNum"); sprintf(handle, "AriesEnergy_Xtal");
  aries_xtal     = H2_BOOK(cfg, handle, title, 80, 0, 80,
                                            E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
    for(i=0; i<N_RCMP_POS; i++){ // Create RCMP DSSD strip spectra
      rcmp_strips[i] = H2_BOOK(cfg, rcmp_strips_handles[i], rcmp_strips_handles[i], 2*N_RCMP_STRIPS, 0, 2*N_RCMP_STRIPS,
                                            E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
    }
  close_folder(cfg);
  close_folder(cfg);
  return(0);
}

// Time difference spectra
#define DT_SPEC_LENGTH 4096
#define N_DT 18
char dt_handles[N_DT][32]={ "dt_ge_ge", "dt_ge_bgo", "dt_ge_sep", "dt_ge_zds",
                            "dt_ge_pac", "dt_ge_labr", "dt_ge_rcmp", "dt_pac_zds",
                            "dt_pac_labr", "dt_rcmp_rcmp", "dt_ge_art", "dt_labr_art",
                            "dt_paces_art", "dt_art_art", "dt_art_tac", "dt_zds_tac", "dt_labr_tac", "dt_labr_zds" };
TH1I  *dt_hist[N_DT];
TH1I  *dt_tacs_hist[N_LABR];

// Angular difference spectra
#define GE_ANG_CORR_SPEC_LENGTH 4096
#define N_GE_ANG_CORR 52
TH2I  *gg_ang_corr_hist[N_GE_ANG_CORR];

// TAC spectra
TH1I *tac_labr_hist[(int)((N_LABR)*(N_LABR-1)/2)]; // this index numbers are the LaBr-LaBr position numbers
TH1I *tac_aries_lbl_hist[N_LABR];  // this index number is the LaBr position number
TH1I *tac_aries_art_hist[N_ARIES];  // this index number is the Aries position number

// En-En Coincidence matrices
TH2I *gg, *gg_ab, *gg_opp, *gg_hit, *bgobgo_hit, *aa_hit, *gea_hit, *lba_hit, *ge_paces, *ge_labr, *ge_rcmp, *labr_labr, *labr_rcmp, *ge_art, *paces_art, *labr_art, *art_art;

char rcmp_hit_handles[N_RCMP_POS][32]={ "RCS01_PN_hits","RCS02_PN_hits","RCS03_PN_hits","RCS04_PN_hits","RCS05_PN_hits","RCS06_PN_hits" };
char rcmp_fb_handles[N_RCMP_POS][32]={ "RCS01_Front_Back","RCS02_Front_Back","RCS03_Front_Back","RCS04_Front_Back","RCS05_Front_Back","RCS06_Front_Back" }; // front-back
TH2I  *rcmp_hit[N_RCMP_POS];
TH2I  *rcmp_fb[N_RCMP_POS];

int init_coinc_histos(Config *cfg)
{
  char title[STRING_LEN], handle[STRING_LEN];
  char tmp[STRING_LEN];
  int i,j,k;

  open_folder(cfg, "Hits_and_Sums");
  open_folder(cfg, "Delta_t");
  for(i=0; i<N_DT; i++){ // Create delta-t spectra
    dt_hist[i] = H1_BOOK(cfg, dt_handles[i], dt_handles[i], DT_SPEC_LENGTH, 0, DT_SPEC_LENGTH);
  }
  for(i=0; i<N_LABR; i++){ // Create delta-t spectra separately for each TAC
  sprintf(title, "dt_labr_tac%d",i+1);
    dt_tacs_hist[i] = H1_BOOK(cfg, title, title, DT_SPEC_LENGTH, 0, DT_SPEC_LENGTH);
  }
  close_folder(cfg);
  close_folder(cfg);
  open_folder(cfg, "Coinc");
  open_folder(cfg, "Coinc");
  sprintf(title, "Addback_GG"); sprintf(handle, "Addback_GG");
  gg_ab     = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GG"); sprintf(handle, "GG");
  gg        = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GePaces"); sprintf(handle, "GePaces");
  ge_paces  = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GeLabr"); sprintf(handle, "GeLabr");
   ge_labr   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                       E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "GeAries"); sprintf(handle, "GeAries");
    ge_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
 		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "LaBrAries"); sprintf(handle, "LaBrAries");
    labr_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
   		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "AriesAries"); sprintf(handle, "AriesAries");
    art_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
     		                                        E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "LaBrLabr"); sprintf(handle, "LaBrLabr");
   labr_labr   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
 		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "GeRCMP"); sprintf(handle, "GeRCMP");
   ge_rcmp  = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
 		                                      E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
   sprintf(title, "LaBrRCMP"); sprintf(handle, "LaBrRCMP");
   labr_rcmp  = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                            E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
   sprintf(title, "GGoppo"); sprintf(handle, "GGoppo");
   gg_opp    = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                       E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   close_folder(cfg);
   open_folder(cfg, "Hits");
   sprintf(title, "GeGeHit"); sprintf(handle, "GGHit");
   gg_hit    = H2_BOOK(cfg, handle, title, 64, 0, 64,
		                                       64, 0, 64);
   sprintf(title, "BgoBgoHit"); sprintf(handle, "BgoBgoHit");
   bgobgo_hit  = H2_BOOK(cfg, handle, title, 512, 0, 512,
		                     	                   512, 0, 512);
   sprintf(title, "GeAriesHit"); sprintf(handle, "GAHit");
   gea_hit    = H2_BOOK(cfg, handle, title, 64, 0, 64,
                                            80, 0, 80);
  sprintf(title, "LaBrAriesHit"); sprintf(handle, "LAHit");
  lba_hit    = H2_BOOK(cfg, handle, title, 16, 0, 16,
                                           80, 0, 80);
   sprintf(title, "AriesAriesHit"); sprintf(handle, "AAHit");
   aa_hit    = H2_BOOK(cfg, handle, title, 80, 0, 80,
                                           80, 0, 80);
     for(i=0; i<N_RCMP_POS; i++){ // Create RCMP DSSD hitpatterns
       rcmp_hit[i] = H2_BOOK(cfg, rcmp_hit_handles[i], rcmp_hit_handles[i], N_RCMP_STRIPS, 0, N_RCMP_STRIPS,
                                                                            N_RCMP_STRIPS, 0, N_RCMP_STRIPS);
     }
for(i=0; i<N_RCMP_POS; i++){ // Create RCMP DSSD front energy vs back energy heatmaps
rcmp_fb[i] = H2_BOOK(cfg, rcmp_fb_handles[i], rcmp_fb_handles[i], E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH,
                                    E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
}
   close_folder(cfg);
   open_folder(cfg, "Ang_Corr");
   for(i=0; i<N_GE_ANG_CORR; i++){ // Create Ge-Ge angular correlation spectra
     sprintf(tmp,"Ge-Ge_angular_bin%d",i);
     gg_ang_corr_hist[i] = H2_BOOK(cfg, tmp, tmp, GE_ANG_CORR_SPEC_LENGTH, 0, GE_ANG_CORR_SPEC_LENGTH,
                                                  GE_ANG_CORR_SPEC_LENGTH, 0, GE_ANG_CORR_SPEC_LENGTH);
   }
   close_folder(cfg);
   open_folder(cfg, "LBL_TACs");
   k=-1;
   for(i=1; i<=N_LABR; i++){ // Create TAC coincidence pair spectra
     for(j=(i+1); j<=N_LABR; j++){
       k++;
     sprintf(tmp,"TAC_%d_%d",i,j);
     tac_labr_hist[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
     }
   }

   close_folder(cfg);
   open_folder(cfg, "ART_TACs");
 sprintf(tmp,"TAC-ARIES-LaBr3-1275keV");
 aries_tac = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
   for(i=0; i<N_LABR; i++){ // Create TAC coincidence with ARIES spectra
     sprintf(tmp,"TAC_ART_LBL%d",i+1);
     tac_aries_lbl_hist[i] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
   }
     for(i=0; i<N_ARIES; i++){ // Create TAC coincidence with ARIES spectra
       sprintf(tmp,"TAC_LBL_ART%d",i+1);
       tac_aries_art_hist[i] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
     }

   close_folder(cfg);
   close_folder(cfg);
   return(0);
}

int fill_singles_histos(Grif_event *ptr)
{
  int i, j, dt, pos, sys, elem, clover, ge_addback_gate = 25, ge_sum_gate = 25;
  char *name, c;
  long ts;

  sys = ptr->subsys;
  if( sys >=0 && sys < MAX_SUBSYS ){
    if( mult_hist[sys] != NULL ){ mult_hist[sys]->Fill(mult_hist[sys], ptr->fold, 1);  }
  }
  //pos = ptr->array_posn;
  pos  = crystal_table[ptr->chan];

  // Here we should not use dtype because the dtype can change between experiments.
  // Here we should use the Subsytem name, which we can get from dtype.
  // The dtype->subsystem mapping is done from the PSC table

   switch (sys){
   case SUBSYS_HPGE: // GRGa
     // Only use GRGa
     if(output_table[ptr->chan] == 1){

       //  ge-crystal-sum
       if( pos >= 0 && pos < 64 ){
         ge_sum->Fill(ge_sum, (int)ptr->ecal, 1);
         ge_xtal->Fill(ge_xtal, pos, (int)ptr->ecal, 1);

	 // Separate sum spectra for upstream and downstream
	 if(pos<33){
	   ge_sum_us->Fill(ge_sum_us, (int)ptr->ecal, 1);
	 }else{
	   ge_sum_ds->Fill(ge_sum_ds, (int)ptr->ecal, 1);
	 }

	 // Ge-Addback
	 clover = (int)(pos/16)+1;
         if( clover >= 0 && clover < NUM_CLOVER && ptr->esum >= 0 ){   // ge addback
           ge_ab_e[clover]->Fill(ge_ab_e[clover],  (int)ptr->esum, 1);
           ge_sum_ab   ->Fill(ge_sum_ab,     (int)ptr->esum, 1);

	   // Separate Addback sum spectra for upstream and downstream
	   if(clover<9){
	     ge_sum_ab_us->Fill(ge_sum_ab_us, (int)ptr->ecal, 1);
	   }else{
	     ge_sum_ab_ds->Fill(ge_sum_ab_ds, (int)ptr->ecal, 1);
	   }
         }
	 }else {
	 fprintf(stderr,"bad ge crystal[%d] for chan %d\n", pos, ptr->chan);
       }
     }
     break;
   case SUBSYS_BGO: // BGOs
      pos  = crystal_table[ptr->chan];
      elem = element_table[ptr->chan];
      if( pos < 0 || pos > 64 ){
         fprintf(stderr,"bad bgo crystal[%d] for chan %d\n", pos, ptr->chan);
      } else if( elem < 1 || elem > 5 ){
         fprintf(stderr,"bad bgo element[%d] for chan %d, %s, subsys %d\n", elem, ptr->chan, chan_name[ptr->chan], ptr->subsys);
      } else {
         pos *= 5; pos += (elem-1);
         bgo_xtal->Fill(bgo_xtal, pos, (int)ptr->ecal, 1);
         if(elem <3){ // front
            pos  = crystal_table[ptr->chan];
            pos *= 2; pos += (elem-1);
            bgof_xtal->Fill(bgof_xtal, pos, (int)ptr->ecal, 1);
         } else if(elem>4){ // back
            pos  = crystal_table[ptr->chan];
            bgob_xtal->Fill(bgob_xtal, pos, (int)ptr->ecal, 1);
         } else{ // side
            pos  = crystal_table[ptr->chan];
            pos *= 2; pos += (elem-3);
            bgos_xtal->Fill(bgos_xtal, pos, (int)ptr->ecal, 1);
         }
      }
      break;
   case SUBSYS_LABR_BGO: // Ancillary BGOs
      pos  = crystal_table[ptr->chan];
      elem = element_table[ptr->chan];
      if( pos < 0 || pos > 8 ){
         fprintf(stderr,"bad ancillary bgo crystal[%d] for chan %d\n", pos, ptr->chan);
      } else if( elem < 1 || elem > 3 ){
         fprintf(stderr,"bad ancillary bgo element[%d] for chan %d\n", elem, ptr->chan);
      } else {
         pos *= 3; pos += (elem-1);
         bgoa_xtal->Fill(bgoa_xtal, pos, (int)ptr->ecal, 1);
      }
      break;
   case SUBSYS_PACES: // PACES
    paces_sum->Fill(paces_sum, (int)ptr->ecal, 1);
    pos  = crystal_table[ptr->chan];
    if( pos < 0 || pos > 5 ){
      fprintf(stderr,"bad ancillary PACES crystal[%d] for chan %d\n", pos, ptr->chan);
    } else {
      paces_xtal->Fill(paces_xtal, pos, (int)ptr->ecal, 1);
    }
    break;
   case SUBSYS_LABR_L: // LaBr3 (LBL)
    labr_sum->Fill(labr_sum, (int)ptr->ecal, 1);
      pos  = crystal_table[ptr->chan];
      if( pos < 0 || pos > 8 ){
         fprintf(stderr,"bad ancillary LaBr3 crystal[%d] for chan %d\n", pos, ptr->chan);
      } else {
         labr_xtal->Fill(labr_xtal, pos, (int)ptr->ecal, 1);
      }
      break;
   case SUBSYS_SCEPTAR:
   break;
   case SUBSYS_LABR_T:
   break;
   case SUBSYS_DESCANT:
   break;
   case SUBSYS_ARIES: // ARIES
    aries_sum->Fill(aries_sum, (int)ptr->ecal, 1);
    pos  = crystal_table[ptr->chan];
    if( pos < 1 || pos > 76 ){
      fprintf(stderr,"bad aries crystal[%d] for chan %d\n", pos, ptr->chan);
    } else {
      aries_xtal->Fill(aries_xtal, pos, (int)ptr->ecal, 1);
    }
   break;
   case SUBSYS_ZDS:
   break;
   case SUBSYS_RCMP:
       rcmp_sum->Fill(rcmp_sum, (int)ptr->ecal, 1);
       if(ptr->esum>0){
         rcmp_fb_sum->Fill(rcmp_fb_sum, (int)ptr->esum, 1);
       }
       pos  = crystal_table[ptr->chan];
       elem = element_table[ptr->chan];

       /*
       if(pos>0 && pos<4){
         elem = reorder_rcmp_DS_strips(elem);
       }else if(pos>3 && pos<7){
         elem = reorder_rcmp_US_strips(elem);
       }
       */
       elem = (int)(elem + (int)(polarity_table[ptr->chan]*N_RCMP_STRIPS)); // polarity_table value is 0 or 1
       if( pos < 1 || pos > 6 ){
         fprintf(stderr,"bad RCMP DSSD[%d] for chan %d, elem%d, pol%d\n", pos, ptr->chan, elem, polarity_table[ptr->chan]);
       } else if( elem < 0 || elem > 63 ){
         fprintf(stderr,"bad RCMP strip[%d] for chan %d, pos%d, pol%d\n", elem, ptr->chan, pos, polarity_table[ptr->chan]);
       } else {
         rcmp_strips[(pos-1)]->Fill(rcmp_strips[(pos-1)], elem, (int)ptr->ecal, 1);
       }
       break;
   default: // Unrecognized or unprocessed dtype
      break;
   }// end of switch

   return(0);
}



int frag_hist[MAX_COINC_EVENTS];
int fill_coinc_histos(int win_idx, int frag_idx)
{
  Grif_event *alt, *ptr = &grif_event[win_idx];
  int dt, abs_dt, i, gg_gate=25;
  //int g_rcmp_lower_gate=22, g_rcmp_upper_gate=68;
  int g_aries_upper_gate=25;
  int g_aries_tac_gate=60;
  int g_rcmp_upper_gate=75;
  int pos, c1, c2, index;
  int global_window_size = 100; // size in grif-replay should be double this
  int corrected_tac_value;
  int tac_offset[8] = {-7300,-5585,-6804,0,-6488,-5682,-5416,0};

  // histogram of coincwin-size
  dt = (frag_idx - win_idx + 2*MAX_COINC_EVENTS) %  MAX_COINC_EVENTS;
  ++frag_hist[dt];
  if( win_idx == frag_idx ){ return(0); } // window size = 1 - no coinc
  if( (i=win_idx+1) == MAX_COINC_EVENTS ){ i = 0; } // wrap
  // check all conicidences in window
  while( 1 ){
    alt = &grif_event[i];
    abs_dt = dt = ptr->ts - alt->ts; if( dt < 0 ){ abs_dt = -1*dt; }
    if( abs_dt > global_window_size ){ break; }

    switch(ptr->subsys){
      case SUBSYS_HPGE:

      // Only use GRGa
      if(output_table[ptr->chan] == 1){

        switch(alt->subsys){
          case SUBSYS_HPGE:
          // Only use GRGa
          if(output_table[alt->chan] == 1){

            dt_hist[0]->Fill(dt_hist[0], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // ge-ge
            // Ge-Ge matrices
            if( alt->dtype == 0 && abs_dt < gg_gate ){                 // ge-ge (and addback)
              gg->Fill(gg, (int)ptr->ecal, (int)alt->ecal, 1);
              if( ptr->esum >= 0 &&  alt->esum >= 0 ){
                gg_ab->Fill(gg_ab, (int)ptr->esum, (int)alt->esum, 1);
              }

              c1 = crystal_table[ptr->chan];
              if( c1 >= 0 && c1 < 64 ){
                c2 = crystal_table[alt->chan];
                if( c2 >= 0 && c2 < 64 ){
                  gg_hit->Fill(gg_hit, c1, c2, 1);

                  // Ge-Ge angular correlations
                  // Fill the appropriate angular bin spectrum
                  index = ge_angles_145mm[c1][c2];
                  gg_ang_corr_hist[index]->Fill(gg_ang_corr_hist[index], (int)ptr->ecal, (int)alt->ecal, 1);

                  // Ge-Ge with 180 degrees between Ge1 and Ge2
                  if( c2 == grif_opposite[c1] ){
                    gg_opp->Fill(gg_opp, (int)ptr->ecal, (int)alt->ecal, 1);
                  }
                }
              }
            }
          } break;
          case SUBSYS_BGO:     // ge-bgo
          dt_hist[1]->Fill(dt_hist[1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          break;
          case SUBSYS_PACES: // ge-paces
          dt_hist[4]->Fill(dt_hist[4], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < gg_gate ){
            ge_paces  ->Fill(ge_paces, (int)ptr->ecal, (int)alt->ecal, 1);
          } break;
          case SUBSYS_LABR_L:  // ge-labr
          dt_hist[5]->Fill(dt_hist[5], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < gg_gate ){
            ge_labr->Fill(ge_labr, (int)ptr->ecal, (int)alt->ecal, 1);
          } break;
          case SUBSYS_SCEPTAR:  // ge-sep
          dt_hist[2]->Fill(dt_hist[2], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < gg_gate ){
            ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
          }
          break;
          case SUBSYS_ARIES: // ge-aries
          dt_hist[10]->Fill(dt_hist[10], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < g_aries_upper_gate ){
            ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
            ge_art->Fill(ge_art, (int)ptr->ecal, (int)alt->esum, 1);
            c1 = crystal_table[ptr->chan];
            if(c1 >= 1 && c1 <=64 ){
              c2 = crystal_table[alt->chan];
              if(c2 >= 1 && c2 <=76 ){
                gea_hit->Fill(gea_hit, c1, c2, 1);
              }
            }
          }
          break;
          case SUBSYS_ZDS: // ge-zds
          dt_hist[3]->Fill(dt_hist[3], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < gg_gate ){
            ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
          }
          break;
          case SUBSYS_RCMP: // ge-rcmp
          dt_hist[6]->Fill(dt_hist[6], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < g_rcmp_upper_gate && alt->esum>0){
            ge_rcmp->Fill(ge_rcmp, (int)ptr->ecal, (int)alt->esum, 1);
          }
          break;

          default: break; // unprocessed coincidence combinations
        } // end of inner switch(ALT) for ptr=HPGe
      }
      break; // outer-switch-case-GE

      case SUBSYS_SCEPTAR:  // sceptar matrices
      switch(alt->subsys){
        case SUBSYS_HPGE: // sep-ge
        // Only use GRGa
        if(output_table[alt->chan] == 1){
          dt_hist[2]->Fill(dt_hist[2], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < gg_gate ){
            ge_sum_b->Fill(ge_sum_b, (int)alt->ecal, 1); // beta-gated Ge sum energy spectrum
          }
        }
        break;
        default: break;
      }
      break;

      case SUBSYS_PACES: // paces matrices
      switch(alt->subsys){
        case SUBSYS_ZDS:
        dt_hist[7]->Fill(dt_hist[7], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-zds
        break;
        case SUBSYS_HPGE:
        // Only use GRGa
        if(output_table[alt->chan] == 1){
          dt_hist[4]->Fill(dt_hist[4], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-ge
        }
        break;
        case SUBSYS_LABR_L:
        dt_hist[8]->Fill(dt_hist[8], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-labr
        break;
        case SUBSYS_ARIES:
        dt_hist[13]->Fill(dt_hist[13], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-aries
        if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0){
          paces_art->Fill(paces_art, (int)ptr->ecal, (int)alt->ecal, 1);
        }
        break;
        default: break; // unprocessed coincidence combinations
      } // end of inner switch(ALT)

      case SUBSYS_LABR_L: // Labr matrices
      switch(alt->subsys){
        case SUBSYS_HPGE:
        // Only use GRGa
        if(output_table[alt->chan] == 1){
          dt_hist[5]->Fill(dt_hist[5], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // labr-ge
        }
        break;
        case SUBSYS_LABR_L:
        dt_hist[8]->Fill(dt_hist[8], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // labr-labr
        break;
        case SUBSYS_ZDS: // labr-zds
        dt_hist[17]->Fill(dt_hist[17], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        break;
        case SUBSYS_ARIES:
        dt_hist[12]->Fill(dt_hist[12], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // laBr-aries
        if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0){
          labr_art->Fill(labr_art, (int)ptr->ecal, (int)alt->ecal, 1);

          c1 = crystal_table[ptr->chan];
          if(c1 >= 1 && c1 <=8 ){
            c2 = crystal_table[alt->chan];
            if(c2 >= 1 && c2 <=76 ){
          lba_hit->Fill(lba_hit, c1, c2, 1);
        }
      }
        }
        break;
        case SUBSYS_LABR_T: // labr-tac
        dt_hist[16]->Fill(dt_hist[16], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        dt_tacs_hist[crystal_table[alt->chan]-1]->Fill(dt_tacs_hist[crystal_table[alt->chan]-1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        if(crystal_table[alt->chan] == 8){ // ARIES TAC

          corrected_tac_value = (int)alt->ecal + tac_offset[crystal_table[ptr->chan]-1];
          tac_aries_lbl_hist[crystal_table[ptr->chan]-1]->Fill(tac_aries_lbl_hist[crystal_table[ptr->chan]-1], corrected_tac_value, 1); // tac spectrum per LBL
          if(ptr->ecal >1225 && ptr->ecal <1315){ // gate on LaBr3 energy 1275keV
            aries_tac->Fill(aries_tac, (int)corrected_tac_value, 1); // tac spectrum gated on 1332keV
          }
        }

        break;
        default: break; // unprocessed coincidence combinations
      } // end of inner switch(ALT)
      break;

      case SUBSYS_BGO: // bgo matrices
      if(alt->subsys == SUBSYS_BGO && abs_dt < gg_gate ){ // bgo-bgo
        c1 = crystal_table[ptr->chan];
        c2 = crystal_table[alt->chan];
        bgobgo_hit->Fill(bgobgo_hit, c1, c2, 1);
      } break;

      case SUBSYS_RCMP: // rcmp matrices
      //fprintf(stdout,"RCMP coinc. with %d, with dt=%d\n",alt->subsys,abs_dt);
      switch(alt->subsys){
        case SUBSYS_HPGE: // rcmp-ge
        // Only use GRGa
        if(output_table[alt->chan] == 1){

          dt_hist[6]->Fill(dt_hist[6], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          // if( (abs_dt > g_rcmp_lower_gate && abs_dt < g_rcmp_upper_gate) && ptr->ecal>0){
          if( ( abs_dt < g_rcmp_upper_gate) && ptr->ecal>0){
            ge_rcmp->Fill(ge_rcmp, (int)alt->ecal, (int)ptr->ecal, 1);
          }
        }
        break;
        case SUBSYS_RCMP: // rcmp-rcmp
        // fprintf(stdout,"RCMP coinc. with %d, with dt=%d\n",alt->subsys,abs_dt);
        dt_hist[9]->Fill(dt_hist[9], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // rcmp-rcmp - This might need an extra condition for Polarity difference
        pos = crystal_table[ptr->chan];
        if(pos == crystal_table[alt->chan] && (polarity_table[ptr->chan] != polarity_table[alt->chan])){ // front and back of same DSSD
          c1 = element_table[ptr->chan];
          c2 = element_table[alt->chan];

          /*
          if(pos>0 && pos<4){
          c1 = reorder_rcmp_DS_strips(c1);
          c2 = reorder_rcmp_DS_strips(c2);
        }else if(pos>3 && pos<7){
        c1 = reorder_rcmp_US_strips(c1);
        c2 = reorder_rcmp_US_strips(c2);
      }
      */
      rcmp_hit[(pos-1)]->Fill(rcmp_hit[(pos-1)], c1, c2, 1);
      rcmp_fb[(pos-1)]->Fill(rcmp_fb[(pos-1)], (int)ptr->ecal, (int)alt->ecal, 1);
    }
    break;
    default: break;
  } // end of inner switch(ALT)
  break; // end of inner switch(ALT) for ptr=rcmp

  case SUBSYS_ARIES: // aries matrices
  switch(alt->subsys){
    case SUBSYS_HPGE: // aries-ge
    // Only use GRGa
    if(output_table[alt->chan] == 1){

      dt_hist[10]->Fill(dt_hist[10], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0 && alt->ecal>0){
        ge_art->Fill(ge_art, (int)alt->ecal, (int)ptr->ecal, 1);
        c1 = crystal_table[ptr->chan];
        if(c1 >= 1 && c1 <=76 ){
          c2 = crystal_table[alt->chan];
          if(c2 >= 1 && c2 <=64 ){
        gea_hit->Fill(gea_hit, c2, c1, 1);
      }
    }
      }
    }
    break;
    case SUBSYS_LABR_L: // aries-labr
    dt_hist[11]->Fill(dt_hist[11], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0){
      labr_art->Fill(labr_art, (int)alt->ecal, (int)ptr->ecal, 1);
      c1 = crystal_table[alt->chan];
      if(c1 >= 1 && c1 <=8 ){
        c2 = crystal_table[ptr->chan];
        if(c2 >= 1 && c2 <=76 ){
      lba_hit->Fill(lba_hit, c1, c2, 1);
    }
  }
    }
    break;
    case SUBSYS_PACES: // aries-paces
    dt_hist[12]->Fill(dt_hist[12], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0){
      paces_art->Fill(paces_art, (int)alt->ecal, (int)ptr->ecal, 1);
    }
    break;
    case SUBSYS_ARIES: // aries-aries
    dt_hist[13]->Fill(dt_hist[13], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    c1 = crystal_table[ptr->chan];
    if(c1 >= 1 && c1 <=76 ){
      c2 = crystal_table[alt->chan];
      if(c2 >= 1 && c2 <=76 ){

        aa_hit->Fill(aa_hit, c1, c2, 1);
        art_art->Fill(art_art, (int)ptr->ecal, (int)alt->ecal, 1);
      }
    }
    break;
    case SUBSYS_LABR_T: // aries-tac
        if(crystal_table[alt->chan] == 8){ // ARIES TAC
        dt_hist[14]->Fill(dt_hist[14], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          tac_aries_art_hist[crystal_table[ptr->chan]-1]->Fill(tac_aries_art_hist[crystal_table[ptr->chan]-1], (int)alt->ecal, 1); // tac spectrum per ART tiles
        }
    break;
    default: break;
  } // end of inner switch(ALT)
  break; // end of inner switch(ALT) for ptr=rcmp


  case SUBSYS_ZDS: // zds matrices
  switch(alt->subsys){
    case SUBSYS_LABR_T: // zds-tac
    dt_hist[15]->Fill(dt_hist[15], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    break;
    case SUBSYS_HPGE: // zds-HPGe
    // Only use GRGa
    if(output_table[alt->chan] == 1){
      dt_hist[3]->Fill(dt_hist[3], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    }
    break;
    case SUBSYS_LABR_L: // zds-labr
    dt_hist[17]->Fill(dt_hist[17], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    break;
    default: break;
  } // end of inner switch(ALT)
  break;

  case SUBSYS_LABR_T: // tac matrices
  switch(alt->subsys){
    case SUBSYS_ZDS: // tac-zds
    dt_hist[15]->Fill(dt_hist[15], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    break;
    case SUBSYS_ARIES: // tac-aries
    if(crystal_table[ptr->chan] == 8){ // ARIES TAC
    dt_hist[14]->Fill(dt_hist[14], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      tac_aries_art_hist[crystal_table[alt->chan]-1]->Fill(tac_aries_art_hist[crystal_table[alt->chan]-1], (int)ptr->ecal, 1); // tac spectrum per ART tiles
    }
    break;
    case SUBSYS_LABR_L: // tac-labr
    dt_hist[16]->Fill(dt_hist[16], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    dt_tacs_hist[crystal_table[ptr->chan]-1]->Fill(dt_tacs_hist[crystal_table[ptr->chan]-1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    if(crystal_table[ptr->chan] == 8){ // ARIES TAC

      corrected_tac_value = (int)ptr->ecal + tac_offset[crystal_table[alt->chan]-1];
      tac_aries_lbl_hist[crystal_table[alt->chan]-1]->Fill(tac_aries_lbl_hist[crystal_table[alt->chan]-1], corrected_tac_value, 1); // tac spectrum per LBL
      if(alt->ecal >1225 && alt->ecal <1315){ // gate on LaBr3 energy 1275keV

        aries_tac->Fill(aries_tac, corrected_tac_value, 1); // tac spectrum gated on 1332keV
      }
    }
    break;
    default: break;
  } // end of inner switch(ALT)
  break;

  default: break; // Unrecognized or unprocessed subsys type
}// end of switch(ptr)
if( i == frag_idx ){ break; }
if( ++i == MAX_COINC_EVENTS ){ i = 0; } // wrap
}// end of while
return(0);
}

int reorder_rcmp_US_strips(int c1){
// Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25966
switch(c1){

  case  0: c1 = 1; break;
  case  1: c1 = 3; break;
  case  2: c1 = 5; break;
  case  3: c1 = 7; break;
  case  4: c1 = 9; break;
  case  5: c1 = 11; break;
  case  6: c1 = 13; break;
  case  7: c1 = 15; break;
  case  8: c1 = 31; break;
  case  9: c1 = 29; break;
  case 10: c1 = 27; break;
  case 11: c1 = 25; break;
  case 12: c1 = 23; break;
  case 13: c1 = 21; break;
  case 14: c1 = 19; break;
  case 15: c1 = 17; break;
  case 16: c1 = 16; break;
  case 17: c1 = 18; break;
  case 18: c1 = 20; break;
  case 19: c1 = 22; break;
  case 20: c1 = 24; break;
  case 21: c1 = 26; break;
  case 22: c1 = 28; break;
  case 23: c1 = 30; break;
  case 24: c1 = 14; break;
  case 25: c1 = 12; break;
  case 26: c1 = 10; break;
  case 27: c1 = 8; break;
  case 28: c1 = 6; break;
  case 29: c1 = 4; break;
  case 30: c1 = 2; break;
  case 31: c1 = 0; break;

  default: break;
}
return(c1);
}

int reorder_rcmp_DS_strips(int c1){
// Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25968
switch(c1){

  case  0: c1 = 0; break;
  case  1: c1 = 2; break;
  case  2: c1 = 4; break;
  case  3: c1 = 6; break;
  case  4: c1 = 8; break;
  case  5: c1 = 10; break;
  case  6: c1 = 12; break;
  case  7: c1 = 14; break;
  case  8: c1 = 30; break;
  case  9: c1 = 28; break;
  case 10: c1 = 26; break;
  case 11: c1 = 24; break;
  case 12: c1 = 22; break;
  case 13: c1 = 20; break;
  case 14: c1 = 18; break;
  case 15: c1 = 16; break;
  case 16: c1 = 17; break;
  case 17: c1 = 19; break;
  case 18: c1 = 21; break;
  case 19: c1 = 23; break;
  case 20: c1 = 25; break;
  case 21: c1 = 27; break;
  case 22: c1 = 29; break;
  case 23: c1 = 31; break;
  case 24: c1 = 15; break;
  case 25: c1 = 13; break;
  case 26: c1 = 11; break;
  case 27: c1 = 9; break;
  case 28: c1 = 7; break;
  case 29: c1 = 5; break;
  case 30: c1 = 3; break;
  case 31: c1 = 1; break;

  default: break;
}
return(c1);
}

//#######################################################################
//###########   READ XML ODB DUMP FROM START OF DATA FILE   #############
//#######################################################################

static char   path[256];
static char dirname[64],value[32],type[32];
extern char midas_runtitle[STRINGLEN];

static void *arrayptr;
int read_odb_items(int len, int *bank_data)
{
   char *path_ptr, *ptr, *str, *odb_data = (char *)bank_data, posn[2];
   int i, c = '<', d = '>', dtype=0, active=0, index=0;
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
      } else if( strncmp(ptr,"</dir>",6) == 0 ){
         while(1){
            if( --path_ptr < path ){ path_ptr = path;  *path_ptr = '\0';  break; }
            if( *path_ptr == '/' ){ *path_ptr = '\0';  break; }
         }
         index=0; // for debugger to stop here
      } else if( strncasecmp(ptr,"<key name=\"Run Title\" type=\"STRING\"", 35) == 0 ){
         ptr = str+1;
         if( (str = strchr(ptr,c)) == NULL ){ break; }
         i = (str-ptr) > STRINGLEN-1 ? STRINGLEN-1 : (str-ptr);
         memcpy( midas_runtitle, ptr, i ); midas_runtitle[i] = 0;
         ptr += i+1;
         if( (str = strchr(ptr,d)) == NULL ){ break; }
      } else if( strncmp(ptr,"</keyarray>",10) == 0 ){ active = 0; arrayptr = '\0';
      } else if( strncmp(ptr,"<keyarray ",10) == 0 ){
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

   // arrays typically around 500 entries [one per "chan"] each entry with ...
   //   daq-address, name, type, gains etc.
   //
   gen_derived_odb_tables();
   return(0);
}

int gen_derived_odb_tables()
{
  char sys_name[64], crystal, polarity, type;
  int i, j, tmp, pos, element;
  // Also require Ge crystal numbers - which cannot be calculated from
  // data-fragment [only contains array-posn, which is clover number]
  // so calculate them here ...
  //
  //
  memset(crystal_table,  0xff, MAX_DAQSIZE*sizeof(int)); // initialise all to -1
  memset(element_table,  0xff, MAX_DAQSIZE*sizeof(int));
  memset(polarity_table, 0xff, MAX_DAQSIZE*sizeof(int));
  memset(output_table,   0xff, MAX_DAQSIZE*sizeof(int));
  memset(subsys_dtype_mat,  0,       16*16*sizeof(int));
  for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
    if( (tmp=sscanf(chan_name[i], "%3c%d%c%c%d%c", &sys_name, &pos,
    &crystal, &polarity, &element, &type)) != 6 ){
      fprintf(stderr,"can't decode name[%s] decoded %d of 6 items\n", chan_name[i], tmp );
      continue;
    }

    // Determine Polarity
    // 1 is N, 0 is P or T or S, -1 is anything else
    if(        polarity == 'N' ){ polarity_table[i] = 1;
    } else if( polarity == 'P' ){ polarity_table[i] = 0;
    } else if( polarity == 'T' ){ polarity_table[i] = 0; // TAC signal
    } else if( polarity == 'S' ){ polarity_table[i] = 1; // ARIES Standard signal
    } else if( polarity == 'X' ){ polarity_table[i] = 0; // XXX type
    } else { fprintf(stderr,"unknown polarity[=%c] in %s\n", polarity, chan_name[i]); }
    polarity_table[i] = polarity=='N' ? 1 : (polarity=='P' ? 0 : -1); // Looks like we are doing this twice each time, this is a repeat of the above

    // Determine Output
    // Some detector elements have more than one output (HPGe A and B)
    // 1 is A, 0 is B, -1 is X or unknown
    output_table[i] = type=='A' ? 1 : (type=='B' ? 0 : -1);

    // Determine crystal and element position numbers for each Subsystem
    //if(dtype_table[i] == 3 || dtype_table[i] == 8 || dtype_table[i] == 5){ // LBL and LBS, LaBr3 and ancillary BGOs, PAC paces
    if((strncmp(sys_name,"ART",3) == 0) || (strncmp(sys_name,"LBL",3) == 0) || (strncmp(sys_name,"LBS",3) == 0) || (strncmp(sys_name,"LBT",3) == 0)){ // LBL and LBS, LaBr3 and ancillary BGOs, PAC paces
      crystal_table[i] = pos;
      if(        crystal == 'A' ){ element_table[i] = 1;
      } else if( crystal == 'B' ){ element_table[i] = 2;
      } else if( crystal == 'C' ){ element_table[i] = 3;
      } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3
      } else {
        fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
      }
      //}else if(dtype_table[i] == 10 || dtype_table[i] == 11){ // RCSn and RCSp, RCMP
    }else if(strncmp(sys_name,"RCS",3) == 0){ // RCSn and RCSp, RCMP
      crystal_table[i] = pos;
      element_table[i] = element;
    }else{ // GRG and BGO
      element_table[i] = element;
      pos -= 1; pos *=4;
      if(        crystal == 'B' ){ crystal_table[i] = pos;
      } else if( crystal == 'G' ){ crystal_table[i] = pos+1;
      } else if( crystal == 'R' ){ crystal_table[i] = pos+2;
      } else if( crystal == 'W' ){ crystal_table[i] = pos+3;
      } else if( crystal == 'X' ){ crystal_table[i] = -1; // crystal undefined
      } else {
        fprintf(stderr,"unknown crystal[=%c] in %s\n", crystal, chan_name[i]);
      }
    }

    // Handle bad detector types
    if( dtype_table[i] < 0 || dtype_table[i] >= 16 ){
      fprintf(stderr,"bad datatype[%d] at table position %d\n", dtype_table[i], i);
      continue;
    }

    // Build map of names and dtypes
    for(j=0; j<MAX_SUBSYS; j++){
      if( strncmp(sys_name, subsys_handle[j], 3) == 0 ){
        ++subsys_dtype_mat[j][dtype_table[i]]; break;
      }
    }
    if( j == MAX_SUBSYS ){
      fprintf(stderr,"Unknown subsystem[%s] in %s\n", sys_name, chan_name[i]);
    }
  }

  // list of addresses. array index is channel number
  memset(address_chan, 0xFF, sizeof(address_chan)); // set to -1
  for(i=0; i<MAX_ADDRESS && i<odb_daqsize; i++){
    address_chan[ chan_address[i] ] = i;
  }

  // check the Subsytem to dtype mapping
  // This method finds the most common datatype for each subsys and warns if there is more then one.
  // This does not allow more than one detector type per subsytem
/*
  for(j=0; j<MAX_SUBSYS; j++){ // j:subsystem
  tmp = -1;
  for(i=0; i<MAX_SUBSYS; i++){ // i:datatype
  if( subsys_dtype_mat[j][i] == 0 ){ continue; }
  if( tmp == -1 ){ tmp = i; continue; }
  fprintf(stderr, "multiple datatypes{type=%d[%d],type=%d[%d]} for subsystem %s ... ",
  i, subsys_dtype_mat[j][i], tmp, subsys_dtype_mat[j][tmp], subsys_handle[j]);
  if( subsys_dtype_mat[j][i] > subsys_dtype_mat[j][tmp] ){ tmp = i; }
}
if( (subsys_dtype[j] = tmp) != -1 ){
dtype_subsys[tmp] = j;
fprintf(stdout,"Subsystem %s[dtype=%d] used in %d channels\n",
subsys_handle[j], tmp, subsys_dtype_mat[j][tmp]);
}
}
*/

// check the Subsytem to dtype mapping
// This method finds the most common subsys for each datatype and warns if there is more then one.
// This naturally allows more than one detector type per subsytem which is required for GRG and RCS.


for(j=0; j<MAX_SUBSYS; j++){ // j:datatype
   tmp = -1; dtype_subsys[j] = MAX_SUBSYS-1;
  for(i=0; i<MAX_SUBSYS; i++){ // i:subsystem
    if( subsys_dtype_mat[i][j] == 0 ){ continue; }
    if( tmp == -1 ){ tmp = i; continue; }
    fprintf(stderr,"ERROR: multiple subsystems{%s[%d],%s[%d]} for datatype %d ... ",
    subsys_handle[i], subsys_dtype_mat[i][j], subsys_handle[tmp], subsys_dtype_mat[tmp][j], j);
    if( subsys_dtype_mat[i][j] > subsys_dtype_mat[tmp][j] ){ tmp = i; }
  }
  if( (subsys_dtype[j] = tmp) != -1 ){
    dtype_subsys[j] = tmp;
    fprintf(stdout,"Datatype %d[subsystem=%s] used in %d channels\n", j, subsys_handle[tmp], subsys_dtype_mat[tmp][j]);
  }
}


//memset(address_clover, 0xFF, sizeof(address_clover)); // set to -1
//for(i=0; i<odb_daqsize; i++){
//   if( chan_address[i] >= 0 && chan_address[i] < MAX_ADDRESS ){
//	   if(strncmp("GRG",chan_name[i],3)==0 && strncmp("A",chan_name[i]+strlen(chan_name[i])-1,1)==0){
//	      strncpy(posn,chan_name[i]+3,2);
//	      address_clover[ chan_address[i] ] = atoi(posn);
//      }
//   }
//}
return(0);
}
