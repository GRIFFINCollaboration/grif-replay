#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grif-format.h"
#include "config.h"
#include "histogram.h"

// example user sort
//    GeE with BGO>0 and dt<100         [BGO rejected GeE]
//    GeE gated on background region of dt
//    GeE-GeE dt [this is already a variable]
//

// get rid of pointers in config structures
//    stuff is only ever added/removed from browser
//    so plenty of time for lists to be kept in order
// get rid of sortvar->histolist
// gates do contain condition pointers *Need to regenerate from names*
// gates also contain sortvar pointers BUT sortvars are fixed until recompile
// histos contain gatelist and gate_name pointers - regen gatelist
//                             gatename - wait and see


//#######################################################################
//#####              USER SORT - custom user histograms             #####
//#######################################################################

// fragment at win_start is leaving coincwin (frag[win_end+1] is not in coinc)
// sort pairs win_strt:win_strt+1 to  win_strt:win_end
extern Grif_event grif_event[PTR_BUFSIZE];
int user_sort(int win_strt, int win_end, int flag)
{
   Grif_event *alt, *tmp, *ptr = &grif_event[win_strt];
   Config *cfg = configs[1]; // ** Sort is using config[1]
   int i, j, abs_dt, dt, win_idx;
   Sortvar *var, *yvar;
   float value, yval;
   Histogram *histo;
   Gate *gate;
   Cond *cond;

   //if( cfg->nuser == 0 ){ return(0); } // put this in caller

   if( win_strt != win_end ){ // multiple fragments in window
      for(i=0; i<cfg->nuser; i++){ cfg->user_histos[i]->done_flag = 0; }
   } 

   if( (win_idx = win_strt+1) == PTR_BUFSIZE ){ win_idx = 0; } 
   if( ++win_end == PTR_BUFSIZE ){ win_end = 0; } // need to include win_end
   while( win_idx != win_end ){
      alt = &grif_event[win_idx];
      if( ++win_idx == PTR_BUFSIZE ){ win_idx = 0; }

      abs_dt = dt = ptr->ts - alt->ts; if( dt < 0 ){ abs_dt = -1*dt; }
      // if( abs_dt > global_window_size ){ break; }

      // clear gates
      for(i=0; i<cfg->ngates; i++){ gate = cfg->gatelist[i];
         if( gate->use_count == 0 ){ continue; }
         // for(j=0; j<gate->nconds; j++){ gate->conds[i]->passed = 0; }
         gate->valid=0;
      }
      for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
         cond->valid = cond->passed = 0;
      }

      for(j=0; j<cfg->nuser; j++){
         histo = cfg->user_histos[j];  var = histo->xvar;
         if( var->subsys1 != ptr->subsys ){ continue; }
         for(i=0; i<histo->num_gates; i++){ gate = histo->gatelist[i];
            if( ! gate->valid ){ test_gate(ptr, alt, gate); }
            if( gate->passed == 0 ){ break; }
         }
         if( i < histo->num_gates ){ continue; } // gates not all passed
         if( histo->type == INT_1D ){
            value = var->get_value(ptr, alt, var->subsys1, var->subsys2);
            histo->Fill(histo, value, 1);
            histo->done_flag = 1;
         } else if( histo->type == INT_2D ){ yvar = histo->yvar;
            // 2d histos - already checked all gates + set valids
            if( yvar->subsys1 != alt->subsys ){ continue; }
            yval = var->get_value(alt, ptr, yvar->subsys1, yvar->subsys2);
            ((TH2I *)histo)->Fill((TH2I *)histo, value, yval, 1);
            // do not set done_flag for true 2d histos
         }
      }
   }
   user_removefrom_window(win_strt, win_end);
   return(0);
}

///////////////////////////////////////////////////////////////////////////
//////////////////    Some Currently Unused Code    ///////////////////////
///////////////////////////////////////////////////////////////////////////
// For every event sorted, the current user-sort re-calculates each condition
//    for every event in the window, which for an average of 10 events in window
//     => all conditions are calculated ~10 times PER EVENT
//  *** BUT there are typically zero user histos and conditions ***
//
// For the (never-realized) case of many user histograms, each with many
// conditions, a more efficient way of doing the user sort would be needed ...
// 
// The code below, allows the conditions to be tracked/updated each time events
// entered or left the main sort window [removing the window-size multiplier]
//  -> maintain current condn-passed-count and
//     on entry[exit] check if condition passed, and inc[dec] relevant counter
//////////////////////////////////////////////////////////////////////////////

// reset window counters, which will not be clear after a previous sort
//   (due to not fully clearing coincwin - may fix this later)
int user_sort_init()
{
return(0);
   Config *cfg = configs[1]; // ** Sort is using config[1]
   Cond *cond;
   int i;
   for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
      if( cond->use_count == 0 ){ continue; }
      cond->pass_count = 0;
   }
   return(0);
}

// one CORRECT *single-var* gate-checking method is (as in electronics)
//    keep running counts of either conditions or detector-types in window
//  - inc on adding, dec on leaving

// fragment new_frag has just entered coincwin (which starts at win_strt)
// update any condition counters here, for use later during sort
int user_addto_window(int win_strt, int new_frag)
{
return(0);
   Grif_event *ptr = &grif_event[new_frag];
   Config *cfg = configs[1]; // ** Sort is using config[1]
   int i, *val = (int *)ptr;
   Sortvar *var;
   Cond *cond;
   for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
      if( cond->use_count == 0 ){ continue; }
      var = cond->var;
      if( var->offset == -1 || var->dtype != ptr->dtype || var->local ){
         continue;
      }
      var->value = val[var->offset];
      switch( cond->op ){
      case GATEOP_LT: cond->pass_count += var->value <  cond->value; continue;
      case GATEOP_LE: cond->pass_count += var->value <= cond->value; continue;
      case GATEOP_GT: cond->pass_count += var->value >  cond->value; continue;
      case GATEOP_GE: cond->pass_count += var->value >= cond->value; continue;
      case GATEOP_EQ: cond->pass_count += var->value == cond->value; continue;
      default: printf("Unknown condition operation:%d\n", cond->op);
      }
   }
   return(0);
}
// fragment win_strt is just leaving coincwin
// update any condition counters
int user_removefrom_window(int win_strt, int new_frag)
{
return(0);
   Grif_event *ptr = &grif_event[win_strt];
   Config *cfg = configs[1]; // ** Sort is using config[1]
   int i, *val = (int *)ptr;
   Sortvar *var;
   Cond *cond;
   for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
      if( cond->use_count == 0 ){ continue; }
      var = cond->var;
      if( var->offset == -1 || var->dtype != ptr->dtype || var->local ){
         continue;
      }
      var->value = val[var->offset];
      switch( cond->op ){
      case GATEOP_LT: cond->pass_count -= var->value <  cond->value; continue;
      case GATEOP_LE: cond->pass_count -= var->value <= cond->value; continue;
      case GATEOP_GT: cond->pass_count -= var->value >  cond->value; continue;
      case GATEOP_GE: cond->pass_count -= var->value >= cond->value; continue;
      case GATEOP_EQ: cond->pass_count -= var->value == cond->value; continue;
      default: printf("Unknown condition operation:%d\n", cond->op);
      }
   }
   return(0);
}

int test_gate(Grif_event *ptr, Grif_event *alt, Gate *gate)
{
   Sortvar *var;
   float value;
   Cond *cond;
   int i, j;

   for(i=0; i<gate->nconds; i++){ cond = gate->conds[i];
      if( cond->valid == 0 ){ // not yet tested for this pair
         var = cond->var;
         value = var->get_value(ptr, alt, var->subsys1, var->subsys2);
         switch( cond->op ){
         case GATEOP_LT: cond->passed = value <  cond->value; break;
         case GATEOP_LE: cond->passed = value <= cond->value; break;
         case GATEOP_GT: cond->passed = value >  cond->value; break;
         case GATEOP_GE: cond->passed = value >= cond->value; break;
         case GATEOP_EQ: cond->passed = value == cond->value; break;
         }
         cond->valid = 1;
      }
      if( ! cond->passed ){ break; }
   }
   gate->valid  = ( i >= gate->nconds-1 ); // tested all conditions
   gate->passed = ( i == gate->nconds   );
   return(0);
}

//
const int ge_angles[64][64];
int  clover_angles[16][16];
int sceptar_angles[20][20];
int   gesep_angles[64][20];
int   gepac_angles[64][ 5];
int   gelbl_angles[64][20];
int  seplbl_angles[20][20];

float grif_crystal_theta[64]={
    36.5,  55.1,  55.1,  36.5,    36.5,  55.1,  55.1,  36.5,
    36.5,  55.1,  55.1,  36.5,    36.5,  55.1,  55.1,  36.5,
    80.6,  99.4,  99.4,  80.6,    80.6,  99.4,  99.4,  80.6,
    80.6,  99.4,  99.4,  80.6,    80.6,  99.4,  99.4,  80.6,
    80.6,  99.4,  99.4,  80.6,    80.6,  99.4,  99.4,  80.6,
    80.6,  99.4,  99.4,  80.6,    80.6,  99.4,  99.4,  80.6,
   124.9, 143.5, 143.5, 124.9,   124.9, 143.5, 143.5, 124.9,
   124.9, 143.5, 143.5, 124.9,   124.9, 143.5, 143.5, 124.9
};
float grif_crystal_phi[64]={
    83.4,  79.0,  56.0,  51.6,   173.4, 169.0, 146.0, 141.6,
   263.4, 259.0, 236.0, 231.6,   353.4, 349.0, 326.0, 321.6,
    32.0,  32.0,  13.0,  13.0,    77.0,  77.0,  58.0,  58.0,
   122.0, 122.0, 103.0, 103.0,   167.0, 167.0, 148.0, 148.0,
   212.0, 212.0, 193.0, 193.0,   257.0, 257.0, 238.0, 238.0,
   302.0, 302.0, 283.0, 283.0,   347.0, 347.0, 328.0, 328.0,
    79.0,  83.4,  51.6,  56.0,   169.0, 173.4, 141.6, 146.0,
   259.0, 263.4, 231.6, 236.0,   349.0, 353.4, 321.6, 326.0
};

/////////////////////////////////////////////////////////////////////////////
extern int crystal_table[MAX_DAQSIZE];
int calc_coincvars(Grif_event *ptr1, Grif_event *ptr2)
{
//////////////////////////////////////////////////////////////////////////////
// invalidate all variables, then fill in any values relevant to this event
//    ALSO SKIP ANY THAT ARE NOT IN USE (in any conditions/histograms)
//     - currently not much gain, unless can group by detector-type
//          [dtype->in_use]
   return(0);
}

//#######################################################################
//#####                     SORT VARIABLES                          #####
//#######################################################################
// to allow simple adding of more below, as required, web interface should request the list of variables

float uval_empty (Grif_event *ptr, Grif_event *alt, int s1, int s2){ printf("Unimplemented variable Value\n"); return( 0 ); }
float uval_ecal  (Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 ? ptr->ecal : -1 ); }
float uval_esum  (Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 ? ptr->esum : -1 ); }
float uval_cfd   (Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 ? ptr->cfd  : -1   ); }
float uval_ts    (Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 ? ptr->ts   : -1   ); }
float uval_apos  (Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 ? ptr->array_posn : -1   ); }
float uval_nhit  (Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 ? ptr->nhit : -1  ); }
float uval_ecal_s(Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 && ptr->suppress == 0 ? ptr->ecal : -1 ); }
float uval_esum_s(Grif_event *ptr, Grif_event *alt, int s1, int s2){ return( ptr->subsys == s1 && ptr->suppress == 0 ? ptr->esum : -1 ); }

float uval_ph    (Grif_event *ptr, Grif_event *alt, int s1, int s2){
   if( ptr->subsys != s1 ){ return(-1); }
   return( ( ptr->integ1 == 0 ) ? ptr->q1 : spread(ptr->q1)/ptr->integ1 );
}
float uval_xtl(Grif_event *ptr, Grif_event *alt, int s1, int s2){
   if( ptr->subsys != s1 || ptr->chan == -1 || ptr->chan >= MAX_DAQSIZE ){ return(-1); } 
   return( crystal_table[ptr->chan] );
}
float uval_dt(Grif_event *ptr, Grif_event *alt, int s1, int s2){
   return( ptr->subsys == s1 && alt->subsys == s2 ? ptr->cfd - alt->cfd : -1000000 );
}
float uval_dts(Grif_event *ptr, Grif_event *alt, int s1, int s2){
   return( ptr->subsys == s1 && alt->subsys == s2 ? ptr->ts - alt->ts : -1000000 );
}

typedef float(*grif_get_value)(Grif_event *, Grif_event *, int, int);
typedef float(*config_get_value)(void *, void *, int, int);

struct sortvar_desc_struct { // float32 good for integers to 16.8M
   //float (*get_value)(Grif_event *, Grif_event *, int, int);
   grif_get_value get_value;
   char name[32];  char desc[256];  int subsys1;  int subsys2;
};
#define NUM_SORTVARS 93
struct sortvar_desc_struct sortvarlist[NUM_SORTVARS]={
 // HPGE
   {&uval_ecal_s,"HPGeE",     "Singles HPGe Energy in keV with Compton suppression",               SUBSYS_HPGE_A, -1},
   {&uval_esum_s,"HPGeA",     "Addback HPGe Energy in keV with Compton suppression",               SUBSYS_HPGE_A, -1},
   {&uval_cfd,   "HPGeT",     "HPGe Time from CFD in 10/16 nanosecond steps",                      SUBSYS_HPGE_A, -1},
   {&uval_ts,    "HPGeTS",    "HPGe Timestamp value from leading-edge in 10 nanoseconds steps",    SUBSYS_HPGE_A, -1},
   {&uval_ph,    "HPGePH",    "Singles HPGe raw Pulse Height in ADC units",                        SUBSYS_HPGE_A, -1},
   {&uval_xtl,    "HPGeC",    "HPGe detector number (1-64)",                                       SUBSYS_HPGE_A, -1},
   {&uval_apos,  "HPGeCL",    "HPGe Clover number (1-16)",                                         SUBSYS_HPGE_A, -1},
   {&uval_nhit,  "HPGePU",    "HPGe Pileup value equal to number of Hits",                         SUBSYS_HPGE_A, -1},
   {&uval_empty, "HPGeIP",    "HPGe Integration period of the Pulse Height evaluation algorithum", SUBSYS_HPGE_A, -1},
   {&uval_empty, "HPGeDT",    "HPGe deadtime accumulated since previous accepted hit",             SUBSYS_HPGE_A, -1},
   {&uval_ecal,  "HPGeEU",    "Singles HPGe Energy in keV without Compton suppression",            SUBSYS_HPGE_A, -1},
   {&uval_esum,  "HPGeAU",    "Addback HPGe Energy in keV without Compton suppression",            SUBSYS_HPGE_A, -1},
   {&uval_empty, "GRGTHETA",  "Theta angle of HPGE crystal with respect to the beam axis",         SUBSYS_HPGE_A, -1},
   {&uval_empty, "GRGPHI",    "Phi angle of HPGE crystal with respect to lab coordinate system",   SUBSYS_HPGE_A, -1},
   {&uval_empty, "CLVTHETA",  "Theta angle of centre of HPGE clover with respect to the beam axis",SUBSYS_HPGE_A, -1},
   {&uval_empty, "CLVPHI",    "Phi angle of centre of HPGE clover wrt lab coordinate system",      SUBSYS_HPGE_A, -1},
 // SCEPTAR
   {&uval_ecal,  "SEPE",     "SCEPTAR Pulse Height in arbitrary units",                            SUBSYS_SCEPTAR, -1},
   {&uval_cfd,   "SEPT",     "SCEPTAR Time from CFD in nanoseconds",                               SUBSYS_SCEPTAR, -1},
   {&uval_ts,    "SEPTS",    "SCEPTAR Timestamp value from leading-edge in 10 nanoseconds steps",  SUBSYS_SCEPTAR, -1},
   {&uval_ph,    "SEPPH",    "SCEPTAR raw Pulse Height in ADC units",                              SUBSYS_SCEPTAR, -1},
   {&uval_nhit,  "SEPPU",    "SCEPTAR Pileup value equal to number of Hits",                       SUBSYS_SCEPTAR, -1},
   {&uval_empty, "SEPTHETA", "Theta angle of SCEPTAR paddle with respect to the beam axis",        SUBSYS_SCEPTAR, -1},
   {&uval_empty, "SEPPHI",   "Phi angle of SCEPTAR paddle with respect to lab coordinate system",  SUBSYS_SCEPTAR, -1},
   {&uval_xtl,   "SEPNUM",   "The number of this SCEPTAR paddle [1-20]",                           SUBSYS_SCEPTAR, -1},
 // PACES
   {&uval_ecal,  "PACE",    "PACES energy in keV",                                                 SUBSYS_PACES, -1},
   {&uval_cfd,   "PACT",    "PACES Time from CFD in nanoseconds",                                  SUBSYS_PACES, -1},
   {&uval_ts,    "PACTS",   "PACES Timestamp value from leading-edge in 10 nanoseconds steps",     SUBSYS_PACES, -1},
   {&uval_ph,    "PACPH",   "PACES raw Pulse Height in ADC units",                                 SUBSYS_PACES, -1},
   {&uval_nhit,  "PACPU",   "PACES Pileup value equal to number of Hits",                          SUBSYS_PACES, -1},
   {&uval_empty, "PACTHETA","Theta angle of PACES crystal with respect to the beam axis",          SUBSYS_PACES, -1},
   {&uval_empty, "PACPHI",  "Phi angle of PACES crystal with respect to lab coordinate system",    SUBSYS_PACES, -1},
   {&uval_xtl,   "PACNUM",  "The number of this PACES crystal [1-5]",                              SUBSYS_PACES, -1},
 // LaBr3
   {&uval_ecal,  "LBLE",    "LaBr3 energy in keV",                                                 SUBSYS_LABR_L, -1},
   {&uval_cfd,   "LBLT",    "LaBr3 Time from CFD in nanoseconds",                                  SUBSYS_LABR_L, -1},
   {&uval_ts,    "LBLTS",   "LaBr3 Timestamp value from leading-edge in 10 nanoseconds steps",     SUBSYS_LABR_L, -1},
   {&uval_ph,    "LBLPH",   "LaBr3 raw Pulse Height in ADC units",                                 SUBSYS_LABR_L, -1},
   {&uval_nhit,  "LBLPU",   "LaBr3 Pileup value equal to number of Hits",                          SUBSYS_LABR_L, -1},
   {&uval_empty, "LBLTHETA","Theta angle of LaBr3 crystal with respect to the beam axis",          SUBSYS_LABR_L, -1},
   {&uval_empty, "LBLPHI",  "Phi angle of LaBr3 crystal with respect to lab coordinate system",    SUBSYS_LABR_L, -1},
   {&uval_xtl,   "LBLNUM",  "The number of this LaBr3 crystal [1-8]",                              SUBSYS_LABR_L, -1},
// TACs
   {&uval_ecal,  "LBTE",   "TAC Pulse Height in arbitrary units",                                  SUBSYS_TAC_LABR, -1},
   {&uval_cfd,   "LBTT",   "TAC Time from CFD in nanoseconds",                                     SUBSYS_TAC_LABR, -1},
   {&uval_ts,    "LBTTS",  "TAC Timestamp value from leading-edge in 10 nanoseconds steps",        SUBSYS_TAC_LABR, -1},
   {&uval_ph,    "LBTPH",  "TAC raw Pulse Height in ADC units",                                    SUBSYS_TAC_LABR, -1},
   {&uval_nhit,  "LBTPU",  "TAC Pileup value equal to number of Hits",                             SUBSYS_TAC_LABR, -1},
   {&uval_xtl,   "LBTNUM", "The number of this TAC module [1-8]",                                  SUBSYS_TAC_LABR, -1},
// Clover BGO Suppression shields
   {&uval_ecal,  "GRSE",    "Clover BGO energy in keV",                                            SUBSYS_BGO, -1},
   {&uval_cfd,   "GRST",    "Clover BGO Time from CFD in nanoseconds",                             SUBSYS_BGO, -1},
   {&uval_ts,    "GRSTS",   "Clover BGO Timestamp value from leading-edge in 10 nanoseconds steps",SUBSYS_BGO, -1},
   {&uval_ph,    "GRSPH",   "Clover BGO raw Pulse Height in ADC units",                            SUBSYS_BGO, -1},
   {&uval_nhit,  "GRSPU",   "Clover BGO Pileup value equal to number of Hits",                     SUBSYS_BGO, -1},
   {&uval_xtl,   "GRSNUM",  "The number of this Clover BGO crystal [1 to 20*16=320]",              SUBSYS_BGO, -1},
   {&uval_empty, "GRSPOS",  "The HPGe Clover number (1-16) to which this BGO belongs",             SUBSYS_BGO, -1},
   {&uval_empty, "GRSTYPE", "The type of this Clover BGO crystal [front, side, back]",             SUBSYS_BGO, -1},
// Ancillary position BGO Suppression shield,
   {&uval_ecal,  "LBSE",    "Ancillary position BGO energy in keV",                                   SUBSYS_LABR_BGO, -1},
   {&uval_cfd,   "LBST",    "Ancillary position BGO Time from CFD in nanoseconds",                    SUBSYS_LABR_BGO, -1},
   {&uval_ts,    "LBSTS",   "Ancillary position BGO Timestamp value from leading-edge in 10 ns steps",SUBSYS_LABR_BGO, -1},
   {&uval_ph,    "LBSPH",   "Ancillary position BGO raw Pulse Height in ADC units",                   SUBSYS_LABR_BGO, -1},
   {&uval_nhit,  "LBSPU",   "Ancillary position BGO Pileup value equal to number of Hits",            SUBSYS_LABR_BGO, -1},
   {&uval_xtl,   "LBSNUM",  "The number of this Ancillary position BGO crystal [1 to 3*8=24]",        SUBSYS_LABR_BGO, -1},
   {&uval_empty, "LBSPOS",  "The ancillary position number (1-8) to which this BGO belongs",          SUBSYS_LABR_BGO, -1},
// Time Differences","\t    {
   {&uval_empty,"MIDAS_Time",  "Time since the beginning of run based on the MIDAS CPU time",                     -1, -1},
   {&uval_dt,   "TD_GRG_GRG",  "Time difference between two HPGe crystals using CFD in nanoseconds",           SUBSYS_HPGE_A, SUBSYS_HPGE_A},
   {&uval_dt,   "TD_SEP_SEP",  "Time difference between two SCEPTAR paddles using CFD in nanoseconds",         SUBSYS_SCEPTAR,SUBSYS_SCEPTAR},
   {&uval_dt,   "TD_PAC_PAC",  "Time difference between two PACES crystals using CFD in nanoseconds",          SUBSYS_PACES,  SUBSYS_PACES},
   {&uval_dt,   "TD_LBL_LBL",  "Time difference between two LaBr3 detectors using CFD in nanoseconds",         SUBSYS_LABR_L, SUBSYS_LABR_L},
   {&uval_dts,  "TD_GRG_GRG",  "Time difference between two HPGe crystals using leading edge in 10 ns units",  SUBSYS_HPGE_A, SUBSYS_HPGE_A},
   {&uval_dts,  "TD_SEP_SEP",  "Time difference between two SCEPTAR paddles using leading edge in 10 ns units",SUBSYS_SCEPTAR,SUBSYS_SCEPTAR},
   {&uval_dts,  "TD_PAC_PAC",  "Time difference between two PACES crystals using leading edge in 10 ns units", SUBSYS_PACES,  SUBSYS_PACES},
   {&uval_dts,  "TD_LBL_LBL",  "Time difference between two LaBr3 detectors using leading edge in 10 ns units",SUBSYS_LABR_L, SUBSYS_LABR_L},
   {&uval_dt,   "TD_GRG_SEP",  "Time difference between HPGe and SCEPTAR using CFD in nanoseconds",            SUBSYS_HPGE_A, SUBSYS_SCEPTAR},
   {&uval_dts,  "TSD_GRG_SEP", "Time difference between HPGe and SCEPTAR using leading edge in 10 ns units",   SUBSYS_HPGE_A, SUBSYS_SCEPTAR},
   {&uval_dt,   "TD_GRG_PAC",  "Time difference between HPGe and PACES using CFD in nanoseconds",              SUBSYS_HPGE_A, SUBSYS_PACES},
   {&uval_dts,  "TSD_GRG_PAC", "Time difference between HPGe and PACES using leading edge in 10 ns units",     SUBSYS_HPGE_A, SUBSYS_PACES},
   {&uval_dt,   "TD_SEP_PAC",  "Time difference between PACES and SCEPTAR using CFD in nanoseconds",           SUBSYS_PACES,  SUBSYS_SCEPTAR},
   {&uval_dts,  "TSD_SEP_PAC", "Time difference between PACES and SCEPTAR using leading edge in 10 ns units",  SUBSYS_PACES,  SUBSYS_SCEPTAR},
   {&uval_dt,   "TD_GRG_LBL",  "Time difference between HPGe and LaBr3 using CFD in nanoseconds",              SUBSYS_HPGE_A, SUBSYS_LABR_L},
   {&uval_dts,  "TSD_GRG_LBL", "Time difference between HPGe and LaBr3 using leading edge in 10 ns units",     SUBSYS_HPGE_A, SUBSYS_LABR_L},
   {&uval_dt,   "TD_SEP_LBL",  "Time difference between LaBr3 and SCEPTAR using CFD in nanoseconds",           SUBSYS_LABR_L, SUBSYS_SCEPTAR},
   {&uval_dts,  "TSD_SEP_LBL", "Time difference between LaBr3 and SCEPTAR using leading edge in 10 ns units",  SUBSYS_LABR_L, SUBSYS_SCEPTAR},
// Angular Differences - for HPGe need to know the distance of 110mm or 145mm
   {&uval_empty, "ANG_GRG_GRG",    "Angular difference between the centre of two HPGe crystals in degrees",                  -1, 0},
   {&uval_empty, "ANG_CLV_CLV",    "Angular difference between the centre of two HPGe clovers in degrees",                   -1, 0},
   {&uval_empty, "ANG_SEP_SEP",    "Angular difference between the centre of two SCEPTAR paddles in degrees",                -1, 0},
   {&uval_empty, "ANG_PAC_PAC",    "Angular difference between the centre of two PACES crystals in degrees",                 -1, 0},
   {&uval_empty, "ANG_LBL_LBL",    "Angular difference between the centre of two LaBr3 detectors in degrees",                -1, 0},
   {&uval_empty, "ANG_GRG_SEP",    "Angular difference between the centre of a HPGe crystal and SCEPTAR paddle in degrees",  -1, 0},
   {&uval_empty, "ANG_GRG_PAC",    "Angular difference between the centre of a HPGe crystal and PACES crystal in degrees",   -1, 0},
   {&uval_empty, "ANG_GRG_LBL",    "Angular difference between the centre of a HPGe crystal and LaBr3 detector in degrees",  -1, 0},
   {&uval_empty, "ANG_PAC_LBL",    "Angular difference between the centre of a PACES crystal and LaBr3 detector in degrees", -1, 0},
   {&uval_empty, "ANG_SEP_LBL",    "Angular difference between the centre of a LaBr3 detector and SCEPTAR paddle in degrees",-1, 0},
// Cycle Timing (PPG events)","\t    {
   {&uval_empty, "Cycle_Num",    "Cycle number since beginning of run",                                                      -1, 0},
   {&uval_empty, "Cycle_Time",    "Time since the beginning of the current Cycle",                                           -1, 0},
   {&uval_empty, "Cycle_Pattern",    "The current PPG pattern",                                                              -1, 0}
};

int init_user_config(Config *cfg) // copy sort variable structure to config
{                           // (have multiple versions of these with varying "in-use" etc)
   Sortvar *s, *ptr;
   char *tmp;
   int i;

   for(i=0; i<NUM_SORTVARS; i++){ ++cfg->nsortvar;
      tmp = sortvarlist[i].name;
      if( strlen(tmp) == 0 ){ continue; }
      memcpy(cfg->varlist[i].name, tmp, strlen(tmp)+1);
      tmp = sortvarlist[i].desc;
      memcpy(cfg->varlist[i].title, tmp, strlen(tmp)+1);
      cfg->varlist[i].subsys1 = sortvarlist[i].subsys1;
      cfg->varlist[i].subsys2 = sortvarlist[i].subsys2;
      cfg->varlist[i].get_value  = (config_get_value)sortvarlist[i].get_value;

      // TEMP to avoid changing config.c too much before merge
      //   reuse old sortvar members for new items
      //cfg->varlist[i].
   }
   return(0);
}
