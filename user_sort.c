#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grif-format.h"
#include "config.h"
#include "histogram.h"

//#######################################################################
//#####              USER SORT - custom user histograms             #####
//#######################################################################

// Both default and user sorts need to loop over coinc-fragments
//    very little gain to having loop external to both
//    also user_sort needs to keep track of histo increments per event
// - just run loop twice (once in default and once in user)
//
// singles user_sort - varlist.offset != -1 => singles variable
//


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

// fragment at win_start is leaving coincwin (frag[win_end+1] is not in coinc)
extern Grif_event grif_event[PTR_BUFSIZE];
int user_sort(int win_strt, int win_end, int flag)
{
return(0);
   Grif_event *alt, *ptr = &grif_event[win_strt];
   Config *cfg = configs[1]; // ** Sort is using config[1]
   int i, j, k, m, yvalid, yval, dt, *val = (int *)ptr, len = win_end-win_strt;
   Sortvar *var, *yvar;
   Histogram *histo;
   Gate *gate;
   Cond *cond;

   if( len < 0 ){ len += PTR_BUFSIZE; }

   // ------------ histogram increment tracking ----------------------
   // only want to increment histograms once per set of coinc-frags
   // (unless histo is 2d with DIFF vars - does not include e.g. EvsXtal)
   // ( also - Need to specify which Xtal e.g. XtalGamma1 or XtalGamma2)
   //          ( can ignore above for now as may never use that)
   //
   // so running though set of coinc-pairs and incrementing histos
   //   most only want to inc ONCE if condition is ever met
   // NEED TO KEEP TRACK HERE IF HISTO DONE FOR FRAG_IDX
   //   initialize "done_flag" in first singles run through
   //   set if histo is ever incremented after this point
   //for(i=0; i<cfg->nsortvar; i++){ var = &cfg->varlist[i];
   for(i=0; i<cfg->nusedvar; i++){ var = cfg->usedvars[i];
      // set valid based on current fragment
      //  - used for all 1d-histos, and for x-axis of 2d
      // (if valid, will loop 2d histos over other frags)
      // NOTE coincvars with offset=-1 will be recalculated below
      if( var->offset == -1 ){ continue; }
      var->valid = var->dtype == ptr->dtype;
      var->value = val[var->offset];
      for(j=0; j<var->use_count_x; j++){
         var->histo_list_x[j]->done_flag = 0;
      }
   }
   //for(i=0; i<cfg->nuser; i++){ cfg->user_histos[i]->done_flag = 0; }

   // ------------ "single-variable" conditions ----------------------
   //       check  once here, not per-coinc-pair below
   for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
      if( cond->use_count == 0 ){ continue; }
      var = cond->var; if( var->offset == -1 ){ continue; }
      if( var->local ){      // NOT a window-wide variable e.g. crystal-number
         var->value = val[var->offset];  // => check against current frag ONLY
         switch( cond->op ){
         case GATEOP_LT: cond->passed = var->value <  cond->value; break;
         case GATEOP_LE: cond->passed = var->value <= cond->value; break;
         case GATEOP_GT: cond->passed = var->value >  cond->value; break;
         case GATEOP_GE: cond->passed = var->value >= cond->value; break;
         case GATEOP_EQ: cond->passed = var->value == cond->value; break;
         default: printf("Unknown condition operation:%d\n", cond->op);
         }
      } else {
         cond->passed = (cond->pass_count > 0);
      }
   }
   // fill histograms ...
   //    there are thousands of histos (not so many user-histos)
   //    there are ~100variables
   //    Each histogram has only a single x-variable
   //     => double loop over vars+dep-histos includes all histograms once
   //        (as long as dep-histo Xvar is this variable (not Yvar))
   if( ptr->chan != -1 ){
      m = win_strt;
      while( 1 ){
         if( (( win_end < m ) ? win_end + PTR_BUFSIZE - m :
                win_end - m ) > len ){
            printf("CORRUPT\n");
         }
         if( ++m == PTR_BUFSIZE ){ m = 0; } // wrap
         alt = &grif_event[m];
         if( alt->chan != -1 && m <= win_end ){
            calc_coincvars(ptr, alt); test_gates(ptr, alt);
         }
         //for(i=0; i<cfg->nsortvar; i++){ var = &cfg->varlist[i];
         ////for(i=0; i<cfg->nusedvar; i++){ var = cfg->usedvars[i];
         //   if( var->valid == 0 ){ continue; }
         //   for(j=0; j<var->use_count_x; j++){ histo = var->histo_list_x[j];
         for(i=0; i<cfg->nuser; i++){ histo = cfg->user_histos[i];
               if( histo->done_flag        ){ continue; }
               var = histo->xvar;
             //if( histo->xvar->valid == 0 ){ continue; }
               if( var->dtype != ptr->dtype ){ continue; }
               for(k=0; k<histo->num_gates; k++){
                  if( histo->gate_passed[k] == 0 ){ break; }
               }
               if( k < histo->num_gates ){ continue; }
               if(        histo->type == INT_1D ){
                  histo->Fill(histo, var->value, 1);
                  histo->done_flag = 1;
               } else if( histo->type == INT_2D &&
                  alt->chan != -1 && m <= win_end ){ yvar = histo->yvar;
                  // 2d histos - already checked all gates + set valids
                  // Do not fill if yaxis variable was not seen
                  if( yvar->offset != -1 ){
                     yvalid = yvar->dtype == alt->dtype;
                     yval = *(((int *)alt)+yvar->offset);
                  } else {
                     yvalid = yvar->valid; // set in calc_coincvars
                     yval = yvar->value;
                  }
                  if( !yvalid ){ continue; }
                  ((TH2I *)histo)->Fill((TH2I *)histo, var->value, yval, 1);
                  // ############# do not set done_flag for true 2d histos
               }
         }
         i = m - win_strt; if( i < 0 ){ i += PTR_BUFSIZE; }
         if( i >= len ){ break; }
      }
   }
   user_removefrom_window(win_strt, win_end);
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

int test_gates(Grif_event *ptr, Grif_event *alt)
{
return(0);
   Config *cfg = configs[1]; // ** Sort is using config[1]
   int i, j, value;
   Sortvar *var;
   Gate *gate;
   Cond *cond;
   for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
      if( cond->use_count == 0 ){ continue; }
      var = cond->var; if( var->offset != -1 ){ continue; } // already done
      cond->passed = 0; // value+valid already set in calc_coincvars
      if( var->valid == 0 ){ continue; }

      switch( cond->op ){
      case GATEOP_LT: cond->passed = var->value <  cond->value; continue;
      case GATEOP_LE: cond->passed = var->value <= cond->value; continue;
      case GATEOP_GT: cond->passed = var->value >  cond->value; continue;
      case GATEOP_GE: cond->passed = var->value >= cond->value; continue;
      case GATEOP_EQ: cond->passed = var->value == cond->value; continue;
      }
   }
   for(i=0; i<cfg->ngates; i++){
      gate = cfg->gatelist[i]; gate->passed = 1;
      for(j=0; j<gate->nconds; j++){ cond = gate->conds[j];
         if( cond->passed == 0 ){ gate->passed = 0; break; }
      }
   }
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
// array of varlist structures is in config
// BUT actual values are scattered in multiple Grif_Evt structures
// => to get values, config structure must contain offsets into Grif_Evt
//    e.g. offset=10, value = *( ((int *)ptr) + offset )
//    int *val = (int *)ptr;   value = val[offset];
//
// => Do not need valid - just check dtype
//
// coinc stuff e.g. time-diff, IS NOT IN Grif_evt
// set offset to be -1 => value is in varlist structure
// (valid is required here - will be set when value calculated from evt-pair)
//
// could also store type (which would allow floating point etc.)
// [and can then have both int and float values within varlist-struct as well]
//////////////////////////////////////////////////////////////////////////////
// invalidate all variables, then fill in any values relevant to this event
//    ALSO SKIP ANY THAT ARE NOT IN USE (in any conditions/histograms)
//     - currently not much gain, unless can group by detector-type
//          [dtype->in_use]
/////////////////////////////////////////////////////////////////////////////
extern int crystal_table[MAX_DAQSIZE];
int calc_coincvars(Grif_event *ptr1, Grif_event *ptr2)
{
   int gege, gesep, gepac, gelbs, sepsep, seppac, seplbs, pacpac, lbslbs;
   Sortvar *s, *ts, *tsd, *ang;
   int i, p1, p2, c1, c2, swap;
   Config *cfg = configs[1];
   Grif_event *tmp;

   return(0);
   for(i=0; i<cfg->nsortvar; i++){
      if( cfg->varlist[i].offset != -1 ){ continue; } // not coinc variable
      cfg->varlist[i].valid = 0;                      // invalidate
   }
   gege = gesep = gepac = gelbs = 0;
   sepsep = seppac = seplbs = pacpac = lbslbs = 0;
   // time differences, and angle between detectors ...
   swap=0;
   while( 1 ){
      ts = tsd = NULL;
      p1 = ptr1->array_posn;           p2 = ptr2->array_posn;
      c1 = crystal_table[ptr1->chan];  c2 = crystal_table[ptr2->chan];
      if(        ptr1->dtype == 0 ){ // Ge
         if(        ptr2->dtype == 0 ){ gege = 1; // Ge-Ge
            ts  = &cfg->varlist[TD_GRG_GRG ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_GRG_GRG]; tsd->valid = 1;
            ang = &cfg->varlist[ANG_CLV_CLV]; ang->valid = 1;
            ang->value = clover_angles[p1][p2];
            ang = &cfg->varlist[ANG_GRG_GRG]; ang->valid = 1;
            ang->value = ge_angles[c1][c2];
         } else if( ptr2->dtype == 2 ){ gesep = 1; // Ge-Sep
            ts  = &cfg->varlist[TD_GRG_SEP ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_GRG_SEP]; tsd->valid = 1;
            ang = &cfg->varlist[ANG_GRG_SEP]; ang->valid = 1;
            ang->value = gesep_angles[c1][p2];
         } else if( ptr2->dtype == 5 ){ gepac = 1; // Ge-PAC
            ts  = &cfg->varlist[TD_GRG_PAC ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_GRG_PAC]; tsd->valid = 1;
            ang = &cfg->varlist[ANG_GRG_PAC]; ang->valid = 1;
            ang->value = gepac_angles[c1][p2];
         } else if( ptr2->dtype == 3 ){ gelbs = 1; // Ge-LaBr
            ts  = &cfg->varlist[TD_GRG_LBL ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_GRG_LBL]; tsd->valid = 1;
            ang = &cfg->varlist[ANG_GRG_LBL]; ang->valid = 1;
            ang->value = gelbl_angles[c1][p2];
         }
      } else if( ptr1->dtype == 2 ){ // Sep
         if(        ptr2->dtype == 2 ){ sepsep = 1; // Sep-Sep
            ts  = &cfg->varlist[TD_SEP_SEP ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_SEP_SEP]; tsd->valid = 1;
            ang = &cfg->varlist[ANG_SEP_SEP]; ang->valid = 1;
            ang->value = sceptar_angles[p1][p2];
         } else if( ptr2->dtype == 5 ){ seppac = 1; // Sep-PAC
            ts  = &cfg->varlist[TD_SEP_PAC ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_SEP_PAC]; tsd->valid = 1;
         } else if( ptr2->dtype == 3 ){ seplbs = 1; // Sep-LaBr
            ts  = &cfg->varlist[TD_SEP_LBL ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_SEP_LBL]; tsd->valid = 1;
            ang = &cfg->varlist[ANG_SEP_LBL]; ang->valid = 1;
            ang->value = seplbl_angles[p1][p2];
         }
      } else if( ptr1->dtype == 5 ){ // Pac
         if(        ptr2->dtype == 5 ){ pacpac = 1; // Pac-Pac
            ts  = &cfg->varlist[TD_PAC_PAC ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_PAC_PAC]; tsd->valid = 1;
         }
      } else if( ptr1->dtype == 3 ){ // Labr
         if(        ptr2->dtype == 3 ){ lbslbs = 1; // Labr-Labr
            ts  = &cfg->varlist[TD_LBL_LBL ]; ts ->valid = 1;
            tsd = &cfg->varlist[TSD_LBL_LBL]; tsd->valid = 1;
         }
      }
      if( ts != NULL && tsd != NULL ){
         ts->value = ptr1->cfd - ptr2->cfd;
         if(ts->value<0){ts ->value *= -1;}
         tsd->value = ptr1->ts  - ptr2->ts;
         if(tsd->value<0){tsd->value *= -1;}
      }
      if( swap || ptr1->dtype == ptr2->dtype ){ break; }
      // repeat with events switched - for reverse combination
      swap = 1; tmp = ptr1; ptr1 = ptr2; ptr2 = tmp;
   }
   // Other coinc stuff ...
   // if(        gege   ){ // Ge-Ge
   // } else if( gesep  ){ // Ge-Sep
   // } else if( gepac  ){ // Ge-PAC
   // } else if( gelbs  ){ // Ge-LaBr
   // } else if( sepsep ){ // Sep-Sep
   // } else if( seppac ){ // Sep-PAC
   // } else if( seplbs ){ // Sep-LaBr
   // } else if( pacpac ){ // Pac-Pac
   // } else if( lbslbs ){ // Labr-Labr
   // }
   return(0);
}

//#######################################################################
//#####                     SORT VARIABLES                          #####
//#######################################################################

#define NUM_SORTVARS 93
char *sortvarlist[2*NUM_SORTVARS]={
 // HPGE
   "HPGeE",    "Singles HPGe Energy in keV with Compton suppression",       // 0
   "HPGeA",    "Addback HPGe Energy in keV with Compton suppression",
   "HPGeT",    "HPGe Time from CFD in nanoseconds",
   "HPGeTS",    "HPGe Timestamp value from leading-edge in 10 nanoseconds steps",
   "HPGePH",    "Singles HPGe raw Pulse Height in ADC units",
   "HPGeC",    "HPGe detector number (1-64)",
   "HPGeCL",    "HPGe Clover number (1-16)",
   "HPGePU",    "HPGe Pileup value equal to number of Hits",
   "HPGeIP",    "HPGe Integration period of the Pulse Height evaluation algorithum",
   "HPGeDT",    "HPGe deadtime accumulated since previous accepted hit",
   "HPGeEU",    "Singles HPGe Energy in keV without Compton suppression",  // 10
   "HPGeAU",    "Addback HPGe Energy in keV without Compton suppression",
   "GRGTHETA",    "Theta angle of HPGE crystal with respect to the beam axis",
   "GRGPHI",    "Phi angle of HPGE crystal with respect to lab coordinate system",
   "CLVTHETA",    "Theta angle of centre of HPGE clover with respect to the beam axis",
   "CLVPHI",    "Phi angle of centre of HPGE clover with respect to lab coordinate system",
 // SCEPTAR"
   "SEPE",    "SCEPTAR Pulse Height in arbitrary units",
   "SEPT",    "SCEPTAR Time from CFD in nanoseconds",
   "SEPTS",    "SCEPTAR Timestamp value from leading-edge in 10 nanoseconds steps",
   "SEPPH",    "SCEPTAR raw Pulse Height in ADC units",
   "SEPPU",    "SCEPTAR Pileup value equal to number of Hits",  // 20
   "SEPTHETA",    "Theta angle of SCEPTAR paddle with respect to the beam axis",
   "SEPPHI",    "Phi angle of SCEPTAR paddle with respect to lab coordinate system",
   "SEPNUM",    "The number of this SCEPTAR paddle [1-20]",
 // PACES
   "PACE",    "PACES energy in keV",
   "PACT",    "PACES Time from CFD in nanoseconds",
   "PACTS",    "PACES Timestamp value from leading-edge in 10 nanoseconds steps",
   "PACPH",    "PACES raw Pulse Height in ADC units",
   "PACPU",    "PACES Pileup value equal to number of Hits",
   "PACTHETA",    "Theta angle of PACES crystal with respect to the beam axis",
   "PACPHI",    "Phi angle of PACES crystal with respect to lab coordinate system",   // 30
   "PACNUM",    "The number of this PACES crystal [1-5]",
 // LaBr3
   "LBLE",    "LaBr3 energy in keV",
   "LBLT",    "LaBr3 Time from CFD in nanoseconds",
   "LBLTS",    "LaBr3 Timestamp value from leading-edge in 10 nanoseconds steps",
   "LBLPH",    "LaBr3 raw Pulse Height in ADC units",
   "LBLPU",    "LaBr3 Pileup value equal to number of Hits",
   "LBLTHETA",    "Theta angle of LaBr3 crystal with respect to the beam axis",
   "LBLPHI",    "Phi angle of LaBr3 crystal with respect to lab coordinate system",
   "LBLNUM",    "The number of this LaBr3 crystal [1-8]",
// TACs
   "LBTE",    "TAC Pulse Height in arbitrary units",  // 40
   "LBTT",    "TAC Time from CFD in nanoseconds",
   "LBTTS",    "TAC Timestamp value from leading-edge in 10 nanoseconds steps",
   "LBTPH",    "TAC raw Pulse Height in ADC units",
   "LBTPU",    "TAC Pileup value equal to number of Hits",
   "LBTNUM",    "The number of this TAC module [1-8]",
// Clover BGO Suppression shields
   "GRSE",    "Clover BGO energy in keV",
   "GRST",    "Clover BGO Time from CFD in nanoseconds",
   "GRSTS",    "Clover BGO Timestamp value from leading-edge in 10 nanoseconds steps",
   "GRSPH",    "Clover BGO raw Pulse Height in ADC units",
   "GRSPU",    "Clover BGO Pileup value equal to number of Hits",  // 50
   "GRSNUM",    "The number of this Clover BGO crystal [1 to 20*16=320]",
   "GRSPOS",    "The HPGe Clover number (1-16) to which this BGO belongs",
   "GRSTYPE",    "The type of this Clover BGO crystal [front, side, back]",
// Ancillary position BGO Suppression shield,
   "LBSE",    "Ancillary position BGO energy in keV",
   "LBST",    "Ancillary position BGO Time from CFD in nanoseconds",
   "LBSTS",    "Ancillary position BGO Timestamp value from leading-edge in 10 nanoseconds steps",
   "LBSPH",    "Ancillary position BGO raw Pulse Height in ADC units",
   "LBSPU",    "Ancillary position BGO Pileup value equal to number of Hits",
   "LBSNUM",    "The number of this Ancillary position BGO crystal [1 to 3*8=24]",
   "LBSPOS",    "The ancillary position number (1-8) to which this BGO belongs",   // 60
// Time Differences","\t    {
   "MIDAS_Time",    "Time since the beginning of run based on the MIDAS CPU time",
   "TD_GRG_GRG",    "Time difference between two HPGe crystals using CFD in nanoseconds",
   "TD_SEP_SEP",    "Time difference between two SCEPTAR paddles using CFD in nanoseconds",
   "TD_PAC_PAC",    "Time difference between two PACES crystals using CFD in nanoseconds",
   "TD_LBL_LBL",    "Time difference between two LaBr3 detectors using CFD in nanoseconds",
   "TD_GRG_GRG",    "Time difference between two HPGe crystals using leading edge in 10 nanosecond units",
   "TD_SEP_SEP",    "Time difference between two SCEPTAR paddles using leading edge in 10 nanosecond units",
   "TD_PAC_PAC",    "Time difference between two PACES crystals using leading edge in 10 nanosecond units",
   "TD_LBL_LBL",    "Time difference between two LaBr3 detectors using leading edge in 10 nanosecond units",
   "TD_GRG_SEP",    "Time difference between HPGe and SCEPTAR using CFD in nanoseconds",    // 70
   "TSD_GRG_SEP",    "Time difference between HPGe and SCEPTAR using leading edge in 10 nanosecond units",
   "TD_GRG_PAC",    "Time difference between HPGe and PACES using CFD in nanoseconds",
   "TSD_GRG_PAC",    "Time difference between HPGe and PACES using leading edge in 10 nanosecond units",
   "TD_SEP_PAC",    "Time difference between PACES and SCEPTAR using CFD in nanoseconds",
   "TSD_SEP_PAC",    "Time difference between PACES and SCEPTAR using leading edge in 10 nanosecond units",
   "TD_GRG_LBL",    "Time difference between HPGe and LaBr3 using CFD in nanoseconds",
   "TSD_GRG_LBL",    "Time difference between HPGe and LaBr3 using leading edge in 10 nanosecond units",
   "TD_SEP_LBL",    "Time difference between LaBr3 and SCEPTAR using CFD in nanoseconds",
   "TSD_SEP_LBL",    "Time difference between LaBr3 and SCEPTAR using leading edge in 10 nanosecond units",
// Angular Differences - for HPGe need to know the distance of 110mm or 145mm
   "ANG_GRG_GRG",    "Angular difference between the centre of two HPGe crystals in degrees",      // 80
   "ANG_CLV_CLV",    "Angular difference between the centre of two HPGe clovers in degrees",
   "ANG_SEP_SEP",    "Angular difference between the centre of two SCEPTAR paddles in degrees",
   "ANG_PAC_PAC",    "Angular difference between the centre of two PACES crystals in degrees",
   "ANG_LBL_LBL",    "Angular difference between the centre of two LaBr3 detectors in degrees",
   "ANG_GRG_SEP",    "Angular difference between the centre of a HPGe crystal and SCEPTAR paddle in degrees",
   "ANG_GRG_PAC",    "Angular difference between the centre of a HPGe crystal and PACES crystal in degrees",
   "ANG_GRG_LBL",    "Angular difference between the centre of a HPGe crystal and LaBr3 detector in degrees",
   "ANG_PAC_LBL",    "Angular difference between the centre of a PACES crystal and LaBr3 detector in degrees",
   "ANG_SEP_LBL",    "Angular difference between the centre of a LaBr3 detector and SCEPTAR paddle in degrees",
// Cycle Timing (PPG events)","\t    {
   "Cycle_Num",    "Cycle number since beginning of run",                               // 90
   "Cycle_Time",    "Time since the beginning of the current Cycle",
   "Cycle_Pattern",    "The current PPG pattern",
};

int init_user_config(Config *cfg) // setup sort variable structures
{
   Sortvar *s, *ptr;
   char *tmp;
   int i;

   return(0);
   for(i=0; i<NUM_SORTVARS; i++){ ++cfg->nsortvar;
      tmp = sortvarlist[2*i];
      memcpy(cfg->varlist[i].name, tmp, strlen(tmp)+1);
      tmp = sortvarlist[2*i+1];
      memcpy(cfg->varlist[i].title, tmp, strlen(tmp)+1);
   }
   // the following set the location of the value in the Event structure
   // also the required datatype to be valid
   s=&cfg->varlist[HPGeE       ]; s->offset = 39; s->dtype = 0;
   s=&cfg->varlist[HPGeA       ]; s->offset = 38; s->dtype = 0;
   s=&cfg->varlist[HPGeT       ]; s->offset = 24; s->dtype = 0;
   s=&cfg->varlist[HPGeTS      ]; s->offset = 34; s->dtype = 0;
   s=&cfg->varlist[HPGePH      ]; s->offset =  8; s->dtype = 0;
   s=&cfg->varlist[HPGeC       ]; s->offset = 40; s->dtype = 0; s->local = 1;
   s=&cfg->varlist[HPGeCL      ]; s->offset =  7; s->dtype = 0; s->local = 1;
   s=&cfg->varlist[HPGePU      ]; s->offset = 27; s->dtype = 0;
   s=&cfg->varlist[HPGeIP      ]; s->offset = 11; s->dtype = 0;
   s=&cfg->varlist[HPGeDT      ]; s->offset = 32; s->dtype = 0;
   s=&cfg->varlist[HPGeEU      ]; s->offset =  4; s->dtype = 0;
   s=&cfg->varlist[HPGeAU      ]; s->offset = 35; s->dtype = 0;
   s=&cfg->varlist[GRGTHETA    ]; s->offset = 44; s->dtype = 0;
   s=&cfg->varlist[GRGPHI      ]; s->offset = 45; s->dtype = 0;
   s=&cfg->varlist[CLVTHETA    ]; s->offset = 46; s->dtype = 0;
   s=&cfg->varlist[CLVPHI      ]; s->offset = 47; s->dtype = 0;
   s=&cfg->varlist[SEPE        ]; s->offset =  4; s->dtype = 2;
   s=&cfg->varlist[SEPT        ]; s->offset = 24; s->dtype = 2;
   s=&cfg->varlist[SEPTS       ]; s->offset = 34; s->dtype = 2;
   s=&cfg->varlist[SEPPH       ]; s->offset =  8; s->dtype = 2;
   s=&cfg->varlist[SEPPU       ]; s->offset = 27; s->dtype = 2;
   s=&cfg->varlist[SEPTHETA    ]; s->offset = 44; s->dtype = 2;
   s=&cfg->varlist[SEPPH       ]; s->offset = 45; s->dtype = 2;
   s=&cfg->varlist[SEPNUM      ]; s->offset =  3; s->dtype = 2; s->local = 1;
   s=&cfg->varlist[PACE        ]; s->offset =  4; s->dtype = 5;
   s=&cfg->varlist[PACT        ]; s->offset = 24; s->dtype = 5;
   s=&cfg->varlist[PACTS       ]; s->offset = 34; s->dtype = 5;
   s=&cfg->varlist[PACPH       ]; s->offset =  8; s->dtype = 5;
   s=&cfg->varlist[PACPU       ]; s->offset = 27; s->dtype = 5;
   s=&cfg->varlist[PACTHETA    ]; s->offset = 44; s->dtype = 5;
   s=&cfg->varlist[PACPHI      ]; s->offset = 45; s->dtype = 5;
   s=&cfg->varlist[PACNUM      ]; s->offset =  3; s->dtype = 5; s->local = 1;
   s=&cfg->varlist[LBLE        ]; s->offset =  4; s->dtype = 3;
   s=&cfg->varlist[LBLT        ]; s->offset = 24; s->dtype = 3;
   s=&cfg->varlist[LBLTS       ]; s->offset = 34; s->dtype = 3;
   s=&cfg->varlist[LBLPH       ]; s->offset =  8; s->dtype = 3;
   s=&cfg->varlist[LBLPU       ]; s->offset = 27; s->dtype = 3;
   s=&cfg->varlist[LBLTHETA    ]; s->offset = 44; s->dtype = 3;
   s=&cfg->varlist[LBLPHI      ]; s->offset = 45; s->dtype = 3;
   s=&cfg->varlist[LBLNUM      ]; s->offset =  3; s->dtype = 3; s->local = 1;
   s=&cfg->varlist[LBT         ]; s->offset =  4; s->dtype = 4;
   s=&cfg->varlist[LBTT        ]; s->offset = 24; s->dtype = 4;
   s=&cfg->varlist[LBTTS       ]; s->offset = 34; s->dtype = 4;
   s=&cfg->varlist[LBTPH       ]; s->offset =  8; s->dtype = 4;
   s=&cfg->varlist[LBTPU       ]; s->offset = 27; s->dtype = 4;
   s=&cfg->varlist[LBTNUM      ]; s->offset =  3; s->dtype = 4; s->local = 1;
   s=&cfg->varlist[GRSE        ]; s->offset =  4; s->dtype = 7;
   s=&cfg->varlist[GRST        ]; s->offset = 24; s->dtype = 7;
   s=&cfg->varlist[GRSTS       ]; s->offset = 34; s->dtype = 7;
   s=&cfg->varlist[GRSPH       ]; s->offset =  8; s->dtype = 7;
   s=&cfg->varlist[GRSPU       ]; s->offset = 27; s->dtype = 7;
   s=&cfg->varlist[GRSNUM      ]; s->offset = 40; s->dtype = 7; s->local = 1;
   s=&cfg->varlist[GRSPOS      ]; s->offset =  3; s->dtype = 7;
   s=&cfg->varlist[GRSTYPE     ]; s->offset = 44; s->dtype = 7;
   s=&cfg->varlist[LBSE        ]; s->offset =  4; s->dtype = 8;
   s=&cfg->varlist[LBST        ]; s->offset = 24; s->dtype = 8;
   s=&cfg->varlist[LBSTS       ]; s->offset = 34; s->dtype = 8;
   s=&cfg->varlist[LBSPH       ]; s->offset =  8; s->dtype = 8;
   s=&cfg->varlist[LBSPU       ]; s->offset = 27; s->dtype = 8;
   s=&cfg->varlist[LBSNUM      ]; s->offset = 40; s->dtype = 8; s->local = 1;
   s=&cfg->varlist[LBSPOS      ]; s->offset =  3; s->dtype = 8; s->local = 1;
   s=&cfg->varlist[MIDAS_Time  ]; s->offset = -1;
   s=&cfg->varlist[TD_GRG_GRG  ]; s->offset = -1;
   s=&cfg->varlist[TD_SEP_SEP  ]; s->offset = -1;
   s=&cfg->varlist[TD_PAC_PAC  ]; s->offset = -1;
   s=&cfg->varlist[TD_LBL_LBL  ]; s->offset = -1;
   s=&cfg->varlist[TD_GRG_SEP  ]; s->offset = -1;
   s=&cfg->varlist[TSD_GRG_GRG ]; s->offset = -1;
   s=&cfg->varlist[TSD_SEP_SEP ]; s->offset = -1;
   s=&cfg->varlist[TSD_PAC_PAC ]; s->offset = -1;
   s=&cfg->varlist[TSD_LBL_LBL ]; s->offset = -1;
   s=&cfg->varlist[TSD_GRG_SEP ]; s->offset = -1;
   s=&cfg->varlist[TSD_GRG_SEP ]; s->offset = -1;
   s=&cfg->varlist[TD_GRG_PAC  ]; s->offset = -1;
   s=&cfg->varlist[TSD_GRG_PAC ]; s->offset = -1;
   s=&cfg->varlist[TD_SEP_PAC  ]; s->offset = -1;
   s=&cfg->varlist[TSD_SEP_PAC ]; s->offset = -1;
   s=&cfg->varlist[TD_GRG_LBL  ]; s->offset = -1;
   s=&cfg->varlist[TSD_GRG_LBL ]; s->offset = -1;
   s=&cfg->varlist[TD_SEP_LBL  ]; s->offset = -1;
   s=&cfg->varlist[TSD_SEP_LBL ]; s->offset = -1;
   s=&cfg->varlist[ANG_GRG_GRG ]; s->offset = -1;
   s=&cfg->varlist[ANG_CLV_CLV ]; s->offset = -1;
   s=&cfg->varlist[ANG_SEP_SEP ]; s->offset = -1;
   s=&cfg->varlist[ANG_PAC_PAC ]; s->offset = -1;
   s=&cfg->varlist[ANG_LBL_LBL ]; s->offset = -1;
   s=&cfg->varlist[ANG_GRG_SEP ]; s->offset = -1;
   s=&cfg->varlist[ANG_GRG_PAC ]; s->offset = -1;
   s=&cfg->varlist[ANG_GRG_LBL ]; s->offset = -1;
   s=&cfg->varlist[ANG_PAC_LBL ]; s->offset = -1;
   s=&cfg->varlist[ANG_SEP_LBL ]; s->offset = -1;
   s=&cfg->varlist[PPG_NUM     ]; s->offset = -1;
   s=&cfg->varlist[PPG_TIME    ]; s->offset = -1;
   s=&cfg->varlist[PPG_PAT     ]; s->offset = -1;
   return(0);
}
