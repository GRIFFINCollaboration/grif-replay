#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dragon-format.h"
#include "config.h"
#include "histogram.h"

//#######################################################################
//#####              USER SORT - custom user histograms             #####
//#######################################################################

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
extern Dragon_event evbuf[EVT_BUFSIZE];
int user_sort(int win_strt, int win_end, int flag)
{
   Dragon_event *alt, *ptr = &evbuf[win_strt];
   Config *cfg = configs[1]; // ** Sort is using config[1]
   int i, j, abs_dt, dt, win_idx, user_cnt = cfg->nuser;
   Sortvar *var, *yvar;
   float value, yval;
   Histogram *histo;
   Gate *gate;
   Cond *cond;

   //if( user_cnt == 0 ){ return(0); } move this to caller

   if( win_strt != win_end ){ // multiple fragments in window
      for(i=0; i<cfg->nuser; i++){ cfg->user_histos[i]->done_flag = 0; }
   } 

   // when a histo is done for this event, dec user_cnt, and stop when zero
   if( (win_idx = win_strt+1) == EVT_BUFSIZE ){ win_idx = 0; } 
   if( ++win_end == EVT_BUFSIZE ){ win_end = 0; } // need to include win_end
   while( win_idx != win_end && user_cnt > 0 ){
      alt = &evbuf[win_idx];
      if( ++win_idx == EVT_BUFSIZE ){ win_idx = 0; }

      abs_dt = dt = ptr->ts - alt->ts; if( dt < 0 ){ abs_dt = -1*dt; }
      // if( abs_dt > global_window_size ){ break; }

      // clear gates and conditions
      for(i=0; i<cfg->ngates; i++){ gate = cfg->gatelist[i];
         if( gate->use_count == 0 ){ continue; }
         //for(j=0; j<gate->nconds; j++){ gate->conds[i]->passed = 0; }
         gate->valid=0;
      }
      for(i=0; i<cfg->nconds; i++){ cond = cfg->condlist[i];
         cond->valid = cond->passed = 0;
      }

      for(j=0; j<cfg->nuser; j++){ histo = cfg->user_histos[j];
         if( histo->done_flag != 0 ){ continue; }
         var = histo->xvar;
         //if( var->subsys1 != ptr->subsys ){ continue; }
         for(i=0; i<histo->num_gates; i++){ gate = histo->gatelist[i];
            if( ! gate->valid ){ test_gate(ptr, alt, gate); }
            if( gate->passed == 0 ){ break; }
         }
         if( i < histo->num_gates ){ continue; } // gates not all passed
         if( histo->type == INT_1D ){
            value = var->get_value(ptr, alt, var->subsys1, var->subsys2);
            histo->Fill(histo, value, 1);
            histo->done_flag = 1;  --user_cnt;
         } else if( histo->type == INT_2D ){ yvar = histo->yvar;
            // 2d histos - already checked all gates + set valids
            //if( yvar->subsys1 != alt->subsys ){ continue; }
            yval = var->get_value(alt, ptr, yvar->subsys1, yvar->subsys2);
            ((TH2I *)histo)->Fill((TH2I *)histo, value, yval, 1);
            // do not set done_flag for true 2d histos
         }
      }
   }
   user_removefrom_window(win_strt, win_end);
   return(0);
}

int test_gate(Dragon_event *ptr, Dragon_event *alt, Gate *gate)
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

int user_addto_window(int win_strt, int new_frag){ return(0); }
int user_removefrom_window(int win_strt, int new_frag){ return(0); }

//#######################################################################
//#####                     SORT VARIABLES                          #####
//#######################################################################
// to allow simple adding of more below, as required, web interface should request the list of variables

float uval_empty    (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){ printf("Unimplemented variable Value\n"); return( 0 ); }
float uval_sb0      (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){ }
float uval_septof   (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){ }
float uval_icsum    (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){ return( ptr->type == s1 && alt->type == s2 ? ptr->head_tail_data.tail_data.ic_sum : -1 ); }
float uval_icanode0 (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){ }
float uval_icanode1 (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){ }
float uval_bgoes0   (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){
   return( ptr->type == s1 && alt->type == s2 ? ptr->head_tail_data.head_data.bgo_e0 : -1 );
}
float uval_bgoz0  (Dragon_event *ptr, Dragon_event *alt, int s1, int s2){
   int ch = ptr->head_tail_data.head_data.bgo_ch0;
   if( ptr->type != s1 || alt->type != s2 || ch < 0 || ch >= BGO_MAXCHAN ){ return (-1); }
   return( bgo_tdc_zposn[ch] );
}

typedef float(*dragon_get_value)(Dragon_event *, Dragon_event *, int, int);
typedef float(*config_get_value)(void *, void *, int, int);

struct sortvar_desc_struct { // float32 good for integers to 16.8M
   dragon_get_value get_value;
   char name[32];  char desc[256];  int subsys1;  int subsys2;
};
#define NUM_SORTVARS 93
struct sortvar_desc_struct sortvarlist[NUM_SORTVARS]={
 // HPGE
   {&uval_sb0,      "SB_E0",       "Surface Barrier Channel 0 Energy",                              HEAD_EVENT, TAIL_EVENT},
   {&uval_septof,   "Sep_TOF",     "Separator Time of Flight between Gamma and HeavyIon Detection", HEAD_EVENT, TAIL_EVENT},
   {&uval_icsum,    "IC_Sum",      "Ion Chamber Summed Energy",                                     HEAD_EVENT, TAIL_EVENT},
   {&uval_icanode0, "IC_E0",       "Ion Chamber Anode#0 Energy",                                    HEAD_EVENT, TAIL_EVENT},
   {&uval_icanode1, "IC_E1",       "Ion Chamber Anode#1 Energy",                                    HEAD_EVENT, TAIL_EVENT},
   {&uval_bgoes0,   "BGO_Esort0",  "BGO Largest Energy",                                            HEAD_EVENT, TAIL_EVENT},
   {&uval_bgoz0,    "BGO_Z0",      "BGO Z position with Largest Energy",                            HEAD_EVENT, TAIL_EVENT},
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
