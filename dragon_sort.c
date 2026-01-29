//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#######################################################################

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "config.h"
#include "dragon-format.h"
#include "histogram.h"
#include "odb.h"
//#include "default_sort.h"

extern Dragon_event evbuf[EVT_BUFSIZE];

int presort_window_width = 200;  // 2us [May not even need a presort]
int sort_window_width    = 200;  // 2us - MAXIMUM (indiv. gates can be smaller)

// Default sort function declarations
extern int init_parameters_from_globals(Config *cfg);
extern int init_chan_histos(Config *cfg);
extern int init_histos(Config *cfg);
extern int fill_chan_histos(Dragon_event *ptr);
extern int fill_singles_histos(Dragon_event *ptr);
extern int fill_coinc_histos(int win_idx, int frag_idx);

float spread(int val){ return( val + rand()/(1.0*RAND_MAX) ); }

/////////////////       ODB TABLES    ////////////////////////
// read all the detector information [adc/tdc module/channel number]
//                         also gains/offsets/pedastals etc.
// there are a number of detector types ...
//    bgo mcp dssd ic db nai ge, also "head/tail?"
//
///////////////////////////////////////////////////////////////////////////
// Odb stuff from run 13187 ...
// dragon/bgo/variables/
//   [0-29] - adc/{channel,pedastal,slope[=0.00415],offset[=0]}
//    [0-29] - tdc/{channel,slope[=0.1],offset[=0],position[=0|4|-4],y,z}
// mcp/variables/adc{chan[18-21].module,slope,offset},tac_adc[ch22],tdc[ch6,7]
// dsssd/variables/adc{chan[0-31],modules,slope,offset}
//                 tdc_frnt,back[1ch-each:4,5]
// ic/variables/adc/channel[25-29] tdc[chan0-3]
// sb/variables/adc[ch16,17]
// nai/variables/adc[ch30,31]
// ge/variables/adc[ch23]
// head/variables/xtdc[ch50],rftdc[ch48],tdc0[ch49]
// tail/variables/xtdc[ch10],rftdc[ch8],tdc0[ch9]
// 
//    HEAD - BGO SB[target-monitor]             [NaI Ge ?]
//
//    TAIL - DSSD[Position,E,Tof], Mcp[Tof], IC[mass]
// 
///////////////////////////////////////////////////////////////////////////
// head detectors ...
int   ge_adc_chan, ge_adc_module, head_xtdc_chan, head_rftdc_chan, head_tdc0_chan;
float ge_adc_slope, ge_adc_offset, head_xtdc_slope, head_xtdc_offset;
float head_rftdc_slope, head_rftdc_offset, head_tdc0_slope, head_tdc0_offset;
int   bgo_adc_chan    [BGO_MAXCHAN];
int   bgo_adc_module  [BGO_MAXCHAN];
float bgo_adc_slope   [BGO_MAXCHAN];
float bgo_adc_offset  [BGO_MAXCHAN];
int   bgo_adc_pedestal[BGO_MAXCHAN];
int   bgo_tdc_chan    [BGO_MAXCHAN];
float bgo_tdc_slope   [BGO_MAXCHAN];
float bgo_tdc_offset  [BGO_MAXCHAN];
float bgo_tdc_xposn   [BGO_MAXCHAN];
float bgo_tdc_yposn   [BGO_MAXCHAN];
float bgo_tdc_zposn   [BGO_MAXCHAN];
// also bgo/hv/channel
int   sb_adc_chan     [SB_MAXCHAN];
int   sb_adc_module   [SB_MAXCHAN];
float sb_adc_slope    [SB_MAXCHAN];
float sb_adc_offset   [SB_MAXCHAN];
int   nai_adc_chan    [NAI_MAXCHAN];
int   nai_adc_module  [NAI_MAXCHAN];
float nai_adc_slope   [NAI_MAXCHAN];
float nai_adc_offset  [NAI_MAXCHAN];
// tail detectors ...
int   tail_xtdc_chan,  tail_rftdc_chan, tail_tdc0_chan;
float tail_xtdc_slope, tail_xtdc_offset, tail_rftdc_slope, tail_rftdc_offset;
float tail_tdc0_slope, tail_tdc0_offset; 
int   dssd_adc_chan    [DSSD_MAXCHAN];
int   dssd_adc_module  [DSSD_MAXCHAN]; // 32 dssd chan are in two 16-chan adcs
float dssd_adc_slope   [DSSD_MAXCHAN];
float dssd_adc_offset  [DSSD_MAXCHAN];
int   dssd_tdc_chan    [DSSD_TDCCHAN];
float dssd_tdc_slope   [DSSD_TDCCHAN];
float dssd_tdc_offset  [DSSD_TDCCHAN];
int   ic_adc_chan      [IC_MAXCHAN];
int   ic_adc_module    [IC_MAXCHAN];
float ic_adc_slope     [IC_MAXCHAN];
float ic_adc_offset    [IC_MAXCHAN];
int   ic_tdc_chan      [IC_MAXCHAN];
float ic_tdc_slope     [IC_MAXCHAN];
float ic_tdc_offset    [IC_MAXCHAN];
int   mcp_adc_chan     [MCP_MAXCHAN];
int   mcp_adc_module   [MCP_MAXCHAN];
float mcp_adc_slope    [MCP_MAXCHAN];
float mcp_adc_offset   [MCP_MAXCHAN];
int   mcptac_adc_chan;
int   mcptac_adc_module;
float mcptac_adc_slope;
float mcptac_adc_offset;
int   mcp_tdc_chan     [MCP_TDCCHAN];
float mcp_tdc_slope    [MCP_TDCCHAN]; // was called tdc0,tdc1
float mcp_tdc_offset   [MCP_TDCCHAN];
// also coinc/variables/window/buffer_time

// derived tables mapping adc/tdc channels to specific detector-type/channel
int head_adc_dettype[2*V792_MAXCHAN]; int head_adc_dstchan[2*V792_MAXCHAN]; 
int head_tdc_dettype[ V1190_MAXCHAN]; int head_tdc_dstchan[ V1190_MAXCHAN]; 
int tail_adc_dettype[2*V792_MAXCHAN]; int tail_adc_dstchan[2*V792_MAXCHAN]; 
int tail_tdc_dettype[ V1190_MAXCHAN]; int tail_tdc_dstchan[ V1190_MAXCHAN];

static char subsys_name[16][8]={
   "   BGO", "    SB", "   NAI", "  DSSD",   "    IC", "   MCP", " MCPTAC", "    GE", "  XTDC",
   " RFTDC", "  TDC0", "", "",   "", "", ""
};

int add_v1190_item(int v1190_chan, int subsys_type, int subsys_chan, int *tdc_dettype, int *tdc_dstchan)
{
   if( v1190_chan < 0 || v1190_chan >= V1190_MAXCHAN ){
      printf("read_dragon_odb: invalid %s v1190 channel number:%d\n",
             subsys_name[subsys_type], v1190_chan); return(-1); 
   }
   if( tdc_dettype[v1190_chan] != -1 ){
      printf("read_dragon_odb: multiple assignments for v1190 chan:%d\n",
             v1190_chan);  return(-1); 
   }
   tdc_dettype[v1190_chan] = subsys_type;  tdc_dstchan[v1190_chan] = subsys_chan;
   return(0);
}
int add_v792_item(int v792_chan, int v792_module, int subsys_type, int subsys_chan, int *adc_dettype, int *adc_dstchan)
{
   printf("Add V792 Chan[%d,%d] Sys[%d] Dstchan[%d]\n", v792_chan, v792_module, subsys_type, subsys_chan);
   if( v792_module < 0 || v792_module > 1 ){
      printf("read_dragon_odb: invalid %s v792 module number:%d\n",
             subsys_name[subsys_type], v792_module); return(-1); 
   }
   if( v792_chan < 0 || v792_chan >= V792_MAXCHAN ){
      printf("read_dragon_odb: invalid %s v792 channel number:%d\n",
             subsys_name[subsys_type], v792_chan); return(-1); 
   }
   v792_chan += V792_MAXCHAN*v792_module;
   if( adc_dettype[v792_chan] != -1 ){
      printf("read_dragon_odb: multiple assignments for v792 chan:%d\n",
             v792_chan);
      // since dragon odb is wrong - allow later assignments to overwrite
      // return(-1); 
   }
   adc_dettype[v792_chan] = subsys_type;  adc_dstchan[v792_chan] = subsys_chan;
   return(0);
}

int read_dragon_odb(int bank_len, int *bank_data)
{
   void *array; int i, chan, mod, type, nval, size, err=0;
   int *i_array; float *f_array;
   
   read_odb_tree(bank_len, bank_data);
   
   ge_adc_chan = ge_adc_module = head_xtdc_chan = head_rftdc_chan = head_tdc0_chan = -1;
   odbval_int  ("dragon/ge/variables/adc/channel",    &ge_adc_chan);
   odbval_int  ("dragon/ge/variables/adc/module",     &ge_adc_module);
   odbval_float("dragon/ge/variables/adc/slope",      &ge_adc_slope);
   odbval_float("dragon/ge/variables/adc/offset",     &ge_adc_offset);
   odbval_int  ("dragon/head/variables/xtdc/channel" ,&head_xtdc_chan);
   odbval_float("dragon/head/variables/xtdc/slope",   &head_xtdc_slope);
   odbval_float("dragon/head/variables/xtdc/offset",  &head_xtdc_offset);
   odbval_int  ("dragon/head/variables/rf_tdc/channel", &head_rftdc_chan);
   odbval_float("dragon/head/variables/rf_tdc/slope",  &head_rftdc_slope);
   odbval_float("dragon/head/variables/rf_tdc/offset", &head_rftdc_offset);
   odbval_int  ("dragon/head/variables/tdc0/channel", &head_tdc0_chan);
   odbval_float("dragon/head/variables/tdc0/slope",   &head_tdc0_slope);
   odbval_float("dragon/head/variables/tdc0/offset",  &head_tdc0_offset);

   memset(bgo_adc_chan, -1, BGO_MAXCHAN*sizeof(int) );
   memset(bgo_tdc_chan, -1, BGO_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/bgo/variables/adc/channel",
                         bgo_adc_chan,      BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/adc/slope",
                         bgo_adc_slope,     BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/adc/offset",
                         bgo_adc_offset,   BGO_MAXCHAN);
   err += odbarray_int  ("dragon/bgo/variables/adc/pedestal",
                         bgo_adc_pedestal, BGO_MAXCHAN);
   err += odbarray_int  ("dragon/bgo/variables/tdc/channel",
                         bgo_tdc_chan,     BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/slope",
                         bgo_adc_slope,    BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/offset",
                         bgo_adc_offset,   BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/position/x",
                         bgo_tdc_xposn,    BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/position/y",
                         bgo_tdc_yposn,    BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/position/z",
                         bgo_tdc_zposn,    BGO_MAXCHAN);
   
   memset(sb_adc_chan,  -1,  SB_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/sb/variables/adc/channel",
                         sb_adc_chan,      SB_MAXCHAN);
   err += odbarray_int  ("dragon/sb/variables/adc/module",
                         sb_adc_module,      SB_MAXCHAN);
   err += odbarray_float("dragon/sb/variables/adc/slope",
                         sb_adc_slope,     SB_MAXCHAN);
   err += odbarray_float("dragon/sb/variables/adc/offset",
                         sb_adc_offset,   SB_MAXCHAN);

   memset(nai_adc_chan,  -1,  NAI_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/nai/variables/adc/channel",
                         nai_adc_chan,      NAI_MAXCHAN);
   err += odbarray_int  ("dragon/nai/variables/adc/module",
                         nai_adc_module,      NAI_MAXCHAN);
   err += odbarray_float("dragon/nai/variables/adc/slope",
                         nai_adc_slope,     NAI_MAXCHAN);
   err += odbarray_float("dragon/nai/variables/adc/offset",
                         nai_adc_offset,   NAI_MAXCHAN);

   tail_xtdc_chan = tail_rftdc_chan = tail_tdc0_chan = -1;
   odbval_int  ("dragon/tail/variables/xtdc/channel" ,&tail_xtdc_chan);
   odbval_float("dragon/tail/variables/xtdc/slope",   &tail_xtdc_slope);
   odbval_float("dragon/tail/variables/xtdc/offset",  &tail_xtdc_offset);
   odbval_int  ("dragon/tail/variables/rf_tdc/channel", &tail_rftdc_chan);
   odbval_float("dragon/tail/variables/rf_tdc/slope",  &tail_rftdc_slope);
   odbval_float("dragon/tail/variables/rf_tdc/offset", &tail_rftdc_offset);
   odbval_int  ("dragon/tail/variables/tdc0/channel", &tail_tdc0_chan);
   odbval_float("dragon/tail/variables/tdc0/slope",   &tail_tdc0_slope);
   odbval_float("dragon/tail/variables/tdc0/offset",  &tail_tdc0_offset);

   memset(dssd_adc_chan,  -1,  DSSD_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/dsssd/variables/adc/channel",
                         dssd_adc_chan,      DSSD_MAXCHAN);
   err += odbarray_int  ("dragon/dsssd/variables/adc/module",
                         dssd_adc_module,    DSSD_MAXCHAN);
   err += odbarray_float("dragon/dsssd/variables/adc/slope",
                         dssd_adc_slope,     DSSD_MAXCHAN);
   err += odbarray_float("dragon/dsssd/variables/adc/offset",
                         dssd_adc_offset,    DSSD_MAXCHAN);
   dssd_tdc_chan[0] = dssd_tdc_chan[1] = -1;
   odbval_int  ("dragon/dsssd/variables/tdc_front/channel", &dssd_tdc_chan  [0]);
   odbval_int  ("dragon/dsssd/variables/tdc_back/channel",  &dssd_tdc_chan  [1]);
   odbval_float("dragon/dsssd/variables/tdc_front/slope",   &dssd_tdc_slope [0]);
   odbval_float("dragon/dsssd/variables/tdc_front/offset",  &dssd_tdc_offset[1]);
   odbval_float("dragon/dsssd/variables/tdc_back/slope",    &dssd_tdc_slope [0]);
   odbval_float("dragon/dsssd/variables/tdc_back/offset",   &dssd_tdc_offset[1]);

   memset(ic_adc_chan,  -1,  IC_MAXCHAN*sizeof(int) );
   memset(ic_tdc_chan,  -1,  IC_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/ic/variables/adc/channel",
                         ic_adc_chan,      IC_MAXCHAN);
   err += odbarray_int  ("dragon/ic/variables/adc/channel",
                         ic_adc_chan,      IC_MAXCHAN);
   err += odbarray_int  ("dragon/ic/variables/tdc/channel",
                         ic_tdc_chan,      IC_MAXCHAN);
   err += odbarray_float("dragon/ic/variables/adc/slope",
                         ic_adc_slope,     IC_MAXCHAN);
   err += odbarray_float("dragon/ic/variables/adc/offset",
                         ic_adc_offset,   IC_MAXCHAN);

   memset(mcp_adc_chan,  -1,  MCP_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/mcp/variables/adc/channel",
                         mcp_adc_chan,      MCP_MAXCHAN);
   err += odbarray_int  ("dragon/mcp/variables/adc/module",
                         mcp_adc_module,    MCP_MAXCHAN);
   err += odbarray_float("dragon/mcp/variables/adc/slope",
                         mcp_adc_slope,    MCP_MAXCHAN);
   err += odbarray_float("dragon/mcp/variables/adc/offset",
                         mcp_adc_offset,    MCP_MAXCHAN);

   memset(mcp_tdc_chan,  -1,  MCP_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/mcp/variables/tdc/channel", mcp_tdc_chan, MCP_TDCCHAN);
   err += odbarray_float("dragon/mcp/variables/tdc/slope",   mcp_tdc_slope, MCP_TDCCHAN);
   err += odbarray_float("dragon/mcp/variables/tdc/offset",  mcp_tdc_offset, MCP_TDCCHAN);

   mcptac_adc_chan = -1;
   odbval_int  ("dragon/mcp/variables/tac_adc/channel", &mcptac_adc_chan);
   odbval_int  ("dragon/mcp/variables/tac_adc/module",  &mcptac_adc_module);
   odbval_float("dragon/mcp/variables/tac_adc/slope",   &mcptac_adc_slope);
   odbval_float("dragon/mcp/variables/tac_adc/offset",  &mcptac_adc_offset);
   
   if( err ){ printf("read_dragon_odb: %d errors\n", err); }

   // generate derived tables ---------------------------------
   memset(head_adc_dettype, -1,   V792_MAXCHAN*sizeof(int) );
   memset(head_adc_dstchan, -1,   V792_MAXCHAN*sizeof(int) );
   memset(head_tdc_dettype, -1,  V1190_MAXCHAN*sizeof(int) );
   memset(head_tdc_dstchan, -1,  V1190_MAXCHAN*sizeof(int) );
   memset(tail_adc_dettype, -1, 2*V792_MAXCHAN*sizeof(int) );
   memset(tail_adc_dstchan, -1, 2*V792_MAXCHAN*sizeof(int) );
   memset(tail_tdc_dettype, -1,  V1190_MAXCHAN*sizeof(int) );
   memset(tail_tdc_dstchan, -1,  V1190_MAXCHAN*sizeof(int) );

   ////////////////////////////////// HEAD ////////////////////////////////////////
   for(i=0; i<BGO_MAXCHAN; i++){ // BGO first due to odb channel assignment errors
      if( (chan = bgo_adc_chan[i]) != -1 ){
         add_v792_item (chan, bgo_adc_module[i], SUBSYS_BGO, i, head_adc_dettype, head_adc_dstchan);
      }
      if( (chan = bgo_tdc_chan[i]) != -1 ){
         add_v1190_item(chan, SUBSYS_BGO, i,head_tdc_dettype,head_tdc_dstchan);
      }
   }
   if( (chan = ge_adc_chan) != -1 ){
      add_v792_item (chan, 0, SUBSYS_GE, 0,head_adc_dettype,head_adc_dstchan);
   }
   if( (chan = head_xtdc_chan) != -1 ){
      add_v1190_item(chan, SUBSYS_XTDC, 0,head_tdc_dettype,head_tdc_dstchan);
   }
   if( (chan = head_rftdc_chan) != -1 ){
      add_v1190_item(chan, SUBSYS_RFTDC, 0,head_tdc_dettype,head_tdc_dstchan);
   }
   if( (chan = head_tdc0_chan) != -1 ){
      add_v1190_item(chan, SUBSYS_TDC0, 0,head_tdc_dettype,head_tdc_dstchan);
   }
   for(i=0; i<SB_MAXCHAN; i++){
      if( (chan = sb_adc_chan[i]) != -1 ){
         add_v792_item (chan, sb_adc_module[i], SUBSYS_SB, i, head_adc_dettype, head_adc_dstchan);
      }
   }
   for(i=0; i<NAI_MAXCHAN; i++){
      if( (chan = nai_adc_chan[i]) != -1 ){
         add_v792_item (chan, nai_adc_module[i], SUBSYS_NAI, i, head_adc_dettype, head_adc_dstchan);
      }
   }
   ////////////////////////////////// TAIL ///////////////////////////////////
   if( (chan = tail_xtdc_chan) != -1 ){
      add_v1190_item(chan, SUBSYS_XTDC, 0,tail_tdc_dettype, tail_tdc_dstchan);
   }
   if( (chan = tail_rftdc_chan) != -1 ){
      add_v1190_item(chan, SUBSYS_RFTDC, 0,tail_tdc_dettype, tail_tdc_dstchan);
   }
   if( (chan = tail_tdc0_chan) != -1 ){
      add_v1190_item(chan, SUBSYS_TDC0, 0,tail_tdc_dettype, tail_tdc_dstchan);
   }
   for(i=0; i<DSSD_MAXCHAN; i++){
      if( (chan = dssd_adc_chan[i]) == -1 ){ continue; }
      add_v792_item (chan, dssd_adc_module[i], SUBSYS_DSSD, i, tail_adc_dettype, tail_adc_dstchan);
   }
   for(i=0; i<DSSD_TDCCHAN; i++){
      if( (chan = dssd_tdc_chan[i]) == -1 ){ continue; }
      add_v1190_item(chan, SUBSYS_DSSD, i, tail_tdc_dettype, tail_tdc_dstchan);
   }
   for(i=0; i<IC_MAXCHAN; i++){
      if( (chan = ic_adc_chan[i]) != -1 ){
         add_v792_item (chan, ic_adc_module[i], SUBSYS_IC, i, tail_adc_dettype, tail_adc_dstchan);
      }     
      if( (chan = ic_tdc_chan[i]) != -1 ){
         add_v1190_item(chan, SUBSYS_IC, i,tail_tdc_dettype,tail_tdc_dstchan);
      }      
   }
   for(i=0; i<MCP_MAXCHAN; i++){
      if( (chan = mcp_adc_chan[i]) == -1 ){ continue; }
      add_v792_item (chan, mcp_adc_module[i], SUBSYS_MCP, i, tail_adc_dettype, tail_adc_dstchan);   
   }
   for(i=0; i<MCP_TDCCHAN; i++){
      if( (chan = mcp_tdc_chan[i]) == -1 ){ continue; }
      add_v1190_item(chan, SUBSYS_MCP, i, tail_tdc_dettype, tail_tdc_dstchan);
   }
   if( (chan = mcptac_adc_chan ) != -1 ){
      add_v792_item(chan, mcptac_adc_module, SUBSYS_MCPTAC, 0, tail_adc_dettype, tail_adc_dstchan);
   }

   return(0);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int init_default_histos(Config *cfg, Sort_status *arg)
{
   Cal_coeff *cal;
   int i, j;

   init_parameters_from_globals(cfg);
   init_chan_histos(cfg);
   init_histos(cfg);

   return(0);
}

int init_parameters_from_globals(Config *cfg)
{
   // Here set parameters from the Globals
   // time-difference conditions for coincidences
   // PRESORT time-difference Conditions
   Global *global;
   char tmp[32];
   int i, j;

   // Initialize all time differences between subsystems to the default 250ns

   // Intialize all PRESORT timing windows to their defaults
    
   // Search globals for time difference settings and overwrite their values
   for(i=0; i<cfg->nglobal; i++){
      global = cfg->globals[i];
      sprintf(tmp,"%s",global->name);
      if(strncmp(tmp,"time_diff_",10) == 0){
        // This global is a time difference value
        // figure out which one and store value
      }
   }// end of for(i=0; i<cfg->nglobal; i++){
   return(0);
}

//#######################################################################
//######## PRESORT(gain corrections, addback, suppression)     ##########
//#######################################################################
// There are two coincidence windows - presort-window and sort-window
// For griffin, the presort was used to sum scattered Gamma Energies
//   to recover their full energy, and also to suppress Gammas that scattered
//   out of detectors into the suppression shields
// The sort-window was used for coincidences, using cleaned up full energies
// - with partial scatters etc. already removed during the presort
// The two window arrangement greatly simplified the code
// 
// The presort only deals with one event per call - either the first
// or last event in the presort-window, in order to check coincidences with
// events preceeding or following the event of interest, so there are separate
// presort_enter/leave functions.
//
// The sort-window is still only processed once, relying on loops, and
// swapping orders of events to deal with all coincidences
// - the code could probably be simplified by also calling twice?
//
//
// pre_sort_enter used to be called apply_gains
// it is the first function to be called when processing an event
// event:frag_idx is the "current event", all other events are BEFORE it
int pre_sort_enter(int start_idx, int frag_idx)
{
   Dragon_event *ptr = &evbuf[frag_idx];
   Head_data *head = &ptr->head_tail_data.head_data;
   Tail_data *tail = &ptr->head_tail_data.tail_data;
   float energy, ecal;
   int dt;
   if(        ptr->type == HEAD_EVENT ){
      ;
   } else if( ptr->type == TAIL_EVENT ){
      ;
   } 
   return(0);
}

// Presort - do Suppression and Addback here
//  - frag_idx is about to leave coinc window (which ends at end_idx)
//    check other frags in window for possible suppression and/or summing
//  also calculate multiplicities[store in frag_idx only]
int pre_sort_exit(int frag_idx, int end_idx)
{
   Dragon_event *alt, *ptr = &evbuf[frag_idx];
   int i, j, dt;

   i = frag_idx;
   while( i != end_idx ){ // need at least two events in window
      if( ++i >= EVT_BUFSIZE ){ i=0; } alt = &evbuf[i]; // WRAP
      // Determine absolute time difference between timestamps
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }
   }// end of while
   return(0);
}

int default_sort(int win_idx, int frag_idx, int flag)
{
   Dragon_event *ptr;
   int i;

   // sort first event, even if window only contains that event
   // (usually require at least two)
   for(i=win_idx; ; i++){ ptr = &evbuf[i];
      if( i >= EVT_BUFSIZE ){ i=0; } // WRAP
      if( i != win_idx && flag == SORT_ONE ){ break; }
      if( i != win_idx && i==frag_idx ){ break; }
      fill_chan_histos(ptr);
      fill_singles_histos(ptr);
      if( i==frag_idx ){ break; }
    }
    fill_coinc_histos(win_idx, frag_idx);
    return(0);
}

//#######################################################################
//########        Individual channel singles HISTOGRAMS        ##########
//#######################################################################
// Hitpatterns
#define N_HITPAT  4
#define HIT_BINS  256
char hit_handles[N_HITPAT][32]={ "q_hit_h","t_hit_h", "q_hit_t","t_hit_t" };
char   hit_names[N_HITPAT][32]={ "ADC_head", "TDC_head", "ADC_tail", "TDC_tail" };

TH1I  *hit_hist[N_HITPAT];

#define ADC_CHAN 64
#define ADC_BINS 8192
TH1I  *adc_hist_head[ADC_CHAN];
TH1I  *adc_hist_tail[ADC_CHAN];

#define TDC_CHAN 64
#define TDC_BINS 8192
TH1I  *tdc_hist_head[ADC_CHAN];
TH1I  *tdc_hist_tail[ADC_CHAN];


int init_chan_histos(Config *cfg)
{  // 1d histograms for Q,E,T,Wf for each channel in odb
   char title[STRING_LEN], handle[STRING_LEN];
   int i, j, k, pos;

   open_folder(cfg, "Singles");
   open_folder(cfg, "ADCs");
   for(i=0; i<ADC_CHAN; i++){ // ADC singles spectra
      sprintf(title,  "ADC_Head_Chan_%d", i );
      adc_hist_head[i] = H1_BOOK(cfg, title, title, ADC_BINS, 0, ADC_BINS);
      sprintf(title,  "ADC_Tail_Chan_%d", i );
      adc_hist_tail[i] = H1_BOOK(cfg, title, title, ADC_BINS, 0, ADC_BINS);
   }
   close_folder(cfg);
   open_folder(cfg, "TDCs");
   for(i=0; i<TDC_CHAN; i++){ // TDC singles spectra
      sprintf(title,  "TDC_Head_Chan_%d", i );
      tdc_hist_head[i] = H1_BOOK(cfg, title, title, TDC_BINS, 0, TDC_BINS);
      sprintf(title,  "TDC_Tail_Chan_%d", i );
      tdc_hist_tail[i] = H1_BOOK(cfg, title, title, TDC_BINS, 0, TDC_BINS);
   }
   close_folder(cfg);
   close_folder(cfg);
   open_folder(cfg, "Hits_and_Sums");
   open_folder(cfg, "Hits");
   for(i=0; i<N_HITPAT; i++){ // Create Hitpattern spectra
      sprintf(title,  "Hitpattern_%s",    hit_names[i] );
      hit_hist[i] = H1_BOOK(cfg, hit_handles[i], title, HIT_BINS, 0, HIT_BINS);
   }
   close_folder(cfg);
   close_folder(cfg);
   return(0);
}

int fill_chan_histos(Dragon_event *ptr)
{
   static int event;
   int i, j, count, chan, val, *data;

   //if( ++event < 16384 ){
   //   //ts_hist -> Fill(ts_hist, event,  (int)(ptr->ts/100));
   //} else if( (event % 1000) == 0 ){
   //   //ts_hist -> Fill(ts_hist, 16367+(int)(event/1000),  (int)(ptr->ts/100));
   //}

   //hit_hist[0]   -> Fill(hit_hist[0],    chan,            1);
   //if( ptr->ecal        >= 1 ){ hit_hist[1] -> Fill(hit_hist[1], chan, 1);
   //hit_hist[6] -> Fill(hit_hist[6], ptr->dtype, 1);

   if( ptr->type == HEAD_EVENT ){
      if( (count = ptr->v792a.d_count) <= MAX_V792_WORDS ){
         data = &ptr->v792a.data[0];
         for(i=0; i<count; i++){
            chan = (data[i] >> 16) & 0x01F;
            val  = (data[i]      ) & 0xFFF;
            if( chan >= 0 && chan < 32 ){
               hit_hist[0] -> Fill(hit_hist[0], chan, 1);
               adc_hist_head[chan] -> Fill(adc_hist_head[chan], val, 1);
            }
         }
      }
      if( (count = ptr->v1190.d_count) <= MAX_V1190_WORDS ){
         data = &ptr->v1190.data[0];
         for(i=0; i<count; i++){
            chan = (data[i] >> 24) & 0x7F;
            val  = (data[i]      ) & 0xFFFFF; val >> 4;
            if( chan >= 0 && chan < 64 ){
               hit_hist[1] -> Fill(hit_hist[1], chan, 1);
               tdc_hist_head[chan] -> Fill(tdc_hist_head[chan], val, 1);
            }
         }
      }
   }
   if( ptr->type == TAIL_EVENT ){
      if( (count = ptr->v792a.d_count) <= MAX_V792_WORDS ){
         data = &ptr->v792a.data[0];
         for(i=0; i<count; i++){
            chan = (data[i] >> 16) & 0x01F;
            val  = (data[i]      ) & 0xFFF;
            if( chan >= 0 && chan < 32 ){
               hit_hist[2] -> Fill(hit_hist[2], chan, 1);
               adc_hist_tail[chan] -> Fill(adc_hist_tail[chan], val, 1);
            }
         }
      }
      if( (count = ptr->v792b.d_count) <= MAX_V792_WORDS ){
      data = &ptr->v792b.data[0];
         for(i=0; i<count; i++){
            chan = (data[i] >> 16) & 0x01F;
            val  = (data[i]      ) & 0xFFF;
            if( chan >= 0 && chan < 32 ){ chan += 32;
               hit_hist[2] -> Fill(hit_hist[2], chan, 1);
               adc_hist_tail[chan] -> Fill(adc_hist_tail[chan], val, 1);
            }
         }
      }
      if( (count = ptr->v1190.d_count) <= MAX_V1190_WORDS ){
         data = &ptr->v1190.data[0];
         for(i=0; i<count; i++){
            chan = (data[i] >> 24) & 0x7F;
            val  = (data[i]      ) & 0xFFFFF; val >> 4;
            if( chan >= 0 && chan < 64 ){
               hit_hist[3] -> Fill(hit_hist[3], chan, 1);
               tdc_hist_tail[chan] -> Fill(tdc_hist_tail[chan], val, 1);
            }
         }
      }
   }
   
   return(0);
}

//#######################################################################
//########             Singless and coinc  HISTOGRAMS          ##########
//#######################################################################


// *** NOTE MOST OF THESE ARE COINCIDENCE HISTOS (gated on various conditions)
//           -> move to proper section

TH1I  *sb_ecal[2]; char *sb_title="";
TH1I  *xtofh;  // gamma to hvy-ion 
TH1I  *ic_sum; 
TH2I  *ic_0v1; 
TH1I  *bgo_zpat; 
TH1I  *bgo_e0; 

int init_histos(Config *cfg)
{
   char hnd[32], title[256];
   static Config *save_cfg;
   int i, j, k;

   if( cfg == NULL ){ cfg = save_cfg; } else { save_cfg = cfg; }
   
   open_folder(cfg, "Singles");
   open_folder(cfg, "Basic");
   for(i=0; i<2; i++){
      sprintf(hnd,  "SB_%d", i);
      sprintf(title,"Surface Barrier (target scattering monitor) %d", i);
      sb_ecal[i] = H1_BOOK(cfg, hnd, title, ADC_BINS, 0, ADC_BINS);
   }
   xtofh    = H1_BOOK(cfg, "XTOFH",    "Separator TOF (gamma->HI)", 5000, -10000, 9999);
   ic_sum   = H1_BOOK(cfg, "IC_SUM",   "Summed Energy Loss in Ion Chamber", ADC_BINS, 0, ADC_BINS);
   bgo_zpat = H1_BOOK(cfg, "BGO_ZPAT", "BGO Z HITPATTERN", 100, -50, 49);
   bgo_e0   = H1_BOOK(cfg, "BGO_E0",   "BGO Energy spectrum [chan0]", ADC_BINS, 0, 100);
   ic_0v1   = H2_BOOK(cfg, "IC_0v1",   "Ion Chamber Energy 0 vs 1", ADC_BINS, 0, ADC_BINS, ADC_BINS, 0, ADC_BINS);
   close_folder(cfg);
   close_folder(cfg);

   return(0);
}

int fill_singles_histos(Dragon_event *ptr)
{
   Head_data *head = &ptr->head_tail_data.head_data;
   Tail_data *tail = &ptr->head_tail_data.tail_data;
   int i, chan, val, cnt;
   float v_cal, sum;

   if( ptr->type == HEAD_EVENT ){
      for(i=0; i<BGO_MAXCHAN; i++){
         if( head->bgo_energy[i] > 1 ){
            if( i == 0 ){ bgo_e0->Fill(bgo_e0, head->bgo_energy[i], 1); }
            bgo_zpat->Fill(bgo_zpat, bgo_tdc_zposn[i], 1);
         }
      }
      for(i=0; i<SB_MAXCHAN; i++){
         sb_ecal[i]->Fill(sb_ecal[i], head->sb_energy[i], 1);
      }
   } else if(  ptr->type == TAIL_EVENT ){
      sum = 0; for(i=0; i<IC_MAXCHAN; i++){ sum += tail->ic_energy[i]; }
      ic_sum->Fill(ic_sum, (int)sum, 1);
      if( tail->ic_energy[0] > 0 && tail->ic_energy[1] > 0 ){
         ic_0v1->Fill(ic_0v1, tail->ic_energy[0], tail->ic_energy[1], 1 );
      }
   }
   return(0);
}

// loop over window and sort all head-tail pairs
int fill_coinc_histos(int win_idx, int frag_idx)
{
   int dt, abs_dt,  pos, c1, c2, index, ptr_swap;
   Dragon_event *alt, *ptr = &evbuf[frag_idx];
   // histogram of coincwin-size
   //dt = (frag_idx - win_idx + 2*EVT_BUFSIZE) %  EVT_BUFSIZE; ++frag_hist[dt];
   while( win_idx != frag_idx ){ // check all conicidences in window
      if( ++win_idx == EVT_BUFSIZE ){ win_idx = 0; } // wrap
      alt = &evbuf[win_idx];
      if( ptr->type == alt->type ){ continue; }// ignore head-head and tail-tail
      if( ptr->type == TAIL_EVENT ){ // swap so ptr is always head and alt is always the tail event
         ptr = &evbuf[win_idx]; alt = &evbuf[frag_idx];
      }
      // apply any time gates (coinditions on dt), and increment head-tail stuff here
      //
      //
      dt = ptr->ts - alt->ts; abs_dt = ( dt < 0 ) ? -1*dt : dt;
      if( abs_dt > 10000 ){ continue; }
      xtofh->Fill(xtofh, dt, 1);
   }
   return(0);
}
