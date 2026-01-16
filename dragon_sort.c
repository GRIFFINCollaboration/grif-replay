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

///////////////////////////////////////////////////////////////////////////
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
int ge_adc_chan;
int head_xtdc_chan, head_rftdc_chan, head_tdc0_chan; 
#define BGO_MAXCHAN 32 // currently 30 are used
int   bgo_adc_chan    [BGO_MAXCHAN];
float bgo_adc_slope   [BGO_MAXCHAN];
float bgo_adc_offset  [BGO_MAXCHAN];
int   bgo_adc_pedastal[BGO_MAXCHAN];
int   bgo_tdc_chan    [BGO_MAXCHAN];
float bgo_tdc_slope   [BGO_MAXCHAN];
float bgo_tdc_offset  [BGO_MAXCHAN];
float bgo_tdc_xposn   [BGO_MAXCHAN];
float bgo_tdc_yposn   [BGO_MAXCHAN];
float bgo_tdc_zposn   [BGO_MAXCHAN];
#define SB_MAXCHAN 2
int   sb_adc_chan    [SB_MAXCHAN];
#define NAI_MAXCHAN 2
int   nai_adc_chan    [NAI_MAXCHAN];
// tail detectors ...
int tail_xtdc_chan, tail_rftdc_chan, tail_tdc0_chan; 
#define DSSD_MAXCHAN 32
int   dssd_adc_chan    [DSSD_MAXCHAN];
int   dssd_adc_module  [DSSD_MAXCHAN]; // there are 32 dssd chan in two 16-chan adcs
float dssd_adc_slope   [DSSD_MAXCHAN];
float dssd_adc_offset  [DSSD_MAXCHAN];
int dssd_front_tdcchan;
int dssd_back_tdcchan;
#define IC_MAXCHAN 4
int   ic_adc_chan    [IC_MAXCHAN];
int   ic_tdc_chan    [IC_MAXCHAN];
#define MCP_MAXCHAN 4
int   mcp_adc_chan    [MCP_MAXCHAN];
int   mcp_adc_module  [MCP_MAXCHAN];
float mcp_adc_slope   [MCP_MAXCHAN];
float mcp_adc_offset  [MCP_MAXCHAN];
int mcp_tac_chan, mcp_tdc_chan0,mcp_tdc_chan1; 

int read_dragon_odb(int bank_len, int *bank_data)
{
   void *array; int i, type, nval, size, err=0;
   int *i_array; float *f_array;
   
   read_odb_tree(bank_len, bank_data);
   
   ge_adc_chan = head_xtdc_chan = head_rftdc_chan = head_tdc0_chan = -1;
   odbval_int("dragon/sb/variables/adc/channel", &ge_adc_chan);
   odbval_int("dragon/head/variables/xtdc", &head_xtdc_chan);
   odbval_int("dragon/head/variables/rftdc", &head_rftdc_chan);
   odbval_int("dragon/head/variables/tdc0", &head_tdc0_chan);

   memset(bgo_adc_chan, -1, BGO_MAXCHAN*sizeof(int) );
   memset(bgo_tdc_chan, -1, BGO_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/bgo/variables/adc/channel",
                         bgo_adc_chan,      BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/adc/slope",
                         bgo_adc_slope,     BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/adc/offset",
                         bgo_adc_offset,   BGO_MAXCHAN);
   err += odbarray_int  ("dragon/bgo/variables/adc/pedastal",
                         bgo_adc_pedastal, BGO_MAXCHAN);
   err += odbarray_int  ("dragon/bgo/variables/tdc/channel",
                         bgo_tdc_chan,     BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/slope",
                         bgo_adc_slope,    BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/offset",
                         bgo_adc_offset,   BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/posn",
                         bgo_tdc_xposn,    BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/y",
                         bgo_tdc_yposn,    BGO_MAXCHAN);
   err += odbarray_float("dragon/bgo/variables/tdc/z",
                         bgo_tdc_zposn,    BGO_MAXCHAN);
   
   memset(sb_adc_chan,  -1,  SB_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/sb/variables/adc/channel",
                         sb_adc_chan,      SB_MAXCHAN);

   memset(nai_adc_chan,  -1,  NAI_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/nai/variables/adc/channel",
                         nai_adc_chan,      NAI_MAXCHAN);

   tail_xtdc_chan = tail_rftdc_chan = tail_tdc0_chan = -1;
   odbval_int("dragon/tail/variables/xtdc", &tail_xtdc_chan);
   odbval_int("dragon/tail/variables/rftdc", &tail_rftdc_chan);
   odbval_int("dragon/tail/variables/tdc0", &tail_tdc0_chan);

   memset(dssd_adc_chan,  -1,  DSSD_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/dssd/variables/adc/channel",
                         dssd_adc_chan,      DSSD_MAXCHAN);
   err += odbarray_int  ("dragon/dssd/variables/adc/module",
                         dssd_adc_module,    DSSD_MAXCHAN);
   err += odbarray_float("dragon/dssd/variables/adc/slope",
                         dssd_adc_slope,     DSSD_MAXCHAN);
   err += odbarray_float("dragon/dssd/variables/adc/offset",
                         dssd_adc_offset,    DSSD_MAXCHAN);
   dssd_front_tdcchan = dssd_back_tdcchan = -1;
   odbval_int("dragon/dssd/variables/tdc_front", &dssd_front_tdcchan);
   odbval_int("dragon/dssd/variables/tdc_back",  &dssd_back_tdcchan);

   memset(ic_adc_chan,  -1,  IC_MAXCHAN*sizeof(int) );
   memset(ic_tdc_chan,  -1,  IC_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/ic/variables/adc/channel",
                         ic_adc_chan,      IC_MAXCHAN);
   err += odbarray_int  ("dragon/ic/variables/tdc/channel",
                         ic_tdc_chan,      IC_MAXCHAN);

   mcp_tac_chan = mcp_tdc_chan0 = mcp_tdc_chan1 = -1;
   odbval_int("dragon/mcp/variables/tac_adc", &mcp_tac_chan);
   odbval_int("dragon/mcp/variables/tdc", &mcp_tdc_chan0);
   odbval_int("dragon/mcp/variables/tdc", &mcp_tdc_chan1);

   memset(mcp_adc_chan,  -1,  MCP_MAXCHAN*sizeof(int) );
   err += odbarray_int  ("dragon/mcp/variables/adc/channel",
                         mcp_adc_chan,      MCP_MAXCHAN);
   err += odbarray_int  ("dragon/mcp/variables/adc/module",
                         mcp_adc_module,    MCP_MAXCHAN);
   err += odbarray_float("dragon/mcp/variables/adc/slope",
                         mcp_adc_slope,    MCP_MAXCHAN);
   err += odbarray_float("dragon/mcp/variables/adc/offset",
                         mcp_adc_offset,    MCP_MAXCHAN);
   
   if( err ){ printf("read_dragon_odb: %d errors\n", err); }
   return(0);
}

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

   // Initialize all time differences between subsystems to be the default 250ns

   // Intialize all PRESORT timing windows to their defaults
    
   // Search the globals for time difference settings and overwrite their values
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

// (used to be called apply_gains) this is the first function to be called
// on processing an event - before any singles/coinc-sorting ...
// ** the current event has just been added and is last in the window
//      => all other window events are BEFORE the current event

//

int pre_sort_enter(int start_idx, int frag_idx)
{
   Dragon_event *alt, *ptr = &evbuf[frag_idx];
   float energy, ecal, psd, correction;
   int i, ppg_index;
   int dt, bin;

   // Calculate the energy and calibrated energies
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
      if( (count = ptr->v792a.d_count) < MAX_V792_WORDS ){
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
      if( (count = ptr->v1190.d_count) < MAX_V1190_WORDS ){
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
      if( (count = ptr->v792a.d_count) < MAX_V792_WORDS ){
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
      if( (count = ptr->v792b.d_count) < MAX_V792_WORDS ){
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
      if( (count = ptr->v1190.d_count) < MAX_V1190_WORDS ){
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
TH1I  *bgo_e; 

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
   xtofh    = H1_BOOK(cfg, "XTOFH",    "Separator TOF (gamma->HI)", ADC_BINS, 0, ADC_BINS);
   ic_sum   = H1_BOOK(cfg, "IC_SUM",   "Summed Energy Loss in Ion Chamber", ADC_BINS, 0, ADC_BINS);
   bgo_zpat = H1_BOOK(cfg, "BGO_ZPAT", "BGO Z HITPATTERN", ADC_BINS, 0, ADC_BINS);
   bgo_e    = H1_BOOK(cfg, "BGO_E",    "BGO Energy spectrum", ADC_BINS, 0, ADC_BINS);
   ic_0v1   = H2_BOOK(cfg, "IC_0v1",   "Ion Chamber Energy 0 vs 1", ADC_BINS, 0, ADC_BINS, ADC_BINS, 0, ADC_BINS);
   close_folder(cfg);
   close_folder(cfg);

   return(0);
}

// NOTE - have to hunt through event-data arrays for each item
//      - quicker to unpack once into aux arrays (but ~triples memory requirements)
//      - test simple method first
//
// to do proper method efficiently, need reverse tables to quickly sort adc/tdc data
//    into appropriate positions, then just run through adc/tdc data once
int fill_singles_histos(Dragon_event *ptr)
{
   int i, j, chan, val, cnt;
   float v_cal;

   if( ptr->type == HEAD_EVENT ){
      for(j=0; j<SB_MAXCHAN; j++){
         if( (chan = sb_adc_chan[j]) != -1 ){
            if( (cnt = ptr->v792a.d_count) > 0 ){
               for(i=0; i<cnt; i++){
                  val = ptr->v792a.data[i];
                  if( chan != ((val >> 24) & 0xff) ){ continue; }
                  v_cal = (val & 0xffffff);
                  sb_ecal[j]->Fill(sb_ecal[j], (int)v_cal, 1);
               }
            }
         }
      }
   }
   return(0);
}

// loop over window and sort all head-tail pairs
int fill_coinc_histos(int win_idx, int frag_idx)
{
   int dt, abs_dt,  pos, c1, c2, index, ptr_swap;
   return(0);
}

