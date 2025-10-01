//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#######################################################################

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "config.h"
#include "grif-format.h"
#include "histogram.h"
#include "grif-angles.h"
#include "default_sort.h"

int          odb_daqsize;// number of daq channels currently defined in the odb
int         subsys_table[MAX_DAQSIZE];
int        crystal_table[MAX_DAQSIZE]; // Ge/BGO have 4 "crystals" per clover
int        element_table[MAX_DAQSIZE]; // BGO have 5 elements per crystal
int       polarity_table[MAX_DAQSIZE]; // 1 is negative, 0 is positive, -1 is unset
short       address_chan[MAX_ADDRESS];
static short  addr_table[MAX_DAQSIZE]; short   *addrs = addr_table;
char    chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
static int   dtype_table[MAX_DAQSIZE]; int    *dtypes = dtype_table;
static float  gain_table[MAX_DAQSIZE]; float   *gains = gain_table;
static float  offs_table[MAX_DAQSIZE]; float *offsets = offs_table;
static float  quad_table[MAX_DAQSIZE]; float   *quads = quad_table;
float  pileupk1[MAX_DAQSIZE][7];
float  pileupk2[MAX_DAQSIZE][7];
float  pileupE1[MAX_DAQSIZE][7];
float  crosstalk[MAX_DAQSIZE][3][16];
static short *chan_address = addr_table;
static int subsys_initialized[MAX_SUBSYS];
int subsys_deadtime_count[MAX_SUBSYS];
int subsys_prg_ddtm[MAX_SUBSYS];
int previous_trig_acc[MAX_DAQSIZE];
extern Grif_event grif_event[PTR_BUFSIZE];

int presort_window_width= 1940;  // 19.4us needed for all crosstalk corrections. 5us needed for pileup corrections.
int sort_window_width   = 200; //  2us - MAXIMUM (indiv. gates can be smaller)

// Default sort function declarations
extern int init_parameters_from_globals(Config *cfg);
extern int perform_pileup_correction(Grif_event *ptr, Grif_event *alt, int dt, int chan, int chan2, int i, int end_idx);
extern int init_chan_histos(Config *cfg);
extern int init_histos(Config *cfg, int subsystem);
extern int fill_chan_histos(Grif_event *ptr);
extern int fill_singles_histos(Grif_event *ptr);
extern int fill_coinc_histos(int win_idx, int frag_idx);

float spread(int val){ return( val + rand()/(1.0*RAND_MAX) ); }

// odb tables need to be transferred into config, which is saved with histos
int init_default_histos(Config *cfg, Sort_status *arg)
{
  Cal_coeff *cal;
  int i, j;

  // Initialize all pileup and crosstalk parameters to unset values
  for(i=0; i<odb_daqsize; i++){
    for(j=0; j<7; j++){
      pileupk1[i][j] = pileupk2[i][j] = pileupE1[i][j] = -1;
    }
    for(j=0; j<16; j++){
      crosstalk[i][0][j] = crosstalk[i][1][j] = crosstalk[i][2][j] = -1;
    }
  }

  cfg->odb_daqsize = odb_daqsize;
  for(i=0; i<odb_daqsize; i++){ // update config with odb info
    edit_calibration(cfg, chan_name[i], offsets[i], gains[i], quads[i], pileupk1[i], pileupk2[i], pileupE1[i],
      crosstalk[i][0], crosstalk[i][1], crosstalk[i][2], chan_address[i],  dtype_table[i], arg->cal_overwrite );
    }
    // ALSO need to transfer config info to the arrays that are used in sort
    for(i=0; i<odb_daqsize; i++){

      cal = cfg->calib[i];
      if( strcmp(chan_name[i], cal->name) != 0 ){ // conf not in odb order
        for(j=0; j<cfg->ncal; j++){ cal = cfg->calib[j];
          if( strcmp(chan_name[i], cal->name) == 0 ){ break; }
        }
        if( j == cfg->ncal ){ continue; } // not found in config
      }

      // overwrite = 0 => USE CONFIG NOT ODB for offset, gain, quads
      if( arg->cal_overwrite == 0 ){
        offsets[i]=cal->offset; gains[i]=cal->gain;  quads[i]=cal->quad;
      }

      // Pileup or crosstalk parameters do not exist in the MIDAS ODB so must always be copied from the config (file or config)
      if(strncmp(chan_name[i],"GRG",3)==0){
        for(j=0; j<7; j++){
          pileupk1[i][j] = (isnan(cal->pileupk1[j])) ? 0.0 : cal->pileupk1[j];
          pileupk2[i][j] = (isnan(cal->pileupk2[j])) ? 0.0 : cal->pileupk2[j];
          pileupE1[i][j] = (isnan(cal->pileupE1[j])) ? 0.0 : cal->pileupE1[j];
        }
        for(j=0; j<16; j++){
          crosstalk[i][0][j] = (isnan(cal->crosstalk0[j])) ? 0.0 : cal->crosstalk0[j];
          crosstalk[i][1][j] = (isnan(cal->crosstalk1[j])) ? 0.0 : cal->crosstalk1[j];
          crosstalk[i][2][j] = (isnan(cal->crosstalk2[j])) ? 0.0 : cal->crosstalk2[j];
        }
      }
    }

    init_parameters_from_globals(cfg);
    init_chan_histos(cfg);
    init_histos(cfg, SUBSYS_HPGE_A); // always create Ge histos

    // Reset the deadtime counters and previous_trig_acc at BOR
    memset(subsys_deadtime_count,0,MAX_SUBSYS*sizeof(int));
    memset(previous_trig_acc,0,MAX_DAQSIZE*sizeof(int));
    ppg_bin_end = ppg_cycles_binning_factor;

    return(0);
  }

  int init_parameters_from_globals(Config *cfg){
    // Here set parameters from the Globals
    // time-difference conditions for subsys-subsys coincidences
    // PRESORT time-difference Conditions
    // TAC energy offset parameters based on each LBL-LBL combination
    int i,j,k,c1,c2,index;
    Global *global;
    char tmp[32];
    const char *ptr;

    // Initialize all time differences between subsystems to be the default 250ns
    for(i=0; i<MAX_SUBSYS; i++){
      for(j=0; j<MAX_SUBSYS; j++){
        time_diff_gate_min[i][j] = 0;  // default is 0 nanoseconds
        time_diff_gate_max[i][j] = 25; // default is 250 nanoseconds
      }
    }

    // Intialize all PRESORT timing windows to their defaults
    bgo_window_min = addback_window_min = rcmp_fb_window_min = lbl_tac_window_min = zds_tac_window_min = art_tac_window_min = desw_beta_window_min = 0;
    bgo_window_max = 20;
    addback_window_max = 20;
    rcmp_fb_window_max = 10;
    lbl_tac_window_max = 25;
    art_tac_window_max = 25;
    zds_tac_window_max = 25;
    desw_beta_window_max = 50;

    // Initalize the cycles gamma-ray energy gate values
    ppg_cycles_gamma_gate_min = 1800;
    ppg_cycles_gamma_gate_max = 1820;
    ppg_cycles_binning_factor = 10000000; // Default of 10,000,000 converts 10ns to 100 millisecond binning

    // Search the globals for time difference settings and overwrite their values
    for(i=0; i<cfg->nglobal; i++){
      global = cfg->globals[i];
      sprintf(tmp,"%s",global->name);
      if(strncmp(tmp,"time_diff_",10) == 0){
        // This global is a time difference value
        // Identify the subsystem types and then save the value in the correct place
        for(j=0; j<MAX_SUBSYS; j++){
          if(strlen(subsys_handle[j])<2){ continue; }
          if((ptr = strstr(tmp,subsys_handle[j])) > 0){
            ptr += strlen(subsys_handle[j]);
            // Identify the second subsystem type
            for(k=0; k<MAX_SUBSYS; k++){
              if(strlen(subsys_handle[k])<2){ continue; }
              if(strstr(ptr,subsys_handle[k]) > 0){
                // save the value in the correct place
                fprintf(stdout,"time_diff_%s_%s [%d,%d] set to %d,%d\n",subsys_handle[j],subsys_handle[k],j,k,global->min,global->max);
                time_diff_gate_min[j][k] = global->min;
                time_diff_gate_max[j][k] = global->max;
                time_diff_gate_min[k][j] = global->min;
                time_diff_gate_max[k][j] = global->max;
                break;
              }
            } break;
          }
        }
      }else if(strncmp(tmp,"TAC-Offset-",11) == 0){
        // TAC energy offset parameters based on each LBL-LBL combination
        // Only use the maximum value, ignore the minimum
        // Derive the index from the name
        sscanf(tmp,"TAC-Offset-%02d-%02d",&c1,&c2);
        c1--; c2--;
        if(c1>=0 && c1<N_LABR && c2>=0 && c2<N_LABR){
          index = tac_labr_hist_index[c1][c2];
          tac_lbl_combo_offset[index] = global->max;
        }else{
          fprintf(stderr,"Problem decoding TAC Offset from Global, %s\n",tmp);
        }
      }else if(strncmp(tmp,"Timestamp-offset-LBT",20) == 0){
        // LBT timestamp offset value
        // Derive the index from the name
        sscanf(tmp,"Timestamp-offset-LBT%02d",&c1);
        c1--;
        if(c1>=0 && c1<N_TACS){
          tac_ts_offset[c1] = global->max;
        }else{
          fprintf(stderr,"Problem decoding LBT timestamp offset from Global, %s\n",tmp);
        }
      }else if(strncmp(tmp,"presort_time_diff_addback",25) == 0){
        addback_window_min = global->min; addback_window_max = global->max;
      }else if(strncmp(tmp,"presort_time_diff_zds_tac",25) == 0){
        zds_tac_window_min = global->min; zds_tac_window_max = global->max;
      }else if(strncmp(tmp,"presort_time_diff_art_tac",25) == 0){
        art_tac_window_min = global->min; art_tac_window_max = global->max;
      }else if(strncmp(tmp,"presort_time_diff_lbl_tac",25) == 0){
        lbl_tac_window_min = global->min; lbl_tac_window_max = global->max;
      }else if(strncmp(tmp,"presort_time_diff_desw_beta",27) == 0){
        desw_beta_window_min = global->min; desw_beta_window_max = global->max;
      }else if(strncmp(tmp,"presort_time_diff_suppression",29) == 0){
        bgo_window_min = global->min; bgo_window_max = global->max;
      }else if(strncmp(tmp,"presort_time_diff_rcmp_front-back",33) == 0){
        rcmp_fb_window_min = global->min; rcmp_fb_window_max = global->max;
      }else if(strncmp(tmp,"cycles_bin_size_in_ms",21) == 0){
        ppg_cycles_binning_factor = (int)(global->max * 100000); // milliseconds to 10ns timestamp unit conversion
      }else if(strncmp(tmp,"cycles_gamma_gate",17) == 0){
        ppg_cycles_gamma_gate_min = global->min; ppg_cycles_gamma_gate_min = global->max;
      }else if(strncmp(tmp,"presort_window_width",20) == 0){
        presort_window_width = global->max;
      }else if(strncmp(tmp,"sort_window_width",17) == 0){
        sort_window_width = global->max;
      }// end of for(i=0; i<cfg->nglobal; i++){
      }
      return(0);
    }

    //#######################################################################
    //######## PRESORT(gain corrections, addback, suppression)     ##########
    //#######################################################################

    // (used to be called apply_gains) this is the first function to be called
    // on processing an event - before any singles/coinc-sorting ...
    // ** the current event has just been added and is last in the window
    //      => all other window events are BEFORE the current event
    int pre_sort_enter(int start_idx, int frag_idx)
    {
      Grif_event *alt, *ptr = &grif_event[frag_idx];
      int caen_ts_offset = -60; // this value (-60) aligns the timestamps of HPGe with ZDS(CAEN)
      float energy, ecal, psd, correction;
      int i, ppg_index;
      int dt, bin, chan2, chan = ptr->chan;
      int clover, ge1, c1,c2, add;
      int ct_index[4][4] = {{-1,0,1,2},{0,-1,1,2},{0,1,-1,2},{0,1,2,-1}};

      // Protect against invalid channel numbers
      if( chan < 0 || chan >= odb_daqsize ){
        if( ptr->address == 0xFFFF ){
          /*
          ppg_index=-1;
          for(i=0; i<N_PPG_PATTERNS; i++){ if( (ptr->master_pattern & 0xFFFF) == ppg_patterns[i] ){ ppg_index = i; break; } }
          if(ppg_index<0){ fprintf(stderr,"unrecognized ppg pattern, 0x%04X\n", (ptr->master_pattern & 0xFFFF)); return(-1); }
          */
          //  fprintf(stdout,"PPG PATTERN: 0x%04X (%s, %s) at time %10.4f seconds\n", (ptr->master_pattern & 0xFFFF), ppg_handles[ppg_index], ppg_names[ppg_index], (double)(ptr->ts/100000000) );
        } else {
          fprintf(stderr,"unpack_event: ignored event in chan:%d [0x%04x]\n", chan, ptr->address );
        }
        return(-1);
      }

      // Calculate the energy and calibrated energies
      energy = ( ptr->integ1 == 0 ) ? ptr->q1 : spread(ptr->q1)/ptr->integ1;
      ptr->ecal = ptr->esum=offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
      // NOBODY CURRENTLY USES e2,e3,e4 ...

      // Assign the subsys type
      if( (ptr->subsys = subsys_table[chan]) == -1 ){ return(-1); }
      if( subsys_initialized[ptr->subsys] == 0 ){
        //init_histos(configs[1], ptr->subsys);
        init_histos(NULL, ptr->subsys);
      }

      // Check this timestamp against the current cycle bin and fill deadtime histograms if it is a new bin
      if(ppg_cycles_active==1 && ppg_last_ptr_ts>ppg_bin_end && ptr->ts>ppg_last_ptr_ts){
        // Fill deadtime histograms and reset the deadtime count
        bin = (int)((ppg_last_ptr_ts-ppg_cycle_start)/ppg_cycles_binning_factor);  // convert 10ns to binning size set as Global
        if(ppg_cycle_number<MAX_CYCLES){
          gea_cycle_num_dt[ppg_cycle_number]->Fill(gea_cycle_num_dt[ppg_cycle_number], bin, subsys_deadtime_count[SUBSYS_HPGE_A]);
          geb_cycle_num_dt[ppg_cycle_number]->Fill(geb_cycle_num_dt[ppg_cycle_number], bin, subsys_deadtime_count[SUBSYS_HPGE_B]);
          cycle_num_vs_ge_dt->Fill(cycle_num_vs_ge_dt, ppg_cycle_number, bin, subsys_deadtime_count[SUBSYS_HPGE_A]);
          cycle_num_vs_ge_b_dt->Fill(cycle_num_vs_ge_b_dt, ppg_cycle_number, bin, subsys_deadtime_count[SUBSYS_HPGE_B]);
        }
        memset(subsys_deadtime_count,0,MAX_SUBSYS*sizeof(int));
        // Calculate the timestamp of the next bin
        ppg_bin_end += ppg_cycles_binning_factor;
      }

      // Increment deadtime counter for fixed deadtime of this event
      // Check if any events were lost since previous event using Accepted Event counter
      subsys_deadtime_count[ptr->subsys] += (ptr->deadtime>0) ? ptr->deadtime : subsys_prg_ddtm[ptr->subsys];
      if(ptr->trig_acc - previous_trig_acc[chan] != 1 && ptr->trig_acc - previous_trig_acc[chan] != 16383){
        if(ptr->trig_acc - previous_trig_acc[chan]>16100){ // Handle the wrap at 14 bits
          add = ((ptr->trig_acc + 16383) - previous_trig_acc[chan]);
          add *= (ptr->deadtime>0) ? ptr->deadtime : subsys_prg_ddtm[ptr->subsys];
          subsys_deadtime_count[ptr->subsys] += add;
        }else if(ptr->trig_acc - previous_trig_acc[chan] > 0){
          add = (ptr->trig_acc - previous_trig_acc[chan]);
          add *= (ptr->deadtime>0) ? ptr->deadtime : subsys_prg_ddtm[ptr->subsys];
          if(add<2400){ // More than 20 missed events is likely an error
            subsys_deadtime_count[ptr->subsys] += add;
          }
        }
      }
      previous_trig_acc[chan] = ptr->trig_acc;

      // Check this timstamp against the cycle to see if the pattern has changed
      if(ppg_cycles_active==1 && ppg_last_ptr_ts>ppg_pattern_end && ptr->ts>ppg_last_ptr_ts){
        // Recalculate PPG cycle variables
        // Here we update the current PPG pattern, cycle number, cycle start timestamp with the latest values.
        // All subsequent events will use these values

        ppg_pattern_start = ppg_pattern_end;                           // Timestamp of the start of the current pattern
        ppg_cycle_step++;                                              // Current pattern number within this cycle. Patterns counted from zero at beginning of cycle
        if(ppg_cycle_step==ppg_cycle_length){ // Move to next cycle
          ppg_cycle_step = 0;
          ppg_cycle_number++;                                          // Current cycle number. Cycles counted from zero at beginning of run
          ppg_cycle_start = ppg_pattern_start;                       // Timestamp of the start of the current cycle
          ppg_cycle_end += ppg_cycle_duration;                         // Timestamp of the end of the current cycle
        }
        ppg_current_pattern = ppg_cycle_pattern_code[ppg_cycle_step]; // Index of the current PPG cycle pattern for use with the ppg_patterns array
        ppg_pattern_end = ppg_pattern_start + ppg_cycle_pattern_duration[ppg_cycle_step]; // Timestamp of the end of the current pattern
        //  fprintf(stdout,"Cycle %04d, start/finish [%ld/%ld]: step %d, %s, start/finish [%ld/%ld]\n",
        //  ppg_cycle_number, ppg_cycle_start, ppg_cycle_end, ppg_cycle_step, ppg_handles[ppg_current_pattern], ppg_pattern_start, ppg_pattern_end);
      }
      ppg_last_ptr_ts = ptr->ts; // Remember this timestamp for checking at the next event. Avoids rare bug where single events are out of order.

      // The TAC module produces its output signal around 2 microseconds later
      // than the start and stop detector signals are processed.
      if( ptr->subsys == SUBSYS_TAC_LABR || ptr->subsys == SUBSYS_TAC_ZDS || ptr->subsys == SUBSYS_TAC_ART){
        ptr->ts -= tac_ts_offset[crystal_table[ptr->chan]-1]; // Subtract some amount from TAC timestamps
      }

      // DESCANT detectors
      // use psd for Pulse Shape Discrimination provides a distinction between neutron and gamma events
      //if( ptr->subsys == SUBSYS_DESCANT || ptr->subsys == SUBSYS_DESWALL){
      if( ptr->subsys == SUBSYS_DESWALL){
        //ptr->ts -= caen_ts_offset; // Subtract from CAEN timestamps to align coincidences
        psd = ( ptr->q1 != 0 ) ? (spread(ptr->cc_short) / ptr->q1) : 0;
        ptr->psd = (int)(psd*1000.0); // psd = long integration divided by short integration
      }

      // HPGe B
      if( ptr->subsys == SUBSYS_HPGE_B){
        ptr->pu_class = PU_OTHER; // Pileup class - default value for all HPGe events
        if(ptr->pileup==1 && ptr->nhit ==1){
          ptr->pu_class = PU_SINGLE_HIT; // Single hit events, no pileup, this is the most common type of HPGe event
        }
      }

      // HPGe A
      if( ptr->subsys == SUBSYS_HPGE_A){
        ptr->pu_class = PU_OTHER; // Pileup class - default value for all HPGe events
        if(ptr->pileup==1 && ptr->nhit ==1){
          ptr->pu_class = PU_SINGLE_HIT; // Single hit events, no pileup, this is the most common type of HPGe event
        }

        // HPGe Clover time-dependant crosstalk corrections within same clover
        i = start_idx;
        while( i != frag_idx ){ // need at least two events in window
          if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
          if( (chan2 = alt->chan)<0 || alt->chan >= odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
            continue;
          }
          if((dt=ptr->ts - alt->ts)>479 || alt->subsys != SUBSYS_HPGE_A){ continue; }

          if(chan2 != chan ){
            if((clover=(int)(crystal_table[chan2]/4)) == (int)(crystal_table[chan]/4)){
              // HPGe Clover time-dependant crosstalk corrections within same clover
              // dt is always positive here
              // The original hit (ptr) came after the crosstalk-inducing hit (alt)
              // Make correction to ptr hit based on energy of alt.
              bin = (int)((1940+dt)/160);
              if(bin<0 || bin>15){ fprintf(stderr,"pre_sort_enter bin [%d] out of bounds for dt %d\n",bin,dt); continue; }
              ge1 = crystal_table[chan];
              c1 = ge1%4;
              c2 = ct_index[c1][(int)(crystal_table[chan2]%4)];
              if(crosstalk[ge1][c2][bin] != -1 ){
                //  correction = crosstalk[ge1][c2][bin] + ((crosstalk[ge1][c2][bin+1] - crosstalk[ge1][c2][bin]) * (float)(((1940+dt)%160)/160));
                correction = crosstalk[ge1][c2][bin];
                //  fprintf(stdout,"CT enter, %d %d: %d %f %f %f %f: %f + %f = %f\n",chan,chan2,bin,crosstalk[ge1][c2][bin+1],crosstalk[ge1][c2][bin],(float)(((1940+dt)%160)/160),alt->ecal,ptr->ecal,(alt->ecal * correction),(ptr->ecal+(alt->ecal * correction)));
                ptr->ecal += alt->ecal * correction;
              }
            }
          }
        } // end of while


        // Fill crosstalk histograms
        i = start_idx;
        while( i != frag_idx ){ // need at least two events in window
          if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
          if( (chan2 = alt->chan)<0 || alt->chan >= odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
            continue;
          }
          if((dt=ptr->ts - alt->ts)>1000 || alt->subsys != SUBSYS_HPGE_A){ continue; }

          if(alt->ecal>1327 && alt->ecal<1337){ // Crosstalk inducing hit was 1332keV
            if(chan2 != chan ){
              if((clover=(int)(crystal_table[chan2]/4)) == (int)(crystal_table[chan]/4)){

                // (c2%4) = Crystal Color [B, G, R, W]
                // Hits with ptr arriving after alt
                switch(crystal_table[chan2]%4){
                  case 0:  ct_e_vs_dt_B[crystal_table[chan]]->Fill(ct_e_vs_dt_B[crystal_table[chan]], ((int)((alt->ts - ptr->ts)/4)+300), (int)ptr->ecal-1100, 1); break;
                  case 1:  ct_e_vs_dt_G[crystal_table[chan]]->Fill(ct_e_vs_dt_G[crystal_table[chan]], ((int)((alt->ts - ptr->ts)/4)+300), (int)ptr->ecal-1100, 1); break;
                  case 2:  ct_e_vs_dt_R[crystal_table[chan]]->Fill(ct_e_vs_dt_R[crystal_table[chan]], ((int)((alt->ts - ptr->ts)/4)+300), (int)ptr->ecal-1100, 1); break;
                  case 3:  ct_e_vs_dt_W[crystal_table[chan]]->Fill(ct_e_vs_dt_W[crystal_table[chan]], ((int)((alt->ts - ptr->ts)/4)+300), (int)ptr->ecal-1100, 1); break;
                }
              }
            }
          }
        } // end of while


      } // end of if( ptr->subsys == SUBSYS_HPGE_A){
        return(0);
      }

      // Presort - do Suppression and Addback here
      //  - frag_idx is about to leave coinc window (which ends at end_idx)
      //    check other frags in window for possible suppression and/or summing
      //  also calculate multiplicities[store in frag_idx only]
      int pre_sort_exit(int frag_idx, int end_idx)
      {
        Grif_event *alt2, *alt, *ptr = &grif_event[frag_idx];
        float desw_median_distance = 1681.8328; // descant wall median source-to-detector distance in mm
        int i, j, dt, dt13, tof;
        float q1,integ2,q12,k1,k2,k12,e1,e2,e12,m,c;
        int chan,chan2,found,pos;
        int clover, ge1, c1,c2,bin;
        float energy,ecal,correction;
        int ct_index[4][4] = {{-1,0,1,2},{0,-1,1,2},{0,1,-1,2},{0,1,2,-1}};

        // Assign chan local variable and check it is a valid channel number
        if( (chan=ptr->chan)<0 || ptr->chan >= odb_daqsize ){
          fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
          return(-1);
        }
        i = frag_idx; ptr->multiplicity = 1;
        while( i != end_idx ){ // need at least two events in window
          if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
          if( (chan2=alt->chan)<0 || alt->chan >= odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
            continue;
          }

          // Determine absolute time difference between timestamps
          dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }

          // Restrict to 2 microseconds presort window for everything except HPGe crosstalk to maintain speed
          if(alt->subsys != SUBSYS_HPGE_A && dt>250){ continue; }

          // Determine multiplicity
          if( alt->subsys == ptr->subsys ){ ++ptr->multiplicity; }

          // SubSystem-specific pre-processing
          switch(ptr->subsys){
            case SUBSYS_HPGE_A:

            // HPGe Clover time-dependant crosstalk corrections
            if(alt->subsys == SUBSYS_HPGE_A && chan2 != chan){
              if((clover=(int)(crystal_table[chan2]/4)) == (int)(crystal_table[chan]/4)){
                dt = ptr->ts - alt->ts; // Use relative time difference. This is always negative

                if(dt>-1940 && dt <= 0){
                  // The original hit (ptr) came earlier then the crosstalk-inducing hit (alt)
                  // Make correction to ptr hit based on energy of alt.
                  //  dt -= 1920; // Correct the timestamp difference so that the right correction is calculated
                  bin = (int)((1940+dt)/160);
                  if(bin<0 || bin>15){ fprintf(stderr,"pre_sort_exit bin [%d] out of bounds for dt %d\n",bin,dt); continue; }
                  ge1 = crystal_table[chan];
                  c1 = ge1%4;
                  c2 = ct_index[c1][(int)(crystal_table[chan2]%4)];
                  if(crosstalk[ge1][c2][bin] != -1 ){
                    //correction = crosstalk[ge1][c2][bin] + ((crosstalk[ge1][c2][bin+1] - crosstalk[ge1][c2][bin]) * (((1940+dt)%160)/160));
                    correction = crosstalk[ge1][c2][bin];
                    //  fprintf(stdout,"CT exit, %d %d: %d %f %f %f %f: %f + %f = %f\n",chan,chan2,bin,crosstalk[ge1][c2][bin+1],crosstalk[ge1][c2][bin],(float)(((1940+dt)%160)/160),alt->ecal,ptr->ecal,(alt->ecal * correction),(ptr->ecal+(alt->ecal * correction)));
                    ptr->ecal += alt->ecal * correction;
                  }
                }
                if( dt < 0 ){ dt = -1*dt; } // Reset the abs time difference for anything that follows
              }
            }

            // HPGe pile-up corrections
            // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
            // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
            // First assign the pileup class type, then correct the energies
            if(alt->subsys == SUBSYS_HPGE_A && chan2 == chan){
              perform_pileup_correction(ptr, alt, dt, chan, chan2, i, end_idx);
            }

            // BGO suppression of HPGe
            if( (dt >= bgo_window_min && dt <= bgo_window_max) && alt->subsys == SUBSYS_BGO && !ptr->suppress ){
              // could alternatively use crystal numbers rather than clover#
              //    (don't currently have this for BGO)
              if( crystal_table[ptr->chan]/16 == crystal_table[alt->chan]/16 ){ ptr->suppress = 1; }
            }
            // Germanium addback -
            //    earliest fragment has the sum energy, others are marked -1
            // Remember the other crystal channel number in alt_chan for use in Compton Polarimetry
            if( (dt >= addback_window_min && dt <= addback_window_max) && alt->subsys == SUBSYS_HPGE_A ){
              if( alt->esum >= 0 && crystal_table[alt->chan]/16 == crystal_table[ptr->chan]/16 ){
                ptr->esum += alt->esum; alt->esum = -1; ptr->alt_chan = alt->chan;
              }
            }
            break;
            case SUBSYS_HPGE_B:
            // HPGe B pile-up corrections
            if(alt->subsys == SUBSYS_HPGE_B && chan2 == chan){
              perform_pileup_correction(ptr, alt, dt, chan, chan2, i, end_idx);
            }
            break;
            case SUBSYS_RCMP:
            // RCMP Front-Back coincidence
            // Ensure its the same DSSD and that the two events are front and back
            // The charged particles enter the P side and this has superior energy resolution
            // Ensure the energy collected in the front and back is similar
            ptr->esum = -1; // Need to exclude any noise and random coincidences.
            if( alt->subsys == SUBSYS_RCMP && (dt >= rcmp_fb_window_min && dt <= rcmp_fb_window_max) && (ptr->ecal>0 && ptr->ecal<32768)){
              if((crystal_table[ptr->chan] == crystal_table[alt->chan]) && (polarity_table[ptr->chan] != polarity_table[alt->chan]) && (alt->ecal > 0  && ptr->ecal<32768)){
                if( ((ptr->ecal / alt->ecal)<=1.1 && (ptr->ecal / alt->ecal)>=0.9)){
                  // Ensure esum comes from P side, but use this timestamp
                  ptr->esum = polarity_table[ptr->chan]==0 ? ptr->ecal : (polarity_table[ptr->chan]==1 ? alt->ecal : -1);
                  ptr->suppress = alt->suppress = 1;
                }
              }
            }
            break;
            case SUBSYS_TAC_ZDS:
            // ZDS TAC spectra
            if( alt->subsys == SUBSYS_ZDS_A ){
              // For TAC08 the start is ZDS and the stop is any of the LaBr3. So this is three detector types.
              // Here in the presort we will remember the ZDS information that is in coincidence with the TAC.
              // In the TAC event we save the ZDS chan as alt_chan, and the ZDS energy as alt_ecal.
              // So later in the main coincidence loop we only need to examine LBL and TAC.
              if( dt >= zds_tac_window_min && dt <= zds_tac_window_max ){
                ptr->alt_chan = alt->chan; ptr->alt_ecal = alt->ecal;
              }
            }
            break;
            case SUBSYS_TAC_ART:
            // ARIES TAC spectra
            if( alt->subsys == SUBSYS_ARIES_A ){
              // For TAC08 the start is ARIES and the stop is any of the LaBr3. So this is three detector types.
              // Here in the presort we will remember the ARIES tile that is in coincidence with the TAC.
              // In the TAC event we save the tile chan as alt_chan, and the tile energy as alt_ecal.
              // So later in the main coincidence loop we only need to examine LBL and TAC.
              if( dt >= art_tac_window_min && dt <= art_tac_window_max ){
                ptr->alt_chan = alt->chan; ptr->alt_ecal = alt->ecal;
              }
            }
            break;
            case SUBSYS_LABR_L:
            // LaBr3 TAC spectra
            if( alt->subsys == SUBSYS_TAC_LABR ){
              // For TAC01-07 we have a LBL-LBL coincidence
              // Here save the LBL Id number and the LBL energy in the TAC event
              // Save LBL channel number into ptr->integ2 or integ3 or integ4
              // Save LBL energy ecal into TAC ptr-q2 or q3 or q4
              if( dt >= lbl_tac_window_min && dt <= lbl_tac_window_max ){
                if( alt->q2 < 1 ){ // First LBL in coincidence with this TAC
                  alt->integ2 = ptr->chan; alt->q2 = ptr->ecal;
                } else if( alt->q3 < 1 ){ // This is the second LBL in coincidence with this TAC
                  if( alt->integ2 < 0 || alt->integ2 >= odb_daqsize ){ break; }
                  if( crystal_table[ptr->chan] < crystal_table[alt->integ2] ){ // Order the LBL by crystal number not timestamp
                    alt->integ3 = alt->integ2;   alt->q3 = alt->q2;
                    alt->integ2 = ptr->chan; alt->q2 = ptr->ecal;
                  } else { // More than two LBL in coincidence with this TAC
                    alt->integ3 = ptr->chan; alt->q3 = ptr->ecal;
                  }
                } else {
                  alt->integ4 = ptr->chan; alt->q4 = ptr->ecal; // If this is set then we have LBL multiplicity >2 for this TAC
                }
              }
            }
            break;
            case SUBSYS_ZDS_B: // CAEN Zds
            if(alt->subsys == SUBSYS_DESWALL){
              if(dt >= desw_beta_window_min && dt <= desw_beta_window_max){
                // Calculate time-of-flight and correct it for this DESCANT detector distance
                tof = (spread(abs(ptr->cfd - alt->cfd))*2.0) + 100; //if( tof < 0 ){ tof = -1*tof; }
                //  fprintf(stdout,"tof: %d - %d = %f\n",ptr->cfd, alt->cfd, tof);
                alt->tof = (int)(tof); // Time of flight
                alt->alt_ecal = (int)(spread(tof) * DSW_tof_corr_factor[crystal_table[alt->chan]-1]); // Corrected Time of Flight
                desw_psd_zdse->Fill(desw_psd_zdse, (int)alt->psd, (int)ptr->ecal, 1); // Fill desw_psd_zdse
              }
            }
            break;
            default: break; // Unrecognized or unprocessed subsys type
          }// end of switch
        }// end of while


        // Fill crosstalk histograms
        if(ptr->subsys == SUBSYS_HPGE_A){
          i = frag_idx;
          while( i != end_idx ){ // need at least two events in window
            if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
            if( (chan2=alt->chan)<0 || alt->chan >= odb_daqsize ){
              fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
              continue;
            }
            if(alt->subsys != SUBSYS_HPGE_A){ continue; }

            if(ptr->ecal>1327 && ptr->ecal<1337){ // Crosstalk inducing hit was 1332keV
              if(chan2 != chan ){
                if((clover=(int)(crystal_table[chan2]/4)) == (int)(crystal_table[chan]/4)){

                  // Fill crosstalk histograms
                  // (c2%4) = Crystal Color [B, G, R, W]
                  // Hits with ptr arriving after alt
                  // Fill higher-channels of the x axis on these plots
                  switch(crystal_table[chan]%4){
                    case 0:  ct_e_vs_dt_B[crystal_table[chan2]]->Fill(ct_e_vs_dt_B[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)/4)+300), (int)alt->ecal-1100, 1); break;
                    case 1:  ct_e_vs_dt_G[crystal_table[chan2]]->Fill(ct_e_vs_dt_G[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)/4)+300), (int)alt->ecal-1100, 1); break;
                    case 2:  ct_e_vs_dt_R[crystal_table[chan2]]->Fill(ct_e_vs_dt_R[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)/4)+300), (int)alt->ecal-1100, 1); break;
                    case 3:  ct_e_vs_dt_W[crystal_table[chan2]]->Fill(ct_e_vs_dt_W[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)/4)+300), (int)alt->ecal-1100, 1); break;
                  }
                }
              }
            }
          }// end of while
        }

        return(0);
      }

      // HPGe pile-up corrections
      // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
      // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
      // First assign the pileup class type, then correct the energies
      int perform_pileup_correction(Grif_event *ptr, Grif_event *alt, int dt, int chan, int chan2, int i, int end_idx)
      {
        Grif_event *alt2;
        int j, dt13;
        float k1,k2, energy, correction, correction12, correction23;

        if(ptr->pileup==1 && ptr->nhit ==1){
          ptr->pu_class = PU_SINGLE_HIT; // no pileup, this is the most common type of HPGe event
          return(0);
        }else if(ptr->pileup==0){
          ptr->pu_class = PU_ERROR; // Pileup class, error
          return(0);
        }
        if(dt>500){ return(0); } // Restrict pileup handling to 5 microseconds. 8 microseconds needed in some early datasets

        // Two hit pileup...
        if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==2 && alt->nhit==1)){
          ptr->pu_class = alt->pu_class = PU_2HIT_ERROR; // Pileup class, error for 2Hits
          ptr->delta_t = alt->delta_t = dt;
          if(ptr->q1>0 && ptr->integ1>0 && ptr->q2>0 && ptr->integ2>0 && alt->q1>0 && alt->integ1>0){
            // 2 Hit, Type A. The (ptr) fragement is the first Hit of a two Hit pile-up event.
            ptr->pu_class = PU_2HIT_A1ST; alt->pu_class = PU_2HIT_A2ND; // Pileup class, 1st and 2nd of 2Hits
            ptr->delta_t = alt->delta_t = dt;  // time difference between hits
          }else{
            // 2 Hit, Type B, where 2nd Hit integration region starts after 1st Hit integration ends but before 2nd Hit CFD has completed
            ptr->pu_class = PU_2HIT_C1ST; alt->pu_class = PU_2HIT_C2ND; // Pileup class, 1st and 2nd of 2Hits
            ptr->delta_t = alt->delta_t = dt;  // time difference between hits
          }
        }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==1 && alt->nhit==1)){
          // 2 Hit, Type C, where 2nd Hit integration region starts after 1st Hit integration and 2nd Hit CFD have ended
          ptr->pu_class = PU_2HIT_B1ST; alt->pu_class = PU_2HIT_B2ND; // Pileup class, 1st and 2nd  of 2Hits
          ptr->delta_t = alt->delta_t = dt; // Save the time difference between pileup hits into both hits
        }
        // Three hit pileup...
        else if((ptr->pileup==1 && ptr->nhit==3) && (alt->pileup==2 && alt->nhit==2)){ // 3Hit pileup
          ptr->pu_class = alt->pu_class = PU_3HIT_ERROR; // Pileup class, error for 3Hits
          if(ptr->q1>0 && ptr->integ1>0 && ptr->q2>0 && ptr->integ2>0 && alt->q1>1 && alt->integ1>0 && alt->q2>0 && alt->integ2>0){
            j=i+1;
            while( j != end_idx ){ // need to find the third events in window associated with this channel
              if( ++j >=  PTR_BUFSIZE ){ break; } alt2 = &grif_event[j]; // WRAP
              if(alt2->chan == chan){ // It must also be a HPGe if the channel number is the same
                if(alt2->pileup==3 && alt2->nhit==1){
                  alt2->pu_class = PU_3HIT_3RD; // Pileup class
                  if(alt2->q1>1 && alt2->integ1>0){
                    // Determine absolute time difference between timestamps for Hit 1 and 3
                    dt13 = ptr->ts - alt2->ts; if( dt13 < 0 ){ dt13 = -1*dt13; }
                    // The Differencitation period of the HPGe Pulse Height evaluation is L = 5000ns.
                    if(dt13>500){
                      // Triple pileup case A ... in which the 3rd pulse occurs more than L samples after the first
                      //                          there are 5 regions, 3 of which are not piled up (1 per pulse)
                      correction23 = (alt->q1/alt->integ1)-((alt->q2/alt->integ2)-(alt2->q1/alt2->integ1));
                      correction12 = (ptr->q1/ptr->integ1)-((ptr->q2/ptr->integ2)-(alt->q1/alt->integ1)-correction23);
                      // Hit 1
                      ptr->pu_class = PU_3HIT_1ST;
                      energy = (spread(ptr->q1)/ptr->integ1) + correction12;
                      ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
                      // Hit 2
                      alt->delta_t = dt; alt->pu_class = PU_3HIT_2ND;
                      energy = (spread(alt->q1)/alt->integ1) - correction12 + correction23;
                      alt->ecal=alt->esum = offsets[chan2]+energy*(gains[chan2]+energy*quads[chan2]);
                      // Hit 3
                      alt2->delta_t = dt13; alt2->pu_class = PU_3HIT_3RD;
                      energy = (spread(alt2->q1)/alt2->integ1) - correction23;
                      alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
                    }else{
                      // Triple pileup case B ... in which the 3rd pulse occurs less than L samples after the first
                      //                          again 5 regions, only 2 of which are not piled up (first and last pulse)
                      //                          There is no region to obtain the height of pulse 2
                      //                          so the event contains K12, the sum of pulse 1+2, in place of pulseheight2
                      correction23 = (alt->q1/alt->integ1)-((alt->q2/alt->integ2)-(alt2->q1/alt2->integ1));
                      correction12 = (ptr->q1/ptr->integ1)-((ptr->q2/ptr->integ2)-(alt->q1/alt->integ1)-correction23);
                      // Hit 1
                      ptr->pu_class = PU_3HIT_1ST;
                      energy = (spread(ptr->q1)/ptr->integ1) + correction12;
                      ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
                      // Hit 2
                      alt->delta_t = dt; alt->pu_class = PU_3HIT_2ND;
                      energy = (spread(alt->q1)/alt->integ1) - correction12 + correction23;
                      alt->ecal=alt->esum = offsets[chan2]+energy*(gains[chan2]+energy*quads[chan2]);
                      // Hit 3
                      alt2->delta_t = dt13; alt2->pu_class = PU_3HIT_3RD;
                      energy = (spread(alt2->q1)/alt2->integ1) - correction23;
                      alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
                    }
                    break; // Break the while if we found the third Hit
                  }
                }
              }
            } // end of while for triple coincidence
          }
        } // end of 3Hit pileup type assignments

        // Now apply hit-specific energy corrections
        if(pileupk1[chan][0] != 1){
          if(ptr->pu_class>=PU_2HIT_A1ST && ptr->pu_class<=PU_2HIT_C2ND){ // 2-Hit pileup events
            // Apply the k1 dependant correction to the energy of the first hit
            // It was already checked that chan for ptr and alt are the same for pileup events
            k1 = ptr->integ1;
            ptr->ecal=ptr->esum = ptr->ecal*( pileupk1[chan][0]+(k1*pileupk1[chan][1])+(k1*k1*pileupk1[chan][2])+(k1*k1*k1*pileupk1[chan][3])
            +(k1*k1*k1*k1*pileupk1[chan][4])+(k1*k1*k1*k1*k1*pileupk1[chan][5])+(k1*k1*k1*k1*k1*k1*pileupk1[chan][6]));

            // Apply the E1 and k2 dependant offset correction to the energy of the second hit
            // Apply the k2 dependant correction to the energy of the second hit
            k2 = alt->integ1;
            correction = ptr->ecal*( pileupE1[chan][0]+(k2*pileupE1[chan][1])+(k2*k2*pileupE1[chan][2])+(k2*k2*k2*pileupE1[chan][3])
            +(k2*k2*k2*k2*pileupE1[chan][4])+(k2*k2*k2*k2*k2*pileupE1[chan][5])+(k2*k2*k2*k2*k2*k2*pileupE1[chan][6]));
            alt->ecal=alt->esum = (alt->ecal*( pileupk2[chan][0]+(k2*pileupk2[chan][1])+(k2*k2*pileupk2[chan][2])+(k2*k2*k2*pileupk2[chan][3])
            +(k2*k2*k2*k2*pileupk2[chan][4])+(k2*k2*k2*k2*k2*pileupk2[chan][5])+(k2*k2*k2*k2*k2*k2*pileupk2[chan][6])))+correction;
          }
        }
        alt->alt_ecal=ptr->ecal; // Remember the ecal of the first Hit in this second Hit. Must be done regardless if a correction is made

        return(0);
      }

      int default_sort(int win_idx, int frag_idx, int flag)
      {
        Grif_event *ptr;
        int i;

        // sort first event, even if window only contains that event
        // (usually require at least two)
        for(i=win_idx; ; i++){ ptr = &grif_event[i];
          if( i >= PTR_BUFSIZE ){ i=0; } // WRAP
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

      //#######################################################################
      //########        Individual channel singles HISTOGRAMS        ##########
      //#######################################################################

      int init_chan_histos(Config *cfg)
      {                      // 1d histograms for Q,E,T,Wf for each channel in odb
        char title[STRING_LEN], handle[STRING_LEN];
        int i, j, k, pos;

        open_folder(cfg, "Hits_and_Sums");
        open_folder(cfg, "Hits");
        for(i=0; i<N_HITPAT; i++){ // Create Hitpattern spectra
          sprintf(title,  "Hitpattern_%s",    hit_names[i] );
          hit_hist[i] = H1_BOOK(cfg, hit_handles[i], title, MAX_DAQSIZE, 0, MAX_DAQSIZE);
        }
        ts_hist = H1_BOOK(cfg, "ts", "Timestamp", 163840, 0, 163840);
        gc_hist = H1_BOOK(cfg, "gc", "ZDS GRIF-CAEN", 16, 0, 16);
        close_folder(cfg);
        open_folder(cfg, "Multiplicities");
        for(i=0; i<MAX_SUBSYS; i++){ mult_hist[i] = NULL;
          if( strncmp(subsys_handle[i],"XXX",3) == 0 ||
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
            if( strncmp(subsys_handle[j],"XXX",3) == 0 ||
            strlen(subsys_handle[j]) == 0 ){ continue; }
            if( memcmp(subsys_handle[j], chan_name[i], 3) == 0 ){ break; }
          }                              // stop at final entry "unknown"
          if(j==(MAX_SUBSYS-1)){ if(strncmp(chan_name[i],"GRS",3) == 0){ j=8; } } // GRS != BGO
          open_folder(cfg, subsys_name[j]);
          open_folder(cfg, "Energy");
          sprintf(title,  "%s_Energy",         chan_name[i] );
          sprintf(handle, "%s_E",              chan_name[i] );
          if( strcmp(subsys_handle[j],"LBT") == 0){
            e_hist[i] = H1_BOOK(cfg, handle, title, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
          }else{
            e_hist[i] = H1_BOOK(cfg, handle, title, E_SPECLEN, 0, E_SPECLEN);
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
            ph_hist[i] = H1_BOOK(cfg, handle, title, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
          }else{
            ph_hist[i] = H1_BOOK(cfg, handle, title, E_SPECLEN, 0, E_SPECLEN);
          }
          close_folder(cfg);
          if( strcmp(subsys_handle[j],"DSW") == 0){
            pos  = crystal_table[i];
            if(pos>0 && pos<=N_DES_WALL){
              open_folder(cfg, "PSD");
              sprintf(title,  "%s_PSD",         chan_name[i] );
              sprintf(handle, "%s_PSD",         chan_name[i] );
              desw_psd[pos] = H1_BOOK(cfg, handle, title, E_PSD_SPEC_LENGTH, 0, E_PSD_SPEC_LENGTH);
              close_folder(cfg);
              open_folder(cfg, "Time_Of_Flight");
              sprintf(title,  "%s_TOF_PSD-gated", chan_name[i] );
              sprintf(handle, "%s_CTOF_PSDn",         chan_name[i] );
              desw_tof_psd[pos] = H1_BOOK(cfg, handle, title, E_TOF_SPEC_LENGTH, 0, E_TOF_SPEC_LENGTH);
              sprintf(title,  "%s_Corrected_TOF", chan_name[i] );
              sprintf(handle, "%s_CTOF",         chan_name[i] );
              desw_tof_corr[pos] = H1_BOOK(cfg, handle, title, E_TOF_SPEC_LENGTH, 0, E_TOF_SPEC_LENGTH);
              sprintf(title,  "%s_TOF",         chan_name[i] );
              sprintf(handle, "%s_TOF",         chan_name[i] );
              desw_tof[pos] = H1_BOOK(cfg, handle, title, E_TOF_SPEC_LENGTH, 0, E_TOF_SPEC_LENGTH);
              close_folder(cfg);
            }
          }
          close_folder(cfg);
          if( strncmp(subsys_handle[j],"XXX",3) == 0 ){     // suppress XXX channels
            ph_hist[i]->suppress = e_hist[i]->suppress =
            /*  cfd_hist[i]->suppress = wave_hist[i]->suppress = */ 1;
          }
        }
        return(0);
      }

      int fill_chan_histos(Grif_event *ptr)
      {
        static int event;
        int chan, sys, pos;

        // Check for unassigned channel numbers
        if( (chan = ptr->chan) == -1 ){
          return(-1);
        }

        // Check for invalid channel numbers, prossibly due to data corruption
        if( chan < 0 || chan > odb_daqsize ){
          fprintf(stderr,"Invalid channel number in fill_chan_histos(), %d\n",chan);
          return(-1);
        }
        if( ++event < 16384 ){
          ts_hist -> Fill(ts_hist, event,  (int)(ptr->ts/100));
        } else if( (event % 1000) == 0 ){
          ts_hist -> Fill(ts_hist, 16367+(int)(event/1000),  (int)(ptr->ts/100));
        }
        ph_hist[chan] -> Fill(ph_hist[chan],  (int)(( ptr->integ1 == 0 ) ? ptr->q1 :
        spread(ptr->q1)/ptr->integ1),  1);
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
        //if( ptr->wf_present  != 0 ){ hit_hist[3] -> Fill(hit_hist[3], chan, 1); }
        //if( ptr->scl_present != 0 ){ hit_hist[4] -> Fill(hit_hist[4], chan, 1); }

        return(0);
      }

      //#######################################################################
      //########               Sums and coinc  HISTOGRAMS            ##########
      //#######################################################################

      // Only create histograms when that subsystem is seen in the data
      // Do not create unused histos -> need to store subsystem types along with definitions ...
      // NOTE: only store one subsystem (even for 2d) - choose the most likely to save space/time
      typedef struct histo_def_struct {
        void **ptr; char title[80]; void *handle;
        int subsys; int xchan; int ychan; int count;
      } Histogram_definition;
      // histo folders have ptr=NULL and xchan=0
      // count != 0 is an array of histos
      //    either - array title will contain %d (padding such as %02d is allowed). The index will run from 0 to (count-1)
      //    or     - title is "" (a blank string) and handle is pointer to array of handles.
      // NOTE that if %d or an array of handles are used then the **ptr should be a pointer and not an address (ie. do not include the leading &)
      // --------------------
      // Set ychan to "SYMMETERIZE" to create a fully symmeterized 2d histogram
      // --------------------
      // default handles derived from title with the following substitutions
      //    "Energy"->"E", "CrystalNum"->"Xtal", "Des_Wall"->"DESW"
      //    "Upstream"->"US", "Downstream"->"DS", ""->"", ""->"",
      #define HISTO_DEF_SIZE 256
      Histogram_definition histodef_array[HISTO_DEF_SIZE] = {
        // Singles
        {NULL,                   "Hits_and_Sums/Sums",                                   },
        {(void **)&ge_sum_ab,    "Addback_Sum_Energy",      "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum,       "Ge_Sum_Energy",           "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b,     "Ge_Sum_En_betaTagged",    "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_sep, "Ge_Sum_En_SceptarTagged", "Ge_Sum_E_B_SEP",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_zds, "Ge_Sum_En_ZdsTagged",     "Ge_Sum_E_B_ZDS",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_art, "Ge_Sum_En_AriesTagged",   "Ge_Sum_E_B_ART",            SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_b_art_brems, "Ge_Sum_En_AriesTagged_Brems",   "Ge_Sum_E_B_ART_Brems",SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&paces_sum,    "PACES_Sum_Energy",        "",                          SUBSYS_PACES,   E_SPECLEN},
        {(void **)&labr_sum,     "LaBr3_Sum_Energy",        "",                          SUBSYS_LABR_L,  E_SPECLEN},
        {(void **)&aries_sum,    "ARIES_Sum_Energy",        "",                          SUBSYS_ARIES_A, E_SPECLEN},
        {(void **)&rcmp_sum,     "RCMP_Sum_Energy",         "",                          SUBSYS_RCMP,    E_SPECLEN},
        {(void **)&rcmp_fb_sum,  "RCMP_Sum_FB_Energy",      "",                          SUBSYS_RCMP,    E_SPECLEN},
        {(void **)&desw_sum_e,   "DES_Wall_Sum_Energy",     "",                          SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof, "DES_Wall_Sum_TOF",        "",                          SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_psd, "DES_Wall_Sum_PSD",        "",                          SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_e_b,      "DES_Wall_Sum_En_betaTagged",  "DESW_Sum_E_B",     SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof_b,    "DES_Wall_Sum_TOF_betaTagged", "DESW_Sum_TOF_B",   SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_e_nn,     "DES_Wall_Sum_En_fold2",       "DESW_Sum_E_nn",    SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof_nn,   "DES_Wall_Sum_TOF_fold2",      "DESW_Sum_TOF_nn",  SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_e_nn_a,   "DES_Wall_Sum_En_fold2_ang60", "DESW_Sum_E_nn_a",  SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&desw_sum_tof_nn_a, "DES_Wall_Sum_TOF_fold2_ang60","DESW_Sum_TOF_nn_a",SUBSYS_DESWALL, E_SPECLEN},
        {(void **)&ge_sum_us,    "Upstream_Ge_Sum_Energy",  "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ds,    "Downstream_Ge_Sum_Energy","",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ab_us, "Upstream_AB_Sum_Energy",  "",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **)&ge_sum_ab_ds, "Downstream_AB_Sum_Energy","",                          SUBSYS_HPGE_A,  E_SPECLEN},
        {(void **) ge_ab_e,      "Addback_%d",               "",                         SUBSYS_HPGE_A,  E_SPECLEN, 0, N_CLOVER},
        {NULL,                   "Hits_and_Sums/Energy",     "",                         },
        {(void **)&ge_xtal,      "GeEnergy_CrystalNum",      "",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&ge_xtal_1hit, "GeEnergy_CrystalNum_1_hit","",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&ge_xtal_2hit, "GeEnergy_CrystalNum_2_hit","",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&ge_xtal_3hit, "GeEnergy_CrystalNum_3_hit","",                         SUBSYS_HPGE_A,   64, E_2D_SPECLEN},
        {(void **)&bgo_xtal,     "BgoEnergy_CrystalNum",     "",                         SUBSYS_BGO,     320, E_2D_SPECLEN},
        {(void **)&bgof_xtal,    "BgoFrontEnergy_CrystalNum","",                         SUBSYS_BGO,     128, E_2D_SPECLEN},
        {(void **)&bgos_xtal,    "BgoSideEnergy_CrystalNum", "",                         SUBSYS_BGO,     128, E_2D_SPECLEN},
        {(void **)&bgob_xtal,    "BgoBackEnergy_CrystalNum", "",                         SUBSYS_BGO,      64, E_2D_SPECLEN},
        {(void **)&bgoa_xtal,    "BgoAncilEnergy_CrystalNum","",                         SUBSYS_LABR_BGO, 32, E_2D_SPECLEN},
        {(void **)&labr_xtal,    "Labr3Energy_CrystalNum",   "LabrE_Xtal",               SUBSYS_LABR_L,   16, E_2D_SPECLEN},
        {(void **)&paces_xtal,   "PacesEnergy_CrystalNum",   "PacesE_Xtal",              SUBSYS_PACES,    16, E_2D_SPECLEN},
        {(void **)&aries_xtal,   "AriesEnergy_CrystalNum",   "AriesE_Xtal",              SUBSYS_ARIES_A,  80, E_2D_SPECLEN},
        // {(void **)&labr_tac_xtal,"TAC_LBL_ART_vs_LBL_Num", "TAC_ART_LBL_LBL_Xtal",       SUBSYS_TAC_ART,   16,  E_2D_SPECLEN},
        {(void **)&art_tac_xtal, "TAC_LBL_ART_vs_ART_Num", "TAC_ART_LBL_ART_Xtal",       SUBSYS_ARIES_A,  80,  E_2D_SPECLEN},
        {(void **)&desw_e_xtal,  "DESWall_En_DetNum",    "DSW_En_Xtal",                 SUBSYS_DESWALL, 64,  E_2D_SPECLEN},
        {(void **)&desw_tof_xtal,"DESWall_TOF_DetNum",   "DSW_TOF_Xtal",                SUBSYS_DESWALL, 64,  E_2D_SPECLEN},
        {(void **)&desw_psd_e,   "DESWall_PSD_En",       "DES_Wall_PSD_En",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_tof, "DES_Wall_PSD_TOF",      "DES_Wall_PSD_TOF",            SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_q,   "DESWall_PSD_q",       "DES_Wall_PSD_q",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_cc,  "DESWall_PSD_cc",       "DES_Wall_PSD_cc",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_q_cc,  "DESWall_q_cc",       "DES_Wall_q_cc",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_q_tof,  "DESWall_q_tof",       "DES_Wall_q_tof",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_cc_tof,  "DESWall_cc_tof",       "DES_Wall_cc_tof",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&desw_psd_zdse, "DES_Wall_PSD_ZDSEn",      "DES_Wall_PSD_ZDSEn",       SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **) rcmp_strips,  "RCS%02d_E_strips",         "",                         SUBSYS_RCMP,   2*N_RCMP_STRIPS, E_2D_RCMP_SPECLEN, N_RCMP_POS},
        {NULL,                   "Hits_and_Sums/Pileup",     "",                         },
        {(void **)&ge_pu_class,  "Pile_up_class",           "",                          SUBSYS_HPGE_A,         64},
        {(void **) ge_sum_class,    "", ge_pu_class_sum_titles[0],                       SUBSYS_HPGE_A,  E_SPECLEN,   0, N_PU_CLASSES},
        {(void **) ge_e_vs_k_class, "", ge_pu_class_2d_titles[0],                        SUBSYS_HPGE_A,       2048, 512, N_PU_CLASSES},
        {(void **)&ge_pu_type,   "Pile_up_type",            "",                          SUBSYS_HPGE_A,         64},
        {(void **)&ge_nhits_type,"nhits_type",              "",                          SUBSYS_HPGE_A,         64},
        {(void **) ge_1hit,      "Ge%02d_Single_hit",       "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {(void **) ge_2hit,      "Ge%02d_2_hit_pileup",     "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {(void **) ge_3hit,      "Ge%02d_3_hit_pileup",     "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
        {NULL,                   "Hits_and_Sums/Pileup_corrections",     "",                         },
        {(void **) ge_e_vs_k_2hit_first,      "Ge%02d_E_vs_k_1st_of_2hit",              "", SUBSYS_HPGE_A,  2048, 720, 64},
        {(void **) ge_e_vs_k_2hit_second,     "Ge%02d_E_vs_k_2nd_of_2hit",              "", SUBSYS_HPGE_A,  2048, 720, 64},
        {(void **) ge_PU2_e2_v_k_gatedxrays,  "Ge%02d_PU2_E2_vs_k2_E1gated_on_Xrays",   "", SUBSYS_HPGE_A,   256, 720, 64},
        {(void **) ge_PU2_e2_v_k_gated1408,   "Ge%02d_PU2_E2_vs_k2_E1gated_on_1408keV", "", SUBSYS_HPGE_A,   256, 720, 64},
        {NULL,                   "Hits_and_Sums/BGO_HV_Alignment",     "",                         },
        {(void **) ge_bgo_gated, "",   ge_bgo_handles[0],  SUBSYS_HPGE_A,  E_2D_SPECLEN, 0, N_BGO},
        // Coinc
        {NULL,                  "Hits_and_Sums/Delta_t"," "},
        {(void **) dt_hist,     "",        dt_handles[0],  SUBSYS_HPGE_A,  DT_SPEC_LENGTH, 0, N_DT }, // leave subsys as GE -> all always defined
        {(void **) dt_tacs_hist,"dt_labr_tac%d",      "",  SUBSYS_TAC_LABR,  DT_SPEC_LENGTH, 0, N_TACS },
        {(void **) tac_lbl_ts_diff,"TAC%02d timestamp offset", "",  SUBSYS_LABR_L,  DT_SPEC_LENGTH, 0, N_TACS },
        {NULL,                  "Coinc/Coinc",        ""},
        {(void **)&gg_ab,       "Addback_GG",         "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&gg,          "GG",                 "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&ge_bgo,      "GeBGO",              "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_paces,    "GePaces",            "",  SUBSYS_PACES,   E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_labr,     "GeLabr",             "",  SUBSYS_LABR_L,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_zds,      "GeZds",              "",  SUBSYS_ZDS_A,   E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_art,      "GeAries",            "",  SUBSYS_ARIES_A, E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&labr_art,    "LaBrAries",          "",  SUBSYS_LABR_L,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&paces_art,   "PacesAries",         "",  SUBSYS_PACES,   E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&art_art,     "AriesAries",         "",  SUBSYS_ARIES_A, E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&labr_labr,   "LaBrLabr",           "",  SUBSYS_LABR_L,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&labr_zds,    "LaBrZds",            "",  SUBSYS_LABR_L,  E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&dsw_dsw,     "DSWDSW",             "",  SUBSYS_DESWALL,E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&ge_dsw,      "GeDSW",              "",  SUBSYS_DESWALL,E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&art_dsw,     "ARTDSW",             "",  SUBSYS_DESWALL,E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&ge_rcmp,     "GeRCMP",             "",  SUBSYS_RCMP,    E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&labr_rcmp,   "LaBrRCMP",           "",  SUBSYS_RCMP,    E_2D_SPECLEN, E_2D_SPECLEN},
        {(void **)&gg_opp,      "GGoppo",             "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {(void **)&gg_ab_opp,   "GG_Addback_oppo",    "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
        {NULL,                  "Coinc/Hits",         ""},
        {(void **)&gg_hit,      "GeGeHit",            "",  SUBSYS_HPGE_A,  64,  64},
        {(void **)&bgobgo_hit,  "BgoBgoHit",          "",  SUBSYS_BGO,    512, 512},
        {(void **)&gea_hit,     "GeAriesHit",         "",  SUBSYS_HPGE_A,  80,  64},
        {(void **)&lba_hit,     "LaBrAriesHit",       "",  SUBSYS_LABR_L,  16,  80},
        {(void **)&aa_hit,      "AriesAriesHit",      "",  SUBSYS_ARIES_A, 80,  80},
        {(void **)&dsw_hit,     "DSWDSWHit",          "",  SUBSYS_DESWALL,64,  64},
        {(void **) rcmp_hit,    "RCS%d_PN_hit",       "",  SUBSYS_RCMP, N_RCMP_STRIPS,     N_RCMP_STRIPS,     N_RCMP_POS},
        {(void **) rcmp_fb,     "RCS%d_Front_Back",  "",   SUBSYS_RCMP, E_2D_RCMP_SPECLEN, E_2D_RCMP_SPECLEN, N_RCMP_POS},
        {(void **)&rcmp_x_ge_hit,"RCS_Xstrips_vs_GeHit","",SUBSYS_RCMP, 192,  64},
        {(void **)&rcmp_y_ge_hit,"RCS_Ystrips_vs_GeHit","",SUBSYS_RCMP, 192,  64},
        {NULL,                  "Ang_Corr/GG_Ang_Corr",         ""},
        {(void **) gg_angcor_110,"Ge_Ge_110mm_angular_bin%02d", "", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  SYMMETERIZE, N_GE_ANG_CORR},
        {(void **) gg_angcor_145,"Ge_Ge_145mm_angular_bin%02d", "", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  SYMMETERIZE, N_GE_ANG_CORR},
        {NULL,                   "Ang_Corr/DSW_DSW_Ang_Corr",   ""},
        {(void **) dsw_angcor,   "DSW_DSW_angular_bin%03d",     "", SUBSYS_DESWALL,DSW_ANGCOR_SPECLEN, SYMMETERIZE, N_DSW_DSW_ANG_CORR },
        {NULL,                   "Ang_Corr/GG_ART_Ang_Corr",    ""},
        {(void **) ge_art_angcor,"Ge_ART_angular_bin%03d",      "", SUBSYS_ARIES_A,  GE_ANGCOR_SPECLEN, GE_ANGCOR_SPECLEN, N_GRG_ART_ANG_CORR },
        {NULL,                   "Fast_Timing/LBL_Walk",        ""},
        {(void **)&tac_labr_CompWalk0,"TAC01_LBL01_00_CompWalk",       "", SUBSYS_TAC_LABR, ECAL_TAC_SPECLEN, 1440},
        {(void **) tac_labr_CompWalk,"TAC00_LBL00_%02d_CompWalk",       "", SUBSYS_TAC_LABR, ECAL_TAC_SPECLEN, 1440, N_LABR },
        {NULL,                   "Fast_Timing/TAC_Gated_LBL_Energy", ""},
        {(void **) tac_gated_lbl,"TAC_gated_LBL%02d",           "", SUBSYS_TAC_LABR, E_SPECLEN, 0, N_LABR },
        {NULL,                   "Fast_Timing/Calibrated_TACs", ""},
        {(void **)&final_tac_sum,"Calibrated_TAC_Sum",          "", SUBSYS_TAC_LABR, E_SPECLEN },
        {(void **) final_tac,    "Calibrated_TAC%02d",          "", SUBSYS_TAC_LABR, E_SPECLEN, 0, N_LABR },
        {(void **)&lbl_lbl_tac,  "LBL_LBL_vs_TAC",              "", SUBSYS_TAC_LABR, E_3D_TAC_SPECLEN, E_3D_LBL_SPECLEN},
        {NULL,                   "Fast_Timing/ART_TACs",        ""},
        {(void **)&tac_aries_lbl_sum, "TAC_ART_LBL_LBLSUM",     "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&tac_aries_art_sum, "TAC_ART_LBL_ARTSUM",     "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&aries_tac,         "TAC_ARIES_LaBr3_1275keV","", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&aries_tac_Egate,   "TAC_ARTE_LaBr3_1275keV", "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
        {(void **)&aries_tac_artEn,   "TAC_ARIES_Energy",       "", SUBSYS_ARIES_A, E_SPECLEN      },
        {(void **) tac_aries_lbl,    "TAC_ART_LBL%d",           "", SUBSYS_ARIES_A, E_TAC_SPECLEN, E_TAC_SPECLEN, N_LABR  },
        {(void **) tac_aries_art,    "TAC_LBL_ART%d",           "", SUBSYS_ARIES_A, E_TAC_SPECLEN, E_TAC_SPECLEN, N_ARIES },
        {NULL,                   "Analysis/Isomer_Spec",        ""},
        {(void **)&gb_dt,      "betaG_time_diff_vs_gamma_energy","", SUBSYS_HPGE_A, DT_SPEC_LENGTH, E_2D_SPECLEN},
        {(void **)&gg_dt,      "GG_time_diff_vs_gamma_energy",   "", SUBSYS_HPGE_A, DT_SPEC_LENGTH, E_2D_SPECLEN},
        {(void **)&ge_isomer_popu, "GammaE_populating",          "", SUBSYS_HPGE_A,     E_SPECLEN},
        {(void **)&ge_isomer_depop,"GammaE_depopulating",        "", SUBSYS_HPGE_A,     E_SPECLEN},
        {NULL,                   "Analysis/Cycles",        ""},
        {(void **)&cycle_num_vs_ge,       "cycle_vs_Ge",         "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_sh_g,  "cycle_vs_GeE_NPU",    "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_dt,    "cycle_vs_Ge_DT",      "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_sh,       "cycle_vs_Ge_NPU",     "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_pu,       "cycle_vs_Ge_PU",      "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_b,     "cycle_vs_GeB",        "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_b_sh_g,"cycle_vs_GeBE_NPU",   "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_ge_b_dt,  "cycle_vs_GeB_DT",     "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_sh_b,     "cycle_vs_GeB_NPU",    "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&cycle_num_vs_pu_b,     "cycle_vs_GeB_PU",     "", SUBSYS_HPGE_A, MAX_CYCLES, CYCLE_SPEC_LENGTH},
        {(void **)&ge_e_vs_cycle_time,    "Time_within_cycle_vs_Ge","",SUBSYS_HPGE_A, CYCLE_SPEC_LENGTH, E_2D_SPECLEN},
        {(void **)&ge_cycle_activity, "HPGe_cycle_activity",     "", SUBSYS_HPGE_A,     CYCLE_SPEC_LENGTH},
        {(void **)&zds_cycle_activity, "ZDS_cycle_activity",     "", SUBSYS_ZDS_A,      CYCLE_SPEC_LENGTH},
        {(void **) ge_cycle_code,    "",    ge_cycle_code_titles[0], SUBSYS_HPGE_A,     E_SPECLEN, 0, N_PPG_PATTERNS},
        {(void **) gg_cycle_code,    "",    gg_cycle_code_titles[0], SUBSYS_HPGE_A,     E_SPECLEN, E_SPECLEN, N_PPG_PATTERNS},
        {(void **) gea_cycle_num,     "HPGeA_cycle%03d",         "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_sh,  "HPGeA_NPU_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_pu,  "HPGeA_PU_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_dt,  "HPGeA_DT_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_g,   "HPGeA_GeE_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) gea_cycle_num_sh_g,"HPGeA_GeE_NPU_cycle%03d", "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num,     "HPGeB_cycle%03d",         "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_sh,  "HPGeB_NPU_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_pu,  "HPGeB_PU_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_dt,  "HPGeB_DT_cycle%03d",      "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_g,   "HPGeB_GeE_cycle%03d",     "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {(void **) geb_cycle_num_sh_g,"HPGeB_GeE_NPU_cycle%03d", "", SUBSYS_HPGE_A,CYCLE_SPEC_LENGTH, 0, MAX_CYCLES },
        {NULL,                   "Analysis/Crosstalk",        ""},
        {(void **) ct_e_vs_dt_B,       "Crosstalk_Blue_E_vs_dt_Ge%02d",  "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {(void **) ct_e_vs_dt_G,       "Crosstalk_Green_E_vs_dt_Ge%02d", "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {(void **) ct_e_vs_dt_R,       "Crosstalk_Red_E_vs_dt_Ge%02d",   "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
        {(void **) ct_e_vs_dt_W,       "Crosstalk_White_E_vs_dt_Ge%02d", "", SUBSYS_HPGE_A, 864, 128, N_HPGE },
      }; // Note initialized array variable is CONST (not same as double-pointer)
      // TH1I *hist;  hist = (TH1I *) 0;   ptr = &hist = (TH1I **)addr;  *ptr =

      void get_histo_handle(char *result, char *title, char *handle )
      {
        if( strlen(handle) != 0 ){ memcpy(result, handle, strlen(handle)); return; }
        while( *title != 0 ){
          if( strcmp(title, "Energy")     == 0){ memcpy(result, "E",    1); title +=  6; result += 1; continue; }
          if( strcmp(title, "CrystalNum") == 0){ memcpy(result, "Xtal", 4); title += 10; result += 4; continue; }
          if( strcmp(title, "Upstream")   == 0){ memcpy(result, "US",   2); title +=  8; result += 2; continue; }
          if( strcmp(title, "Downstream") == 0){ memcpy(result, "DS",   2); title += 10; result += 2; continue; }
          if( strcmp(title, "Des_Wall")   == 0){ memcpy(result, "DESW", 4); title +=  8; result += 4; continue; }
          *result = *title; ++title; ++result;
        }
        *result = 0; return;
      }

      // only create histograms when that subsystem is seen in the data
      // this function is called everytime a new sybsystem is seen
      int init_histos(Config *cfg, int subsystem)
      {
        Histogram_definition *hptr;
        static Config *save_cfg;
        char tmp[STRING_LEN], **tptr;
        int i, j, k;

        if( cfg == NULL ){ cfg = save_cfg; } else { save_cfg = cfg; }
        subsys_initialized[subsystem] = 1;
        for(i=0; i<HISTO_DEF_SIZE; i++){ hptr = &histodef_array[i];
          if( hptr->ptr == NULL ){ // new folder or empty definition
            if( strlen(hptr->title) != 0 ){ memcpy(cfg->current_path, hptr->title, strlen(hptr->title)+1 ); }
            continue;
          }
          if( subsystem != hptr->subsys ){ continue; } // skip
          if( hptr->count == 0 ){ // single histogram
            get_histo_handle(tmp, hptr->title, (char *)hptr->handle );
            if( hptr->ychan == 0 ){ // 1d
              *(TH1I **)(hptr->ptr) = H1_BOOK(cfg, tmp, hptr->title, hptr->xchan, 0, hptr->xchan );
            } else {
              *(TH2I **)(hptr->ptr) = H2_BOOK(cfg, tmp, hptr->title, hptr->xchan, 0, hptr->xchan,
                hptr->ychan, 0, hptr->ychan );
              }
            } else { // array of histos (title and handle both the same)
              for(j=0; j<hptr->count; j++){
                sprintf(tmp, (strlen(hptr->title) == 0) ? hptr->handle+j*HANDLE_LENGTH : hptr->title, j);
                if( hptr->ychan == 0 ){ // 1d
                  *(TH1I **)(hptr->ptr+j) = H1_BOOK(cfg, tmp, tmp, hptr->xchan, 0, hptr->xchan );
                } else {
                  *(TH2I **)(hptr->ptr+j) = H2_BOOK(cfg, tmp, tmp, hptr->xchan, 0, hptr->xchan,
                    hptr->ychan, 0, hptr->ychan );
                  }
                }
              }
            } *cfg->current_path=0; // empty path

            // custom definitions that don't fit usual scheme
            if( subsystem == SUBSYS_TAC_LABR ){ // TAC coincidence pair spectra
              open_folder(cfg, "Fast_Timing");
              open_folder(cfg, "LBL_TAC_Combos");
              k=0; memset(tac_labr_hist_index, -1, N_LABR*N_LABR*sizeof(int));
              for(i=0; i<N_LABR; i++){
                for(j=(i+1); j<N_LABR; j++){
                  tac_labr_hist_index[i][j] = k;
                  sprintf(tmp,"TAC_%02d_%02d", i, j);
                  tac_labr_hist[k++] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
                }
              }
              sprintf(tmp,"TAC_%02d_%02d", 1, 0); // Add additional histogram (2_1) needed for Compton Walk corrections
              tac_labr_hist_index[1][0] = k;
              tac_labr_hist[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
              // and the uncalibrated versions
              k=0;
              for(i=0; i<N_LABR; i++){
                for(j=(i+1); j<N_LABR; j++){
                  sprintf(tmp,"uncalibrated_TAC_%02d_%02d", i, j);
                  tac_labr_hist_uncal[k++] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
                }
              }
              sprintf(tmp,"uncalibrated_TAC_%02d_%02d", 1, 0); // Add additional histogram (2_1) needed for Compton Walk corrections
              tac_labr_hist_uncal[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
              close_folder(cfg);
              close_folder(cfg);
            }
            // fill in subsys EvsE and Dt table pointers [** [X][<=Y]
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_HPGE_A ] = gg;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_BGO    ] = ge_bgo;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_PACES  ] = ge_paces;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_LABR_L ] = ge_labr;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_RCMP   ] = ge_rcmp;
            subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_ZDS_A  ] = ge_zds;
            subsys_e_vs_e[SUBSYS_PACES  ][SUBSYS_ARIES_A] = paces_art;
            subsys_e_vs_e[SUBSYS_LABR_L ][SUBSYS_LABR_L ] = labr_labr;
            subsys_e_vs_e[SUBSYS_LABR_L ][SUBSYS_ARIES_A] = labr_art;
            subsys_e_vs_e[SUBSYS_LABR_L ][SUBSYS_ZDS_A  ] = labr_zds;
            subsys_e_vs_e[SUBSYS_ARIES_A][SUBSYS_ARIES_A] = art_art;
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_HPGE_A  ] = dt_hist[ 0];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_PACES   ] = dt_hist[ 4];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_LABR_L  ] = dt_hist[ 5];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_RCMP    ] = dt_hist[ 6];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_ZDS_A   ] = dt_hist[ 3];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_ARIES_A ] = dt_hist[10];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_BGO     ] = dt_hist[ 1];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_SCEPTAR ] = dt_hist[ 2];
            subsys_dt[SUBSYS_HPGE_A ][SUBSYS_DESWALL ] = dt_hist[19];
            subsys_dt[SUBSYS_PACES  ][SUBSYS_LABR_L  ] = dt_hist[ 8];
            subsys_dt[SUBSYS_PACES  ][SUBSYS_ARIES_A ] = dt_hist[12];
            subsys_dt[SUBSYS_PACES  ][SUBSYS_ZDS_A   ] = dt_hist[ 7];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_LABR_L  ] = dt_hist[26];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_ARIES_A ] = dt_hist[11];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_ZDS_A   ] = dt_hist[17];
            subsys_dt[SUBSYS_LABR_L ][SUBSYS_TAC_LABR  ] = dt_hist[16];
            subsys_dt[SUBSYS_RCMP   ][SUBSYS_RCMP    ] = dt_hist[ 9];
            subsys_dt[SUBSYS_ARIES_A][SUBSYS_ARIES_A ] = dt_hist[13];
            subsys_dt[SUBSYS_ARIES_A][SUBSYS_TAC_ART ] = dt_hist[14];
            subsys_dt[SUBSYS_ZDS_A  ][SUBSYS_TAC_ZDS ] = dt_hist[15];
            subsys_dt[SUBSYS_ZDS_A  ][SUBSYS_ZDS_B   ] = dt_hist[22];
            subsys_dt[SUBSYS_DESWALL][SUBSYS_DESWALL ] = dt_hist[18];
            subsys_dt[SUBSYS_DESWALL][SUBSYS_ARIES_A ] = dt_hist[20];
            subsys_dt[SUBSYS_DESWALL][SUBSYS_ZDS_A   ] = dt_hist[21];
            return(0);
          }

          int fill_singles_histos(Grif_event *ptr)
          {
            int i, j, dt, pu, nhits, chan, pos, sys, elem, clover, c1,c2, index, offset, bin;
            char *name, c;

            chan = ptr->chan;
            // Check for invalid channel numbers, prossibly due to data corruption
            if( chan < 0 || chan > odb_daqsize ){
              fprintf(stderr,"Invalid channel number in fill_singles_histos(), %d\n",chan);
              return(-1);
            }
            sys = ptr->subsys;
            // Check this is a valid susbsytem type
            if( sys <0 || sys > MAX_SUBSYS ){
              if( mult_hist[sys] != NULL && sys>=0 && sys<MAX_SUBSYS){ mult_hist[sys]->Fill(mult_hist[sys], ptr->multiplicity, 1);  }
              return(-1);
            }
            // Get the position for this fragment
            pos  = crystal_table[ptr->chan];

            switch (sys){
              case SUBSYS_HPGE_A: // GRGa
              if( pos >= 0 && pos < 64 ){
                ge_sum->Fill(ge_sum, (int)ptr->ecal, 1);
                ge_xtal->Fill(ge_xtal, pos, (int)ptr->ecal, 1);

                // Separate sum spectra for upstream and downstream
                if(pos<32){
                  ge_sum_ds->Fill(ge_sum_ds, (int)ptr->ecal, 1);
                }else{
                  ge_sum_us->Fill(ge_sum_us, (int)ptr->ecal, 1);
                }

                // Pile-up
                pu = ptr->pileup;
                ge_pu_type->Fill(ge_pu_type, (int)pu, 1);
                nhits = ptr->nhit;
                ge_nhits_type->Fill(ge_nhits_type, (int)nhits, 1);

                // The PU class is assigned in the presort, use ptr->pu_class for pileup class
                ge_pu_class->Fill(ge_pu_class, ptr->pu_class, 1);  // pileup class value
                ge_sum_class[ptr->pu_class]->Fill(ge_sum_class[ptr->pu_class], (int)ptr->ecal, 1);  // energy spectrum of pileup class value
                ge_e_vs_k_class[ptr->pu_class]->Fill(ge_e_vs_k_class[ptr->pu_class], (int)ptr->ecal, (int)ptr->integ1, 1);  // energy spectrum of pileup class value

                if(ptr->pu_class == PU_SINGLE_HIT){  // single hit
                  ge_1hit[pos]->Fill(ge_1hit[pos], (int)ptr->ecal, 1);
                  ge_xtal_1hit->Fill(ge_xtal_1hit, pos, (int)ptr->ecal, 1);
                }
                if(ptr->pu_class >= PU_3HIT_1ST && ptr->pu_class <= PU_3HIT_3RD){ // 3-hit pileup
                  ge_3hit[pos]->Fill(ge_3hit[pos], (int)ptr->ecal, 1);
                  ge_xtal_3hit->Fill(ge_xtal_3hit, pos, (int)ptr->ecal, 1);
                }

                if(ptr->pu_class == PU_2HIT_A1ST || ptr->pu_class == PU_2HIT_B1ST || ptr->pu_class == PU_2HIT_C1ST){ // Select first Hit of two pileup events
                  ge_2hit[pos]->Fill(ge_2hit[pos], (int)ptr->ecal, 1); // 2-hit pileup
                  ge_xtal_2hit->Fill(ge_xtal_2hit, pos, (int)ptr->ecal, 1);
                  // The following used for mapping the k2 dependant correction
                  ge_e_vs_k_2hit_first[pos]->Fill(ge_e_vs_k_2hit_first[pos], (int)ptr->ecal, (int)ptr->integ1, 1);  // energy1 vs k1 spectrum of 1st Hit of 2hit pileup events
                }
                if(ptr->pu_class == PU_2HIT_A2ND || ptr->pu_class == PU_2HIT_B2ND || ptr->pu_class == PU_2HIT_C2ND){ // Select second Hit of two pileup events
                  ge_2hit[pos]->Fill(ge_2hit[pos], (int)ptr->ecal, 1); // 2-hit pileup
                  ge_xtal_2hit->Fill(ge_xtal_2hit, pos, (int)ptr->ecal, 1);
                  // The following is not used for mapping the corrections but is a useful diagnostic
                  ge_e_vs_k_2hit_second[pos]->Fill(ge_e_vs_k_2hit_second[pos], (int)ptr->ecal, (int)ptr->integ1, 1);  // energy2 vs k2 spectrum of 2nd Hit of 2hit pileup events

                  // The following 1408keV matrix used for mapping the E1 offset correction
                  if(ptr->alt_ecal > 1380 && ptr->alt_ecal < 1420){ // Require 152Eu 1408keV as E1 for mapping the E1 offset
                    ge_PU2_e2_v_k_gated1408[pos]->Fill(ge_PU2_e2_v_k_gated1408[pos], (int)ptr->ecal, (int)ptr->integ1, 1);  // e2 vs k2 for fixed e1
                  }

                  // The following x-rays matrix used for mapping the k2 dependant correction for E2
                  // E2 vs k2 gated on fixed x-ray energies
                  // The E2 energy has 1272keV subtracted from it to put the 1408keV peak around 136keV to allow a smaller matrix side and easier processing in the app
                  if(ptr->alt_ecal > 5 && ptr->alt_ecal < 130){ // Require 152Eu x rays (or 121keV because x rays are attenuated in some channels) as E1 for mapping the k2 dependance
                    ge_PU2_e2_v_k_gatedxrays[pos]->Fill(ge_PU2_e2_v_k_gatedxrays[pos], (int)(ptr->ecal - 1272), (int)ptr->integ1, 1);  // e2 vs k2 for fixed e1
                  }
                }

                clover = (int)(pos/16)+1;
                if( clover >= 0 && clover < N_CLOVER && ptr->esum >= 0 ){   // ge addback
                  ge_ab_e[clover]->Fill(ge_ab_e[clover],  (int)ptr->esum, 1);
                  ge_sum_ab   ->Fill(ge_sum_ab,     (int)ptr->esum, 1);

                  if( clover < 9 ){ // Separate Addback sum for upstream and downstream
                    ge_sum_ab_us->Fill(ge_sum_ab_us, (int)ptr->ecal, 1);
                  } else {
                    ge_sum_ab_ds->Fill(ge_sum_ab_ds, (int)ptr->ecal, 1);
                  }
                }

                // PPG Cycles histograms
                if(ppg_cycles_active==1){
                  ge_cycle_code[ppg_current_pattern]->Fill(ge_cycle_code[ppg_current_pattern], (int)ptr->ecal, 1);
                  bin = (int)((ptr->ts-ppg_cycle_start)/ppg_cycles_binning_factor);  // convert 10ns to binning size set as Global
                  ge_cycle_activity->Fill(ge_cycle_activity, bin, 1);
                  ge_e_vs_cycle_time->Fill(ge_e_vs_cycle_time, bin, (int)ptr->ecal, 1);
                  if(ppg_cycle_number<MAX_CYCLES){
                    gea_cycle_num[ppg_cycle_number]->Fill(gea_cycle_num[ppg_cycle_number], bin, 1);
                    cycle_num_vs_ge->Fill(cycle_num_vs_ge, ppg_cycle_number, bin, 1);
                    if(ptr->pu_class == PU_SINGLE_HIT){
                      gea_cycle_num_sh[ppg_cycle_number]->Fill(gea_cycle_num_sh[ppg_cycle_number], bin, 1);
                      cycle_num_vs_sh->Fill(cycle_num_vs_sh, ppg_cycle_number, bin, 1);
                    }else{
                      gea_cycle_num_pu[ppg_cycle_number]->Fill(gea_cycle_num_pu[ppg_cycle_number], bin, 1);
                      cycle_num_vs_pu->Fill(cycle_num_vs_pu, ppg_cycle_number, bin, 1);
                    }
                    if((ptr->ecal>=ppg_cycles_gamma_gate_min) && (ptr->ecal<=ppg_cycles_gamma_gate_max)){
                      gea_cycle_num_g[ppg_cycle_number]->Fill(gea_cycle_num_g[ppg_cycle_number], bin, 1);
                      if(ptr->pu_class == PU_SINGLE_HIT){
                        gea_cycle_num_sh_g[ppg_cycle_number]->Fill(gea_cycle_num_sh_g[ppg_cycle_number], bin, 1);
                        cycle_num_vs_ge_sh_g->Fill(cycle_num_vs_ge_sh_g, ppg_cycle_number, bin, 1);
                      }
                    }
                  }
                }

              }else {
                fprintf(stderr,"bad ge crystal[%d] for chan %d\n", pos, ptr->chan);
              } break;
              case SUBSYS_HPGE_B: // GRGb
              if( pos >= 0 && pos < 64 ){
                // PPG Cycles histograms
                if(ppg_cycles_active==1){
                  bin = (int)((ptr->ts-ppg_cycle_start)/ppg_cycles_binning_factor);  // convert 10ns to binning size set as Global
                  if(ppg_cycle_number<MAX_CYCLES){
                    geb_cycle_num[ppg_cycle_number]->Fill(geb_cycle_num[ppg_cycle_number], bin, 1);
                    cycle_num_vs_ge_b->Fill(cycle_num_vs_ge_b, ppg_cycle_number, bin, 1);
                    if(ptr->pu_class == PU_SINGLE_HIT){
                      geb_cycle_num_sh[ppg_cycle_number]->Fill(geb_cycle_num_sh[ppg_cycle_number], bin, 1);
                      cycle_num_vs_sh_b->Fill(cycle_num_vs_sh_b, ppg_cycle_number, bin, 1);
                    }else{
                      geb_cycle_num_pu[ppg_cycle_number]->Fill(geb_cycle_num_pu[ppg_cycle_number], bin, 1);
                      cycle_num_vs_pu_b->Fill(cycle_num_vs_pu_b, ppg_cycle_number, bin, 1);
                    }
                    if((ptr->ecal>=ppg_cycles_gamma_gate_min) && (ptr->ecal<=ppg_cycles_gamma_gate_max)){
                      geb_cycle_num_g[ppg_cycle_number]->Fill(geb_cycle_num_g[ppg_cycle_number], bin, 1);
                      if(ptr->pu_class == PU_SINGLE_HIT){
                        geb_cycle_num_sh_g[ppg_cycle_number]->Fill(geb_cycle_num_sh_g[ppg_cycle_number], bin, 1);
                        cycle_num_vs_ge_b_sh_g->Fill(cycle_num_vs_ge_b_sh_g, ppg_cycle_number, bin, 1);
                      }
                    }
                  }
                }
              } break;
              case SUBSYS_BGO: // BGOs
              pos  = crystal_table[ptr->chan];
              elem = element_table[ptr->chan];
              if( pos < 0 || pos > 63 ){
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
              }  break;
              case SUBSYS_LABR_BGO: // Ancillary BGOs
              pos  = crystal_table[ptr->chan];
              elem = element_table[ptr->chan];
              if( pos < 1 || pos > 8 ){
                fprintf(stderr,"bad ancillary bgo crystal[%d] for chan %d\n", pos, ptr->chan);
              } else if( elem < 1 || elem > 3 ){
                fprintf(stderr,"bad ancillary bgo element[%d] for chan %d\n", elem, ptr->chan);
              } else {
                pos *= 3; pos += (elem-1);
                bgoa_xtal->Fill(bgoa_xtal, pos, (int)ptr->ecal, 1);
              } break;
              case SUBSYS_PACES: // PACES
              paces_sum->Fill(paces_sum, (int)ptr->ecal, 1);
              pos  = crystal_table[ptr->chan];
              if( pos < 1 || pos > 5 ){
                fprintf(stderr,"bad PACES crystal[%d] for chan %d\n", pos, ptr->chan);
              } else {
                paces_xtal->Fill(paces_xtal, pos, (int)ptr->ecal, 1);
              } break;
              case SUBSYS_LABR_L: // LaBr3 (LBL)
              labr_sum->Fill(labr_sum, (int)ptr->ecal, 1);
              pos  = crystal_table[ptr->chan];
              if( pos < 1 || pos > 8 ){
                fprintf(stderr,"bad LaBr3 crystal[%d] for chan %d\n", pos, ptr->chan);
              } else {
                labr_xtal->Fill(labr_xtal, pos, (int)ptr->ecal, 1);
              } break;
              case SUBSYS_SCEPTAR: break;
              case SUBSYS_TAC_LABR:

              // Save LBL channel number into ptr->integ2 or integ3 or integ4
              // Save LBL energy ecal into TAC ptr-ecal2 or ecal3 or ecal4
              if( ptr->q4 > 0 ){ break; } // more than two LaBr3 in coincidence with this TAC event so reject
              if( ptr->integ2 >=          0 && ptr->integ3 >=           0 &&
                ptr->integ2 < MAX_DAQSIZE && ptr->integ3  < MAX_DAQSIZE ){
                  c1 = crystal_table[ptr->integ2]-1; // c1 runs from 0 to 7. c1 is position of first LBL in coincidence with this TAC.
                  c2 = crystal_table[ptr->integ3]-1; // c2 runs from 0 to 7. c2 is position of second LBL in coincidence with this TAC.
                } else { c1 = c2 = -1; }
                if(c1>=0 && c1<N_LABR){
                  tac_gated_lbl[c1]->Fill(tac_gated_lbl[c1], (int)(ptr->q2), 1); // First LBL energy spectrum, requiring a TAC coincidence
                  if(c2>=0 && c2<N_LABR){
                    index = tac_labr_hist_index[c1][c2];
                    if(index>=0 && index<(int)((N_LABR)*(N_LABR-1)/2)+2){
                      offset = tac_lbl_combo_offset[index];
                      tac_labr_hist[index]->Fill(tac_labr_hist[index], (int)(ptr->ecal)+offset, 1);
                      tac_labr_hist_uncal[index]->Fill(tac_labr_hist_uncal[index], (int)(( ptr->integ1 == 0 ) ? ptr->q1 : spread(ptr->q1)/ptr->integ1), 1);
                    }
                    // Calibrated TAC spectra
                    final_tac[crystal_table[ptr->chan]-1]->Fill(final_tac[crystal_table[ptr->chan]-1], (int)(ptr->ecal)+offset, 1);
                    final_tac_sum->Fill(final_tac_sum, (int)(ptr->ecal)+offset, 1);

                    // A 3d histogram of first LBL energy vs second LBL energy vs TAC
                    // lbl_lbl_tac->Fill(lbl_lbl_tac, (int)((ptr->e2cal/10)*(ptr->e3cal/10)), (int)(ptr->ecal)+offset, 1); // LBL energy vs LBL energy vs TAC
                    if(ptr->ecal>5 && ptr->q2>5 && ptr->q3>5){
                      bin = (int)(((ptr->q2/10)*400)+(ptr->q3/10));
                      lbl_lbl_tac->Fill(lbl_lbl_tac, (int)(ptr->ecal)+offset-250, bin, 1); // LBL energy vs LBL energy vs TAC
                    }

                    // Compton Walk matrix for calibrations
                    // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
                    if(ptr->ecal>5 && crystal_table[ptr->chan] == 1){ // Use the First TAC (TAC01)
                      if(c1 == 0 && c2>0 && ptr->q2>1252 && ptr->q2<1412 && ptr->q3>5){ // LBL01 gated on 1332keV
                        tac_labr_CompWalk[c2]->Fill(tac_labr_CompWalk[c2], (int)(ptr->ecal)+offset, (int)ptr->q3, 1); // TAC01 and other LBL energy
                      }
                    }else if(ptr->ecal>5 && crystal_table[ptr->chan] == 2){
                      if(c1 == 1 && c2 == 2 && ptr->q2>1252 && ptr->q2<1412 && ptr->q3>5){ // LBL02 gated on 1332keV
                        tac_labr_CompWalk0->Fill(tac_labr_CompWalk0, (int)(ptr->ecal)+offset, (int)ptr->q3, 1); // TAC02 to check LBL01
                      }
                    }
                  }
                } break;
                case SUBSYS_DESCANT: break;
                case SUBSYS_DESWALL: // DESCANT Wall
                pos  = crystal_table[chan];
                if(pos>0 && pos<=N_DES_WALL){
                  desw_psd[pos]       -> Fill(desw_psd[pos],   (int)ptr->psd,       1);
                  if(ptr->tof>0){
                    desw_tof[pos]       -> Fill(desw_tof[pos],   (int)ptr->tof,       1);
                    desw_tof_corr[pos]  -> Fill(desw_tof_corr[pos],   (int)ptr->alt_ecal,       1);
                    if(ptr->psd>10 && ptr->psd<710){
                      desw_tof_psd[pos]  -> Fill(desw_tof_psd[pos],   (int)ptr->alt_ecal,       1);
                    }
                  }
                }
                desw_sum_e->Fill(desw_sum_e, (int)ptr->ecal, 1);
                if(ptr->psd>0){ desw_sum_psd->Fill(desw_sum_psd, (int)ptr->psd, 1); }
                if(ptr->alt_ecal>0){ desw_sum_tof->Fill(desw_sum_tof, (int)ptr->alt_ecal, 1); } // alt_ecal = corrected time-of-flight
                pos  = crystal_table[ptr->chan];
                if( pos < 1 || pos > 60 ){
                  fprintf(stderr,"bad descant wall detector[%d] for chan %d\n", pos, ptr->chan);
                } else {
                  if(ptr->ecal>5){
                    desw_e_xtal->Fill(desw_e_xtal, pos, (int)ptr->ecal, 1);
                    if(ptr->psd>5){
                      desw_psd_e->Fill(desw_psd_e, (int)ptr->psd, (int)ptr->ecal, 1);
                      desw_psd_q->Fill(desw_psd_q, (int)ptr->psd, (int)ptr->q1, 1);
                      desw_psd_cc->Fill(desw_psd_cc, (int)ptr->psd, (int)ptr->cc_short, 1);
                      desw_q_cc->Fill(desw_q_cc, (int)ptr->q1, (int)ptr->cc_short, 1);
                    }
                  }
                  if(ptr->alt_ecal>5){ // DESCANT Wall ptr->alt_ecal=corrected-TOF (ptr->tof is TOF)
                    desw_tof_xtal->Fill(desw_tof_xtal, pos, (int)ptr->alt_ecal, 1);
                    if(ptr->psd>5){ desw_psd_tof->Fill(desw_psd_tof, (int)ptr->psd, (int)ptr->alt_ecal, 1); }
                  }
                }
                break;
                case SUBSYS_ARIES_A: // ARIES Standard Output
                aries_sum->Fill(aries_sum, (int)ptr->ecal, 1);
                pos  = crystal_table[ptr->chan];
                if( pos < 1 || pos > 76 ){
                  fprintf(stderr,"bad aries tile[%d] for chan %d\n", pos, ptr->chan);
                } else {
                  aries_xtal->Fill(aries_xtal, pos, (int)ptr->ecal, 1);
                } break;
                case SUBSYS_ZDS_A:
                gc_hist->Fill(gc_hist, 2, 1);

                // PPG Cycles histograms
                if(ppg_cycles_active==1){
                  bin = (int)((ptr->ts-ppg_cycle_start)/ppg_cycles_binning_factor); // binning size set in Global
                  zds_cycle_activity->Fill(zds_cycle_activity, bin, 1);
                }
                break;
                case SUBSYS_ZDS_B: gc_hist->Fill(gc_hist, 1, 1); break;
                case SUBSYS_RCMP:
                rcmp_sum->Fill(rcmp_sum, (int)ptr->ecal, 1);
                if(ptr->esum>0){
                  rcmp_fb_sum->Fill(rcmp_fb_sum, (int)ptr->esum, 1);
                }
                pos  = crystal_table[ptr->chan];
                elem = (int)(element_table[ptr->chan] + (int)(polarity_table[ptr->chan]*N_RCMP_STRIPS)); // polarity_table value is 0 or 1
                if( pos < 1 || pos > 6 ){
                  fprintf(stderr,"bad RCMP DSSD[%d] for chan %d, elem%d, pol%d\n", pos, ptr->chan, elem, polarity_table[ptr->chan]);
                } else if( elem < 0 || elem > 63 ){
                  fprintf(stderr,"bad RCMP strip[%d] for chan %d, pos%d, pol%d\n", elem, ptr->chan, pos, polarity_table[ptr->chan]);
                } else {
                  rcmp_strips[(pos-1)]->Fill(rcmp_strips[(pos-1)], elem, (int)ptr->ecal, 1);
                }
                break;
                default: break; // Unrecognized or unprocessed dtype
              }// end of switch
              return(0);
            }

            int fill_ge_coinc_histos(Grif_event *ptr, Grif_event *alt, int abs_dt)
            {
              int c1, c2, pos, bin, angle_idx;
              switch(alt->subsys){
                case SUBSYS_HPGE_A:
                gg_dt->Fill(gg_dt, (int)((ptr->ts - alt->ts)+DT_SPEC_LENGTH/2), (int)ptr->ecal, 1); // This dt result is always negative
                gg_dt->Fill(gg_dt, (int)((alt->ts - ptr->ts)+DT_SPEC_LENGTH/2), (int)alt->ecal, 1); // This dt result is always positive
                if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) ){
                  if( ptr->esum >= 0 &&  alt->esum >= 0 ){ // addback energies
                    gg_ab->Fill(gg_ab, (int)ptr->esum, (int)alt->esum, 1);
                  }
                  // PPG Cycles histograms
                  if(ppg_cycles_active==1){
                    gg_cycle_code[ppg_current_pattern]->Fill(gg_cycle_code[ppg_current_pattern], (int)ptr->ecal, (int)alt->ecal, 1);
                  }
                  c1 = crystal_table[ptr->chan];
                  c2 = crystal_table[alt->chan];
                  if( c1 >= 0 && c1 < 64 && c2 >= 0 && c2 < 64 ){
                    gg_hit->Fill(gg_hit, c1, c2, 1); // 2d crystal hitpattern
                    if( c2 == grif_opposite[c1] ){
                      // 180 degree coinc matrix for summing corrections
                      gg_opp->Fill(gg_opp, (int)ptr->ecal, (int)alt->ecal, 1);
                      gg_ab_opp->Fill(gg_ab_opp, (int)ptr->esum, (int)alt->esum, 1);
                    }
                    // Ge-Ge angular correlations
                    // Fill the appropriate angular bin spectrum
                    // c1 and c2 run from 0 to 63 for ge_angles_145mm.
                    angle_idx = ge_angles_110mm[c1][c2];
                    gg_angcor_110[angle_idx]->Fill(gg_angcor_110[angle_idx], (int)ptr->ecal, (int)alt->ecal, 1);
                  //  fprintf(stdout,"%d %d have angular difference of %lf\n",c1,c2,angular_diff_GeGe(c1,c2,110));

                    angle_idx = ge_angles_145mm[c1][c2];
                    gg_angcor_145[angle_idx]->Fill(gg_angcor_145[angle_idx], (int)ptr->ecal, (int)alt->ecal, 1);

                  }
                }
                break;
                case SUBSYS_BGO:
                if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_BGO]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_BGO]) ){
                  c1 = crystal_table[ptr->chan];
                  c2 = crystal_table[alt->chan];
                  if( c1 >= 0 && c1 < 64 && c1==c2 && ptr->ecal>5 && alt->ecal>5){
                    bin = (c2*5)+(element_table[alt->chan]-1);
                    if(bin>=0 && bin<N_BGO){
                      ge_bgo_gated[bin]->Fill(ge_bgo_gated[bin], (int)alt->ecal, 1);
                    }
                  }
                }
                break;
                case SUBSYS_SCEPTAR:
                gb_dt->Fill(gb_dt, (int)((ptr->ts - alt->ts)+DT_SPEC_LENGTH/2), (int)ptr->esum, 1);
                if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) ){
                  ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
                  ge_sum_b_sep->Fill(ge_sum_b_sep, (int)ptr->ecal, 1); // Sceptar-gated Ge sum energy spectrum
                }else if((ptr->ts - alt->ts)<-25){
                  ge_isomer_popu->Fill(ge_isomer_popu, (int)ptr->ecal, 1); // Early gamma rays appearing earlier in time than the prompt
                }else if((ptr->ts - alt->ts)>25){
                  ge_isomer_depop->Fill(ge_isomer_depop, (int)ptr->ecal, 1); // Delayed gamma rays appearing later in time than the prompt
                }
                break;
                case SUBSYS_ARIES_A:
                gb_dt->Fill(gb_dt, (int)((ptr->ts - alt->ts)+DT_SPEC_LENGTH/2), (int)ptr->esum, 1);
                if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) ){
                  if(ptr->ecal > 10 && alt->ecal > 10){
                    ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1);         // beta-gated Ge sum energy spectrum
                    ge_sum_b_art->Fill(ge_sum_b_art, (int)ptr->ecal, 1); // Aries-gated Ge sum energy spectrum
                    ge_art->Fill(ge_art, (int)ptr->ecal, (int)alt->esum, 1);
                    c1 = crystal_table[ptr->chan];
                    c2 = crystal_table[alt->chan];
                    if( c1 >= 0 && c1 < 64 && c2 >= 1 && c2 <= 76 ){
                      gea_hit->Fill(gea_hit, c2, c1, 1); // c's start at zero
                      // Ge-ARIES angular correlations
                      // Fill the appropriate angular bin spectrum
                      angle_idx = GRG_ART_angles_110mm[c1][c2-1];
                      ge_art_angcor[angle_idx]->Fill(ge_art_angcor[angle_idx], (int)ptr->ecal, (int)alt->ecal, 1);

                      // Angle veto method for Bremmstrahlung veto. Angle >30 degrees
                      if(angle_idx>8){
                        ge_sum_b_art_brems->Fill(ge_sum_b_art_brems, (int)ptr->ecal, 1); // Aries-gated Ge sum energy spectrum with Bremmstrahlung veto
                      }
                    }
                  }
                }else if((ptr->ts - alt->ts)<0){
                  ge_isomer_popu->Fill(ge_isomer_popu, (int)ptr->ecal, 1); // Early gamma rays appearing earlier in time than the prompt
                }else if((ptr->ts - alt->ts)>0){
                  ge_isomer_depop->Fill(ge_isomer_depop, (int)ptr->ecal, 1); // Delayed gamma rays appearing later in time than the prompt
                }
                break;
                case SUBSYS_ZDS_A:
                gb_dt->Fill(gb_dt, (int)((ptr->ts - alt->ts)+DT_SPEC_LENGTH/2), (int)ptr->esum, 1);
                if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ZDS_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ZDS_A]) ){
                  ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1);         // beta-gated Ge sum energy spectrum
                  ge_sum_b_zds->Fill(ge_sum_b_zds, (int)ptr->ecal, 1); // Zds-gated Ge sum energy spectrum
                }else if((ptr->ts - alt->ts)<0){
                  ge_isomer_popu->Fill(ge_isomer_popu, (int)ptr->ecal, 1); // Early gamma rays appearing earlier in time than the prompt
                }else if((ptr->ts - alt->ts)>0){
                  ge_isomer_depop->Fill(ge_isomer_depop, (int)ptr->ecal, 1); // Delayed gamma rays appearing later in time than the prompt
                }
                break;
                case SUBSYS_DESWALL: // ge-DSW
                ge_dsw->Fill(ge_dsw, (int)ptr->ecal, (int)alt->alt_ecal, 1); // alt_ecal = DSW corrected time-of-flight
                break;
                case SUBSYS_RCMP: // ge-RCMP
                c1 = crystal_table[ptr->chan]; // HPGe crystal number
                pos  = crystal_table[alt->chan]; // RCMP DSSD number
                c2 = (int)(element_table[alt->chan] + (int)((pos-1)*N_RCMP_STRIPS)); // RCMP strip number
                if( c1 >= 0 && c1 < 64 && c2 >= 0 && c2 <= 192 && ptr->ecal>5 && alt->ecal>5 ){
                  if(polarity_table[alt->chan]==0){ rcmp_x_ge_hit->Fill(rcmp_x_ge_hit, c2, c1, 1); }
                  else{ rcmp_y_ge_hit->Fill(rcmp_y_ge_hit, c2, c1, 1); }
                }
                break;
                default: break;
              }
              return(0);
            }

            int fill_labr_coinc_histos(Grif_event *ptr, Grif_event *alt, int abs_dt)
            {
              int lbl_tac_gate=15, tac_offset[8] = {-7300,-5585,-6804,0,-6488,-5682,-5416,0};
              int g_aries_upper_gate=25, c1, c2, dt, corrected_tac_value;
              switch(alt->subsys){
                case SUBSYS_ARIES_A:
                if( (abs_dt >= time_diff_gate_min[SUBSYS_LABR_L][SUBSYS_ARIES_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_LABR_L][SUBSYS_ARIES_A]) ){
                  c1 = crystal_table[ptr->chan];
                  c2 = crystal_table[alt->chan];
                  if(c1 >= 1 && c1 <=8 && c2 >= 1 && c2 <=76 ){
                    lba_hit->Fill(lba_hit, c1, c2, 1);
                  }
                } break;
                case SUBSYS_TAC_LABR:
                c1=crystal_table[alt->chan]-1;  // assign c1 as TAC number
                if(c1 >= 0 && c1 < N_TACS ){ // 8 LBL + 1 ZDS + 4 ARIES
                  dt_tacs_hist[c1]->Fill(dt_tacs_hist[c1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
                  tac_lbl_ts_diff[c1]->Fill(tac_lbl_ts_diff[c1], (int)((alt->ts-ptr->ts)+DT_SPEC_LENGTH/2), 1);
                  if(((abs_dt >= time_diff_gate_min[SUBSYS_LABR_L][SUBSYS_TAC_LABR]) && (abs_dt <= time_diff_gate_max[SUBSYS_LABR_L][SUBSYS_TAC_LABR])) && c1 == 8){ // labr-tac with the ARIES TAC
                    c2 = crystal_table[ptr->chan] - 1; // assign c2 as LBL number
                    if(c2 >= 0 && c2 < 8 ){ // 8 LBL detectors
                      corrected_tac_value = (int)alt->ecal + tac_offset[c2];
                      tac_aries_lbl[c2]->Fill(tac_aries_lbl[c2], corrected_tac_value, 1); // tac spectrum per LBL
                      tac_aries_lbl_sum->Fill(tac_aries_lbl_sum, corrected_tac_value, 1); // sum tac spectrum including all LBL
                      if(ptr->ecal >1225 && ptr->ecal <1315){ // gate on LaBr3 energy 1275keV
                        aries_tac->Fill(aries_tac, (int)corrected_tac_value, 1); // tac spectrum gated on 1275keV
                        aries_tac_artEn->Fill(aries_tac_artEn, alt->alt_ecal, 1); // ARIES energy spectrum in coincidence with TAC
                        if(alt->alt_ecal >24 && alt->alt_ecal <36){ // gate on ARIES energy
                          aries_tac_Egate->Fill(aries_tac_Egate, corrected_tac_value, 1); // tac spectrum gated on 1275keV
                        }
                      }
                    }
                  }
                } break;
                case SUBSYS_TAC_ZDS:
                c1=crystal_table[alt->chan]-1;  // assign c1 as TAC number
                if(c1 == 7){ // 8 LBL + 1 ZDS + 4 ARIES
                  tac_lbl_ts_diff[c1]->Fill(tac_lbl_ts_diff[c1], (int)((alt->ts-ptr->ts)+DT_SPEC_LENGTH/2), 1);
                } break;
                case SUBSYS_TAC_ART:
                c1=crystal_table[alt->chan]-1;  // assign c1 as TAC number
                if(c1 >= 9 && c1 < N_TACS ){ // 8 LBL + 1 ZDS + 4 ARIES
                  tac_lbl_ts_diff[c1]->Fill(tac_lbl_ts_diff[c1], (int)((alt->ts-ptr->ts)+DT_SPEC_LENGTH/2), 1);
                } break;
              }
              return(0);
            }

            // This lookup table reorders strips that have already been reordered in the ODB...
            // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25966
            // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25968
            int reorder_rcmp_strips[7][2][32] = {
              {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
              {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS1 X
              {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // RCS1 Y
              {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS2 X
              {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // RCS2 X
              {{1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30},  // RCS3 X
              {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // RCS3 Y
              {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS4 X
              {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}, // RCS4 Y
              {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS5 X
              {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30}}, // RCS5 Y
              {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},  // RCS6 X
              {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}}  // RCS6 Y
            };

            int frag_hist[PTR_BUFSIZE];
            int fill_coinc_histos(int win_idx, int frag_idx)
            {
              int global_window_size = (int)(sort_window_width/2); // size in grif-replay should be double this
              Grif_event *alt, *ptr = &grif_event[win_idx], *tmp;
              int dt, abs_dt,  pos, c1, c2, index, ptr_swap;
              TH2I *hist_ee; TH1I *hist_dt;

              // histogram of coincwin-size
              dt = (frag_idx - win_idx + 2*PTR_BUFSIZE) %  PTR_BUFSIZE; ++frag_hist[dt];

              while( win_idx != frag_idx ){ // check all conicidences in window
                if( ++win_idx == PTR_BUFSIZE ){ win_idx = 0; } // wrap
                alt = &grif_event[win_idx];
                if( ptr->subsys > alt->subsys ){ tmp = ptr; ptr = alt; alt = tmp; ptr_swap = 1; }

                abs_dt = dt = ptr->ts - alt->ts; if( dt < 0 ){ abs_dt = -1*dt; }
                if( abs_dt > global_window_size ){ break; }

                // the usual subsys-vs-subsys 1d-time-diff and 2d-EvsE
                if( (hist_dt = subsys_dt[ptr->subsys][alt->subsys]) != NULL ){
                  hist_dt->Fill(hist_dt, (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
                }
                if( (hist_ee = subsys_e_vs_e[ptr->subsys][alt->subsys]) != NULL ){
                  if((abs_dt >= time_diff_gate_min[ptr->subsys][alt->subsys]) && (abs_dt <= time_diff_gate_max[ptr->subsys][alt->subsys]) ){
                    hist_ee->Fill(hist_ee, (int)ptr->ecal, (int)alt->ecal, 1);
                  }
                }
                switch(ptr->subsys){ // No Nested switch - use separate functions if needed
                  case SUBSYS_HPGE_A: fill_ge_coinc_histos(ptr,   alt, abs_dt); break;
                  case SUBSYS_LABR_L: fill_labr_coinc_histos(ptr, alt, abs_dt); break;
                  case SUBSYS_BGO:
                  if(alt->subsys == SUBSYS_BGO){
                    if((abs_dt >= time_diff_gate_min[SUBSYS_BGO][SUBSYS_BGO]) && (abs_dt <= time_diff_gate_max[SUBSYS_BGO][SUBSYS_BGO]) ){
                      c1 = crystal_table[ptr->chan];
                      c2 = crystal_table[alt->chan];
                      bgobgo_hit->Fill(bgobgo_hit, c1, c2, 1);
                    }}
                    break;
                    case SUBSYS_RCMP:
                    if(alt->subsys == SUBSYS_RCMP){
                      if( (abs_dt >= time_diff_gate_min[SUBSYS_RCMP][SUBSYS_RCMP]) && (abs_dt <= time_diff_gate_max[SUBSYS_RCMP][SUBSYS_RCMP]) ){
                        if((pos = crystal_table[ptr->chan]) == crystal_table[alt->chan] &&
                        polarity_table[ptr->chan] != polarity_table[alt->chan] ){ // front and back of same DSSD
                          c1 = element_table[ptr->chan];
                          c2 = element_table[alt->chan];
                          rcmp_fb[(pos-1)]->Fill(rcmp_fb[(pos-1)], (int)ptr->ecal, (int)alt->ecal, 1);
                          if(polarity_table[ptr->chan]==0){ rcmp_hit[(pos-1)]->Fill(rcmp_hit[(pos-1)], c1, c2, 1);
                          }else{
                            rcmp_hit[(pos-1)]->Fill(rcmp_hit[(pos-1)], c2, c1, 1);
                          }
                        }}} break;
                        case SUBSYS_ARIES_A:
                        if( alt->subsys == SUBSYS_ARIES_A && ptr->ecal>0 && alt->ecal>0 ){
                          if((abs_dt >= time_diff_gate_min[SUBSYS_ARIES_A][SUBSYS_ARIES_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_ARIES_A][SUBSYS_ARIES_A])){
                            c1 = crystal_table[ptr->chan];
                            c2 = crystal_table[alt->chan];
                            if( c1 >= 1 && c1 <=76 && c2 >= 1 && c2 <=76 ){
                              aa_hit->Fill(aa_hit, c1, c2, 1);
                            }
                          }}
                          if( alt->subsys == SUBSYS_TAC_LABR && crystal_table[alt->chan] == 8 ){ // ARIES TAC
                            // sum tac spectrum including all art
                            tac_aries_art_sum->Fill(tac_aries_art_sum, (int)alt->ecal, 1);
                            c2 = crystal_table[ptr->chan]-1;
                            if( c2>=0 && c2<N_ARIES ){ // tac spectrum per ART tiles
                              tac_aries_art[c2]->Fill(tac_aries_art
                                [c2], (int)alt->ecal, 1);
                              }
                            } break;
                            case SUBSYS_ARIES_B:// ARIES Fast Output in CAEN
                            if(alt->subsys == SUBSYS_DESWALL){ // aries-DSW
                              art_dsw->Fill(art_dsw, (int)ptr->ecal, (int)alt->alt_ecal, 1);
                              desw_sum_e_b->Fill(desw_sum_e_b, (int)alt->ecal, 1);
                              desw_sum_tof_b->Fill(desw_sum_tof_b, (int)alt->alt_ecal, 1); // alt_ecal = corrected time-of-flight
                            } break;
                            case SUBSYS_ZDS_A: // grif16 zds
                            if(alt->subsys == SUBSYS_ZDS_B ){ // ZDS GRIF-CAEN coincidence
                              gc_hist->Fill(gc_hist, 3, 1);
                              gc_hist->Fill(gc_hist, 5, 1);
                              dt_hist[23]->Fill(dt_hist[23], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
                            } break;
                            case SUBSYS_DESWALL:
                            if( alt->subsys == SUBSYS_DESWALL ){
                              dt_hist[24]->Fill(dt_hist[24], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
                              c1 = crystal_table[ptr->chan]-1; c2 = crystal_table[alt->chan]-1;
                              if( c1 >= 0 && c1 < 60 && c2 >= 0 && c2 < 60 ){
                                dsw_hit->Fill(dsw_hit, c1, c2, 1);
                                dsw_dsw->Fill(dsw_dsw, (int)ptr->alt_ecal, (int)alt->alt_ecal, 1);
                                // Fold 2 sum spectra
                                desw_sum_e_nn->Fill(desw_sum_e_nn, (int)ptr->ecal, 1);
                                desw_sum_tof_nn->Fill(desw_sum_tof_nn, (int)ptr->alt_ecal, 1);// alt_ecal = corrected time-of-flight
                                // DSW-DSW angular correlations
                                // Fill the appropriate angular bin spectrum with the corrected time-of-flight value
                                index = DSW_DSW_angles[c1][c2];
                                dsw_angcor[index]->Fill(dsw_angcor
                                  [index],(int)ptr->alt_ecal,(int)alt->alt_ecal, 1);
                                  // Fold 2, angle greater than 60 degrees, sum spectra
                                  // index 13 = 58.555, index 14 = 61.535
                                  if( index > 13 ){                                         // alt_ecal = corrected time-of-flight
                                    desw_sum_e_nn_a->Fill(desw_sum_e_nn_a, (int)ptr->ecal, 1);
                                    desw_sum_tof_nn_a->Fill(desw_sum_tof_nn_a, (int)ptr->alt_ecal, 1);
                                  }
                                }
                              }
                              if( alt->subsys == SUBSYS_ZDS_B ){ // ZDS-DSW
                                dt_hist[25]->Fill(dt_hist[25], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
                                desw_sum_e_b->Fill(desw_sum_e_b, (int)ptr->ecal, 1);
                                desw_sum_tof_b->Fill(desw_sum_tof_b, (int)ptr->alt_ecal, 1); // alt_ecal = corrected time-of-flight

                                desw_q_tof->Fill(desw_q_tof, (int)ptr->q1, (int)ptr->alt_ecal, 1);
                                desw_cc_tof->Fill(desw_cc_tof, (int)ptr->cc_short, (int)ptr->alt_ecal, 1);
                              } break;
                            } // end switch
                          }
                          return(0);
                        }

                        //#######################################################################
                        //###########   READ XML ODB DUMP FROM START OF DATA FILE   #############
                        //#######################################################################

                        // Note - the odb names do not distinguish between subtypes of detectors
                        // e.g Ge A and B channels
                        // the subsystem names will be extended to include this information
                        // (and the odb-specific names are only used below)

                        #define MAX_ODB_SUBSYS 24
                        #define ODBHANDLE_GRG   0
                        #define ODBHANDLE_GRS   1
                        #define ODBHANDLE_SEP   2
                        #define ODBHANDLE_PAC   3
                        #define ODBHANDLE_LBS   4
                        #define ODBHANDLE_LBT   5
                        #define ODBHANDLE_LBL   6
                        #define ODBHANDLE_DSC   7
                        #define ODBHANDLE_ART   8
                        #define ODBHANDLE_ZDS   9
                        #define ODBHANDLE_RCS  10
                        #define ODBHANDLE_XXX  11
                        #define ODBHANDLE_DSW  12
                        #define ODBHANDLE_DSG  13
                        #define ODBHANDLE_DAL  14
                        #define ODBHANDLE_DAT  15
                        #define ODBHANDLE_UNK  23
                        static char odb_handle[MAX_ODB_SUBSYS][8] = {
                          "GRG", "GRS", "SEP", "PAC",  //  0- 3
                          "LBS", "LBT", "LBL", "DSC",  //  4- 7
                          "ART", "ZDS", "RCS", "XXX",  //  8-11
                          "DSW", "DSG", "DAL", "DAT",  //  12-15
                          "",    "",    "",    "",
                          "",    "",    "",    "UNK"
                        };

                        static char   path[256];
                        static char dirname[64],value[32],type[32];
                        extern char midas_runtitle[SYS_PATH_LENGTH];

                        static void *arrayptr;
                        int read_odb_items(int len, int *bank_data)
                        {
                          char *path_ptr, *ptr, *str, *odb_data = (char *)bank_data, posn[2], odb_ppg_current[128];
                          int i, j, c = '<', d = '>', dtype=0, active=0, index=0, ppg_index, odb_ppg_cycle_count=-1;

                          // The currently set cycle is included in the ODB dump after all cycles are defined.
                          // So we need to unpack them all and then select the values we need.
                          ppg_cycles odb_ppg_cycle[MAX_ODB_PPG_CYCLES];

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
                              if( strncmp(path,"/PPG/Cycles/",12) == 0 ){
                                odb_ppg_cycle_count++;
                                strcpy(odb_ppg_cycle[odb_ppg_cycle_count].name, dirname);
                              }
                            } else if( strncmp(ptr,"</dir>",6) == 0 ){
                              while(1){
                                if( --path_ptr < path ){ path_ptr = path;  *path_ptr = '\0';  break; }
                                if( *path_ptr == '/' ){ *path_ptr = '\0';  break; }
                              }
                              index=0; // for debugger to stop here
                            } else if( strncasecmp(ptr,"<key name=\"Run Title\" type=\"STRING\"", 35) == 0 ){
                              ptr = str+1;
                              if( (str = strchr(ptr,c)) == NULL ){ break; }
                              i = (str-ptr) > SYS_PATH_LENGTH-1 ? SYS_PATH_LENGTH-1 : (str-ptr);
                              memcpy( midas_runtitle, ptr, i ); midas_runtitle[i] = 0;
                              ptr += i+1;
                              if( (str = strchr(ptr,d)) == NULL ){ break; }
                            } else if( strncasecmp(ptr,"<key name=\"Current\"", 19) == 0 &&
                            strncmp(path,"/PPG",4) == 0 ){
                              *str = ' '; if( (str = strchr(str,c)) == NULL ){ break; }
                              *str = ' '; if( (str = strchr(str,d)) == NULL ){ break; }
                              if( sscanf(ptr,"<key name=\"Current\" type=\"STRING\" size=\"%d\" %s /key>", &dtype, value) < 2 ){
                                fprintf(stderr,"can't read key name for Current PPG cycle\n"); ptr=str+1; continue;
                              }
                              strcpy(odb_ppg_current,value);
                              ptr += i+1;
                            } else if( strncasecmp(ptr,"<key name=\"prg_ddtm\"", 20) == 0 &&
                            strncmp(path,"/DAQ/params/grif16/template/0",29) == 0 ){
                              if( sscanf(ptr,"<key name=\"prg_ddtm\" type=\"DWORD\">%d</key>", &subsys_prg_ddtm[SUBSYS_HPGE_A]) < 1 ){
                                fprintf(stderr,"can't read key value for /DAQ/params/grif16/template/0/prg_ddtm\n"); ptr=str+1; continue;
                              }
                              fprintf(stdout,"Read in Det type 0 (HPGE A) prg_ddtm as %d\n",subsys_prg_ddtm[SUBSYS_HPGE_A]);
                              while( *(ptr) != '/' ){ ++ptr; } while( *(ptr) != '<' ){ ++ptr; }
                            } else if( strncasecmp(ptr,"<key name=\"prg_ddtm\"", 20) == 0 &&
                            strncmp(path,"/DAQ/params/grif16/template/1",29) == 0 ){
                              if( sscanf(ptr,"<key name=\"prg_ddtm\" type=\"DWORD\">%d</key>", &subsys_prg_ddtm[SUBSYS_HPGE_B]) < 1 ){
                                fprintf(stderr,"can't read key value for /DAQ/params/grif16/template/1/prg_ddtm\n"); ptr=str+1; continue;
                              }
                              fprintf(stdout,"Read in Det type 1 (HPGE B) prg_ddtm as %d\n",subsys_prg_ddtm[SUBSYS_HPGE_B]);
                              while( *(ptr) != '/' ){ ++ptr; } while( *(ptr) != '<' ){ ++ptr; }
                            } else if( strncmp(ptr,"</keyarray>",10) == 0 ){ active = 0; arrayptr = (void *)('\0');
                            if( strncmp(path,"/PPG/Cycles/",12) == 0 ){ odb_ppg_cycle[odb_ppg_cycle_count].length = index+1; }
                          } else if( strncmp(ptr,"<keyarray ",10) == 0 ){
                            if( strncmp(path,"/PPG/Cycles/",12) == 0 ){
                              if( sscanf(ptr,"<keyarray name=\"%s", value) < 1 ){
                                fprintf(stderr,"can't read PPG keyarray entry\n"); ptr=str+1; continue;
                              }
                              if( value[strlen(value)-1]=='\"' ){ value[strlen(value)-1]='\0'; }
                              if( strcmp(value,"PPGcodes") == 0 ){
                                active = 1; arrayptr = (void *)odb_ppg_cycle[odb_ppg_cycle_count].codes; dtype=0;
                              }
                              if( strcmp(value,"durations") == 0 ){
                                active = 1; arrayptr = (void *)odb_ppg_cycle[odb_ppg_cycle_count].durations; dtype=0;
                              }
                              ptr=str+1;
                              continue;
                            }
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

                        // PPG: Now all cycles were unpacked and we identified Current
                        // Copy the relevant ODB PPG pattern into the global variables for this run
                        index=-1;
                        for(i=0; i<MAX_ODB_PPG_CYCLES; i++){
                          if( strncmp(odb_ppg_cycle[i].name,odb_ppg_current,strlen(odb_ppg_current)) == 0 ){
                            index=i;
                            break;
                          }
                        }
                        if(index<0){
                          fprintf(stderr,"Failed to locate Current PPG Cycle in ODB Cycles\n");
                        }else{
                          strcpy(ppg_cycle_name,odb_ppg_cycle[index].name);
                          fprintf(stdout,"PPG cycle for this run is named %s:\n",ppg_cycle_name);
                          ppg_cycle_duration=0;
                          for(i=0; i<odb_ppg_cycle[index].length; i++){
                            ppg_index=-1;
                            for(j=0; j<N_PPG_PATTERNS; j++){ if( (odb_ppg_cycle[index].codes[i] & 0xFFFF) == ppg_patterns[j] ){ ppg_index = j; break; } }
                            if(ppg_index<0){ fprintf(stderr,"unrecognized ppg pattern, 0x%04X\n", (odb_ppg_cycle[index].codes[i] & 0xFFFF)); return(-1); }
                            ppg_cycle_pattern_code[i] = ppg_index;

                            if(odb_ppg_cycle[index].durations[i] == -1){
                              // Infinte duration
                              odb_ppg_cycle[index].length = i+1;
                              ppg_cycle_length = odb_ppg_cycle[index].length;
                              ppg_cycle_duration = 0;
                              ppg_cycles_active = 0;
                            }else{
                              ppg_cycles_active = 1;
                              ppg_cycle_length = odb_ppg_cycle[index].length;
                              ppg_cycle_pattern_duration[i] = (long)odb_ppg_cycle[index].durations[i]*100; // Convert from ODB microseconds to timestamp 10 nanosecond units
                              ppg_cycle_duration += ppg_cycle_pattern_duration[i];
                            }
                            fprintf(stdout,"PPG PATTERN %d: %s (%s) for %10.4f milliseconds (%015ld timestamps)\n", i,
                            ppg_handles[ppg_cycle_pattern_code[i]], ppg_names[ppg_cycle_pattern_code[i]],(double)(ppg_cycle_pattern_duration[i]/100000),ppg_cycle_pattern_duration[i]);
                          }
                          if(ppg_cycle_duration == 0){
                            fprintf(stdout,"PPG cycle duration is infinite, ie. no cycles\n");
                          }else{
                            fprintf(stdout,"PPG cycle duration is %10.4f seconds\n",(double)(ppg_cycle_duration/100000000));
                          }
                          // Set the initial cycle settings
                          ppg_current_pattern = ppg_cycle_pattern_code[0];  // Index of the current PPG cycle pattern for use with the ppg_patterns array
                          ppg_cycle_number = 0;                             // Current cycle number. Cycles counted from zero at beginning of run
                          ppg_cycle_start = 0;                              // Timestamp of the start of the current cycle
                          ppg_cycle_end = ppg_cycle_duration;               // Timestamp of the end of the current cycle
                          ppg_cycle_step = 0;                               // Current pattern number within this cycle. Patterns counted from zero at beginning of cycle
                          ppg_pattern_start = 0;                            // Timestamp of the start of the current pattern
                          ppg_pattern_end = ppg_cycle_pattern_duration[0];  // Timestamp of the end of the current pattern
                          fprintf(stdout,"Cycle %04d, start/finish [%ld/%ld]: step %d, %s, start/finish [%ld/%ld]\n",
                          ppg_cycle_number, ppg_cycle_start, ppg_cycle_end, ppg_cycle_step, ppg_handles[ppg_current_pattern], ppg_pattern_start, ppg_pattern_end);
                        }

                        // arrays typically around 500 entries [one per "chan"] each entry with ...
                        //   daq-address, name, type, gains etc.
                        //
                        gen_derived_odb_tables();

                        return(0);
                      }

                      // original odb arrays were read into {addr_table,chan_name,dtype_table(+gains)}
                      // extract extra details stored in channel names (and record for later)
                      // (these details include crystal/element numbers and polarities)
                      // [use above for subsystem (no longer use datatype to determine subsystems)]
                      extern int read_caen_odb_addresses(int odb_daqsize, unsigned short *addr_table);
                      int gen_derived_odb_tables()
                      {
                        int i, j, tmp, subsys, pos, element, output_type;
                        char sys_name[64], crystal, polarity, type;

                        read_caen_odb_addresses(odb_daqsize, (unsigned short *)addrs);

                        // generate reverse mapping of address to channel number
                        //  (most of this array is undefined and stays at -1)
                        memset(address_chan, 0xFF, sizeof(address_chan)); // set to -1
                        for(i=0; i<MAX_ADDRESS && i<odb_daqsize; i++){
                          address_chan[ (unsigned short)chan_address[i] ] = i;
                        }

                        memset(crystal_table,  0xff, MAX_DAQSIZE*sizeof(int)); // initialise all to -1
                        memset(element_table,  0xff, MAX_DAQSIZE*sizeof(int));
                        memset(polarity_table, 0xff, MAX_DAQSIZE*sizeof(int));
                        memset(subsys_table,   0xff, MAX_DAQSIZE*sizeof(int));
                        for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
                          if( (tmp=sscanf(chan_name[i], "%3c%d%c%c%d%c", sys_name, &pos, &crystal, &polarity, &element, &type)) != 6 ){
                            fprintf(stderr,"can't decode name[%s] decoded %d of 6 items\n", chan_name[i], tmp );
                            continue;
                          }
                          for(j=0; j<MAX_ODB_SUBSYS; j++){
                            if( strncmp(sys_name, odb_handle[j], 3) == 0 ){ subsys = j; break; }
                          }
                          if( j == MAX_ODB_SUBSYS ){ subsys = j-1; // use final entry: "unknown"
                          fprintf(stderr,"Unknown subsystem[%s] in %s\n", sys_name, chan_name[i]);
                        }

                        // Mention bad detector types (no longer relied on for subsystem id)
                        if( dtype_table[i] < 0 || dtype_table[i] >= 16 ){
                          fprintf(stderr,"bad datatype[%d] at table position %d\n", dtype_table[i], i);
                        }

                        // Some detector elements have more than one output (HPGe A and B)
                        // 1 is A, 0 is B, -1 is X or unknown
                        output_type = type=='A' ? 1 : (type=='B' ? 0 : -1);
                        if(subsys == ODBHANDLE_ZDS){ if(type=='X'){ output_type = 1; } } // Older ZDS convention

                        // Polarity: 1 is N, 0 is P or T or S, -1 is anything else
                        if(        polarity == 'N' ){ polarity_table[i] = 1;
                        } else if( polarity == 'P' ){ polarity_table[i] = 0;
                        } else if( polarity == 'T' ){ polarity_table[i] = 0; // TAC signal
                        } else if( polarity == 'S' ){ polarity_table[i] = 1; // ARIES Standard Ouput signal
                        } else if( polarity == 'F' ){ polarity_table[i] = 0; // ARIES Fast Output signal
                        } else if( polarity == 'X' ){ polarity_table[i] = 0; // XXX type
                        } else { fprintf(stderr,"unknown polarity[=%c] in %s\n", polarity, chan_name[i]); }

                        // Record crystal and element numbers [** Naming schemes are subsystem-dependant **]
                        switch(subsys){
                          case ODBHANDLE_LBL: case ODBHANDLE_LBS: // LaBr,Paces, Ares and Zds
                          case ODBHANDLE_LBT: case ODBHANDLE_ART:
                          case ODBHANDLE_DAL: case ODBHANDLE_DAT:
                          case ODBHANDLE_PAC: case ODBHANDLE_ZDS: case ODBHANDLE_DSW:
                          crystal_table[i] = pos;
                          if(        crystal == 'A' ){ element_table[i] = 1;
                          } else if( crystal == 'B' ){ element_table[i] = 2;
                          } else if( crystal == 'C' ){ element_table[i] = 3;
                          } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3, ZDS, ART, LBT
                          } else {
                            fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
                          } break;
                          case ODBHANDLE_RCS:
                          crystal_table[i] = pos;
                          element_table[i] = reorder_rcmp_strips[pos][polarity_table[i]][element];
                          break;
                          case ODBHANDLE_GRG: case ODBHANDLE_GRS:
                          element_table[i] = element;
                          pos -= 1; pos *=4;
                          if(        crystal == 'B' ){ crystal_table[i] = pos;
                          } else if( crystal == 'G' ){ crystal_table[i] = pos+1;
                          } else if( crystal == 'R' ){ crystal_table[i] = pos+2;
                          } else if( crystal == 'W' ){ crystal_table[i] = pos+3;
                          } else if( crystal == 'X' ){ crystal_table[i] = -1; // crystal undefined
                          } else {
                            fprintf(stderr,"unknown crystal[=%c] in %s\n", crystal, chan_name[i]);
                          } break;
                          default: break;
                        }

                        // set full subsystem id (including polarity/output-type etc)
                        switch(subsys){
                          case ODBHANDLE_GRS: subsys_table[i] = SUBSYS_BGO;       break;
                          case ODBHANDLE_SEP: subsys_table[i] = SUBSYS_SCEPTAR;   break;
                          case ODBHANDLE_PAC: subsys_table[i] = SUBSYS_PACES;     break;
                          case ODBHANDLE_LBS: subsys_table[i] = SUBSYS_LABR_BGO;  break;
                          case ODBHANDLE_LBL: subsys_table[i] = SUBSYS_LABR_L;    break;
                          case ODBHANDLE_DAL: subsys_table[i] = SUBSYS_LABR_L;    break;
                          case ODBHANDLE_DSC: subsys_table[i] = SUBSYS_DESCANT;   break;
                          case ODBHANDLE_RCS: subsys_table[i] = SUBSYS_RCMP;      break;
                          case ODBHANDLE_DSW: subsys_table[i] = SUBSYS_DESWALL;  break;
                          case ODBHANDLE_DSG: subsys_table[i] = SUBSYS_DSG;  break;
                          case ODBHANDLE_GRG: subsys_table[i] = (output_type == 1) ? SUBSYS_HPGE_A :SUBSYS_HPGE_B; break;
                          case ODBHANDLE_ZDS: subsys_table[i] = (output_type == 1) ? SUBSYS_ZDS_A  :SUBSYS_ZDS_B;  break;
                          case ODBHANDLE_ART: subsys_table[i] = (polarity_table[i] == 1) ? SUBSYS_ARIES_A:SUBSYS_ARIES_B;break;
                          case ODBHANDLE_XXX: subsys_table[i] = SUBSYS_IGNORE;    break;
                          case ODBHANDLE_UNK: subsys_table[i] = SUBSYS_UNKNOWN;   break;
                          case ODBHANDLE_DAT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
                          }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
                          break;
                          case ODBHANDLE_LBT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
                          }else if(crystal_table[i]>8){ subsys_table[i] = SUBSYS_TAC_ART;
                          }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
                          break;
                        }
                      }
                      memset(subsys_initialized, 0, sizeof(int)*MAX_SUBSYS );

                      return(0);
                    }
