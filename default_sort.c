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
static short *chan_address = addr_table;
static int subsys_initialized[MAX_SUBSYS];
extern Grif_event grif_event[MAX_COINC_EVENTS];

// Default sort function declarations
extern int init_time_diff_gates(Config *cfg);
extern int init_chan_histos(Config *cfg);
extern int init_histos(Config *cfg, int subsystem);
extern int fill_chan_histos(Grif_event *ptr);
extern int fill_singles_histos(Grif_event *ptr);
extern int fill_coinc_histos(int win_idx, int frag_idx);

// odb tables need to be transferred into config, which is saved with histos
int init_default_histos(Config *cfg, Sort_status *arg)
{
   Cal_coeff *cal;
   int i, j;

   // Initialize all pileup parameters to unset values
   for(i=0; i<odb_daqsize; i++){
     for(j=0; j<7; j++){
       pileupk1[i][j] = pileupk2[i][j] = pileupE1[i][j] = -1;
     }
   }

   cfg->odb_daqsize = odb_daqsize;
   for(i=0; i<odb_daqsize; i++){ // update config with odb info
     edit_calibration(cfg, chan_name[i], offsets[i], gains[i], quads[i], pileupk1[i], pileupk2[i], pileupE1[i],
                      chan_address[i],  dtype_table[i], arg->cal_overwrite );
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

       // Pileup parameters do not exist in the MIDAS ODB so must always be copied from the config
       for(j=0; j<7; j++){
         pileupk1[i][j] = (isnan(cal->pileupk1[j])) ? 0.0 : cal->pileupk1[j];
         pileupk2[i][j] = (isnan(cal->pileupk2[j])) ? 0.0 : cal->pileupk2[j];
         pileupE1[i][j] = (isnan(cal->pileupE1[j])) ? 0.0 : cal->pileupE1[j];
       }
     }

   init_time_diff_gates(cfg);
   init_chan_histos(cfg);
   init_histos(cfg, SUBSYS_HPGE_A); // always create Ge histos

   return(0);
}

//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#######################################################################

float spread(int val){ return( val + rand()/(1.0*RAND_MAX) ); }
int GetIDfromAddress(unsigned short addr){ // address must be an unsigned short
  return(address_chan[addr]);
}

int init_time_diff_gates(Config *cfg){
  int i,j,k;
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
  bgo_window_min = addback_window_min = rcmp_fb_window_min = lbl_tac_window_min = art_tac_window_min = desw_beta_window_min = 0;
  bgo_window_max = 20;
  addback_window_max = 20;
  rcmp_fb_window_max = 10;
  lbl_tac_window_max = 25;
  art_tac_window_max = 25;
  desw_beta_window_max = 80;

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
    }else if(strncmp(tmp,"presort_time_diff_addback",25) == 0){
      addback_window_min = global->min; addback_window_max = global->max;
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
    }
  }// end of for(i=0; i<cfg->nglobal; i++){
  return(0);
}

// this is the first function to be called on processing an event -
// before any presort/singles/coinc-sorting
int apply_gains(Grif_event *ptr)
{
  int tac_ts_offset[12] = {60,60,60,60,60,60,60,60,60,60,60,60}; // From Dec 2024
  int caen_ts_offset = -60; // this value (-60) aligns the timestamps of HPGe with ZDS(CAEN)
   float energy,psd;
   int chan = ptr->chan;

    // Protect against invalid channel numbers
    if( chan < 0 || chan >= odb_daqsize ){
       fprintf(stderr,"unpack_event: ignored event in chan:%d [0x%04x]\n", chan, ptr->address );
       return(-1);
    }

   // Calculate the energy and calibrated energies
   ptr->energy = energy = ( ptr->integ == 0 ) ? ptr->q : spread(ptr->q)/ptr->integ;
   ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);

      // NOBODY CURRENTLY USES e2,e3,e4 ...
      if( ptr->integ2 != 0 ){
         energy = ptr->energy2 = spread(ptr->q2)/ptr->integ2;
         ptr->e2cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
      }
      if( ptr->integ3 != 0 ){
         energy = ptr->energy3 = spread(ptr->q3)/ptr->integ3;
         ptr->e3cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
      }
      if( ptr->integ4 != 0 ){
         energy = ptr->energy4 = spread(ptr->q4)/ptr->integ4;
         ptr->e4cal = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
      }

      // Assign the subsys type
         if( (ptr->subsys = subsys_table[chan]) == -1 ){ return(-1); }
         if( subsys_initialized[ptr->subsys] == 0 ){
            //init_histos(configs[1], ptr->subsys);
	      init_histos(NULL, ptr->subsys);
         }

   // HPGe pileup
   if( ptr->subsys == SUBSYS_HPGE_A){
     ptr->psd = 14; // Pileup class - default value of 12 for all HPGe events
     if(ptr->pileup==1 && ptr->nhit ==1){
       // Single hit events
       // no pileup, this is the most common type of HPGe event
       ptr->psd = 1; // Pileup class, default for single hit events
     }
   }

   // The TAC module produces its output signal around 2 microseconds later
   // than the start and stop detector signals are processed.
   if( ptr->subsys == SUBSYS_LABR_T){
    ptr->ts -= tac_ts_offset[crystal_table[ptr->chan]-1]; // Subtract 2 microseconds from TAC timestamps plus a more precise offset
   }

   // DESCANT detectors
   // use psd for Pulse Shape Discrimination provides a distinction between neutron and gamma events
   //if( ptr->subsys == SUBSYS_DESCANT || ptr->subsys == SUBSYS_DESWALL){
   if( ptr->subsys == SUBSYS_DESWALL){
     //ptr->ts -= caen_ts_offset; // Subtract from CAEN timestamps to align coincidences
     psd = ( ptr->q != 0 ) ? (spread(ptr->cc_short) / ptr->q) : 0;
     ptr->psd = (int)(psd*1000.0); // psd = long integration divided by short integration
   }


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
  Grif_event *alt2, *alt, *ptr = &grif_event[frag_idx];
  float desw_median_distance = 1681.8328; // descant wall median source-to-detector distance in mm
  int i, j, dt, dt13, tof;
  float q1,q2,q12,k1,k2,k12,e1,e2,e12,m,c;
  int chan,found,pos;
  float energy,correction;
  float correction12, correction23;

  // Assign chan local variable and check it is a valid channel number
  if( (chan=ptr->chan)<0 || ptr->chan >= odb_daqsize ){
    fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
    return(-1);
  }
  i = frag_idx; ptr->fold = 1;
  while( i != end_idx ){ // need at least two events in window
    if( ++i >=  MAX_COINC_EVENTS ){ i=0; } alt = &grif_event[i]; // WRAP
    if( alt->chan<0 || alt->chan >= odb_daqsize ){
      fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
      return(-1);
    }

    // Determine fold
    if( alt->subsys == ptr->subsys ){ ++ptr->fold; }

    // Determine absolute time difference between timestamps
    dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }

    // SubSystem-specific pre-processing
    switch(ptr->subsys){
      case SUBSYS_HPGE_A:

      // HPGe pile-up corrections
      // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
      // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
      // First assign the pileup class type, then correct the energies
      if(alt->subsys == SUBSYS_HPGE_A && alt->chan == ptr->chan){
        if(ptr->pileup==1 && ptr->nhit ==1){
          // no pileup, this is the most common type of HPGe event
          ptr->psd = 1; // Pileup class
        }else if(ptr->pileup==0){
          // pileup error
          ptr->psd = 0; // Pileup class, error
        }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==2 && alt->nhit==1)){
          ptr->psd = alt->psd = 9; // Pileup class, error for 2Hits
          ptr->ts_int = alt->ts_int = dt;
          if(ptr->q>0 && ptr->integ>0 && ptr->q2>0 && ptr->integ2>0 && alt->q>0 && alt->integ>0){
            // 2 Hit, Type A
            // The (ptr) fragement is the first Hit of a two Hit pile-up event.
            // Assign the pileup class numbers to the two hits and calculate the time difference between them
            ptr->psd = 3; // Pileup class, 1st of 2Hits
            alt->psd = 4; // Pileup class, 2nd of 2Hits
            ptr->ts_int = alt->ts_int = dt;

          }else{
            // 2 Hit, Type B
            // 2Hit pileup where 2nd Hit integration region starts after 1st Hit integration ends but before 2nd Hit CFD has completed
            // Assign the pileup class numbers to the two hits and calculate the time difference between them
            ptr->psd = 7; // Pileup class, 1st of 2Hits
            alt->psd = 8; // Pileup class, 2nd of 2Hits
            ptr->ts_int = alt->ts_int = dt;

          }
        }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==1 && alt->nhit==1)){
          // 2 Hit, Type C
          // 2Hit pileup where 2nd Hit integration region starts after 1st Hit integration and 2nd Hit CFD have ended
          // Assign the pileup class numbers to the two hits and calculate the time difference between them
          ptr->psd = 5; // Pileup class, 1st of 2Hits
          alt->psd = 6; // Pileup class, 2nd of 2Hits
          ptr->ts_int = alt->ts_int = dt; // Save the time difference between pileup hits into both hits

        }else if((ptr->pileup==1 && ptr->nhit==3) && (alt->pileup==2 && alt->nhit==2)){ // 3Hit pileup
          ptr->psd = alt->psd = 13; // Pileup class, error for 3Hits
          if(ptr->q>0 && ptr->integ>0 && ptr->q2>0 && ptr->integ2>0 && alt->q>1 && alt->integ>0 && alt->q2>0 && alt->integ2>0){
            /*
            found=0;
            //  if(ptr->ecal > 1160 && ptr->ecal < 1180 && alt->ecal > 1325 && alt->ecal < 1350){
            fprintf(stdout,"Found a pileup 3 hit group for chan %d with dt %d\n",ptr->chan,dt);
            fprintf(stdout,"ptr: %ld %d, PU=%d, nhits=%d, q: %d %d %d %d, k: %d %d %d %d, ecal: %d %d %d %d, %lf %lf %lf\n",ptr->ts,ptr->cfd,ptr->pileup,ptr->nhit,ptr->q,ptr->q2,ptr->q3,ptr->q4,ptr->integ,ptr->integ2,ptr->integ3,ptr->integ4,ptr->ecal,ptr->e2cal,ptr->e3cal,ptr->e4cal,offsets[ptr->chan],gains[ptr->chan],quads[ptr->chan]);
            fprintf(stdout,"alt: %ld %d, PU=%d, nhits=%d, q: %d %d %d %d, k: %d %d %d %d, ecal: %d %d %d %d, %lf %lf %lf\n",alt->ts,alt->cfd,alt->pileup,alt->nhit,alt->q,alt->q2,alt->q3,alt->q4,alt->integ,alt->integ2,alt->integ3,alt->integ4,alt->ecal,alt->e2cal,alt->e3cal,alt->e4cal,offsets[alt->chan],gains[alt->chan],quads[alt->chan]);
            fprintf(stdout,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",ptr->q,ptr->q2,alt->q,ptr->integ,ptr->integ2,alt->integ,ptr->energy,ptr->energy2,alt->energy,ptr->ecal,ptr->e2cal,alt->ecal);
            found=1;
            //    }
            */
            j=i+1;
            while( j != end_idx ){ // need to find the third events in window associated with this channel
              if( ++j >=  MAX_COINC_EVENTS ){ break; } alt2 = &grif_event[j]; // WRAP

              if(alt2->chan == ptr->chan){ // It must also be a HPGe if the channel number is the same
                /*
                fprintf(stdout,"alt2: %ld %d, PU=%d, nhits=%d, q: %d %d %d %d, k: %d %d %d %d, ecal: %d %d %d %d, %lf %lf %lf\n",alt2->ts,alt2->cfd,alt2->pileup,alt2->nhit,alt2->q,alt2->q2,alt2->q3,alt2->q4,alt2->integ,alt2->integ2,alt2->integ3,alt2->integ4,alt2->ecal,alt2->e2cal,alt2->e3cal,alt2->e4cal,offsets[alt2->chan],gains[alt2->chan],quads[alt2->chan]);
                */
                if(alt2->pileup==3 && alt2->nhit==1){
                  alt2->psd = 12; // Pileup class
                  if(alt2->q>1 && alt2->integ>0){

                    // Determine absolute time difference between timestamps for Hit 1 and 3
                    dt13 = ptr->ts - alt2->ts; if( dt13 < 0 ){ dt13 = -1*dt13; }


                    // Three-Hit pile-up case...
                    // there are two types depending on the relative timing of the third Hit...
                    // Triple pileup case A ... in which the 3rd pulse occurs more than L samples after the first
                    //                          there are 5 regions, 3 of which are not piled up (1 per pulse)
                    //                      ____                ____
                    //    |   |        |  /|    |\  |      |  /|    |\  |      |          :
                    //    |   |        | / |    | \ |      | / |    | \ |      |          :
                    //    |   *________|/__|____|  \|______|/__|____|  \|______|          :
                    //    |  /                   \                  \          \          :
                    //    | /   K1           K12  \    K2        K23 \     K3   \         :
                    //  __*/                       \_____             \______    \________:
                    //    0            S                   X
                    //
                    // ------------------------------------------------------------------------------------------
                    // Triple pileup case B ... in which the 3rd pulse occurs less than L samples after the first
                    //                          again 5 regions, only 2 of which are not piled up (first and last pulse)
                    //                          There is no region to obtain the height of pulse 2
                    //                          so the event contains K12, the sum of pulse 1+2, in place of pulseheight2
                    //                                    ________
                    //    |   |        |   |         |  /|        |\  |      |   |      |   :
                    //    |   |        |   |         | / |        | \ |      |   |      |   :
                    //    |   |        |   |_________|/__|________|  \|______|   |      |   :
                    //    |   |        |  /|          :           |\  |      |\  |      |   :
                    //    |   |        | / |          :           | \ |      | \ |      |   :
                    //    |   *________|/__|__________:___________|  \|______|__\|______|   :
                    //    |  /                        :           \                     \   :
                    //    | /     K1            K12   :    K123    \     K23        K3   \  :
                    //  __*/                          :             \_____                \_:
                    //    0            S              X           L
                    //

                    // The Differencitation period of the HPGe Pulse Height evaluation is L = 5000ns.
                    if(dt13>500){
                      // Triple pileup case A ... in which the 3rd pulse occurs more than L samples after the first
                      //                          there are 5 regions, 3 of which are not piled up (1 per pulse)
                      correction23 = (alt->q/alt->integ)-((alt->q2/alt->integ2)-(alt2->q/alt2->integ));
                      correction12 = (ptr->q/ptr->integ)-((ptr->q2/ptr->integ2)-(alt->q/alt->integ)-correction23);
                      // Hit 1
                      ptr->psd = 10; // Pileup class
                      ptr->energy = energy = (spread(ptr->q)/ptr->integ) + correction12;
                      ptr->ecal=ptr->esum = offsets[ptr->chan]+energy*(gains[ptr->chan]+energy*quads[ptr->chan]);
                      // Hit 2
                      alt->ts_int = dt;
                      alt->psd = 11; // Pileup class
                      alt->energy = energy = (spread(alt->q)/alt->integ) - correction12 + correction23;
                      alt->ecal=alt->esum = offsets[alt->chan]+energy*(gains[alt->chan]+energy*quads[alt->chan]);
                      // Hit 3
                      alt2->ts_int = dt13;
                      alt2->psd = 12; // Pileup class
                      alt2->energy = energy = (spread(alt2->q)/alt2->integ) - correction23;
                      alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
                    }else{
                      // Triple pileup case B ... in which the 3rd pulse occurs less than L samples after the first
                      //                          again 5 regions, only 2 of which are not piled up (first and last pulse)
                      //                          There is no region to obtain the height of pulse 2
                      //                          so the event contains K12, the sum of pulse 1+2, in place of pulseheight2
                      correction23 = (alt->q/alt->integ)-((alt->q2/alt->integ2)-(alt2->q/alt2->integ));
                      correction12 = (ptr->q/ptr->integ)-((ptr->q2/ptr->integ2)-(alt->q/alt->integ)-correction23);
                      // Hit 1
                      ptr->psd = 10; // Pileup class
                      ptr->energy = energy = (spread(ptr->q)/ptr->integ) + correction12;
                      ptr->ecal=ptr->esum = offsets[ptr->chan]+ptr->energy*(gains[ptr->chan]+ptr->energy*quads[ptr->chan]);
                      // Hit 2
                      alt->ts_int = dt;
                      alt->psd = 11; // Pileup class
                      alt->energy = energy = (spread(alt->q)/alt->integ) - correction12 + correction23;
                      alt->ecal=alt->esum = offsets[alt->chan]+energy*(gains[alt->chan]+energy*quads[alt->chan]);
                      // Hit 3
                      alt2->ts_int = dt13;
                      alt2->psd = 12; // Pileup class
                      alt2->energy = energy = (spread(alt2->q)/alt2->integ) - correction23;
                      alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
                    }
                    /*
                    fprintf(stdout,"Complete 3Hit PU event, dt13=%d: %d,%d,%d,%d,%d, %d,%d,%d,%d,%d, %d,%d,%d,%d,%d, %d,%d,%d,%d,%d\n\n",dt13,ptr->q,ptr->q2,alt->q,alt->q2,alt2->q,ptr->integ,ptr->integ2,alt->integ,alt->integ2,alt2->integ,ptr->energy,ptr->energy2,alt->energy,alt->energy2,alt2->energy,ptr->ecal,ptr->e2cal,alt->ecal,alt->e2cal,alt2->ecal);
                    */
                    break; // Break the while if we found the third Hit
                  }
                }
              }
            } // end of while for triple coincidence
          }
        } // end of 3Hit pileup type assignments

        // Now apply hit-specific energy corrections
        if(pileupk1[chan][0] != 1){
          if(ptr->psd>=3 && ptr->psd<=8){ // 2-Hit pileup events
            // Apply the k1 dependant correction to the energy of the first hit
            // It was already checked that chan for ptr and alt are the same for pileup events
            pos  = crystal_table[ptr->chan];
            k1 = ptr->integ;
            ptr->ecal=ptr->esum = ptr->ecal*( pileupk1[chan][0]+(k1*pileupk1[chan][1])+(k1*k1*pileupk1[chan][2])+(k1*k1*k1*pileupk1[chan][3])
            +(k1*k1*k1*k1*pileupk1[chan][4])+(k1*k1*k1*k1*k1*pileupk1[chan][5])+(k1*k1*k1*k1*k1*k1*pileupk1[chan][6]));
            alt->e4cal=ptr->ecal; // Remember the ecal of the first Hit in this second Hit

            // Apply the E1 and k2 dependant offset correction to the energy of the second hit
            // Apply the k2 dependant correction to the energy of the second hit
            k2 = alt->integ;
            correction = ptr->ecal*( pileupE1[chan][0]+(k2*pileupE1[chan][1])+(k2*k2*pileupE1[chan][2])+(k2*k2*k2*pileupE1[chan][3])
            +(k2*k2*k2*k2*pileupE1[chan][4])+(k2*k2*k2*k2*k2*pileupE1[chan][5])+(k2*k2*k2*k2*k2*k2*pileupE1[chan][6]));
            alt->ecal=alt->esum = (alt->ecal*( pileupk2[chan][0]+(k2*pileupk2[chan][1])+(k2*k2*pileupk2[chan][2])+(k2*k2*k2*pileupk2[chan][3])
            +(k2*k2*k2*k2*pileupk2[chan][4])+(k2*k2*k2*k2*k2*pileupk2[chan][5])+(k2*k2*k2*k2*k2*k2*pileupk2[chan][6])))+correction;
          }
        }
      } // end of if(alt->subsys == SUBSYS_HPGE_A && alt->chan == ptr->chan)

      // BGO suppression of HPGe
      if( (dt >= bgo_window_min && dt <= bgo_window_max) && alt->subsys == SUBSYS_BGO && !ptr->suppress ){
        // could alternatively use crystal numbers rather than clover#
        //    (don't currently have this for BGO)
        if( crystal_table[ptr->chan]/16 == crystal_table[alt->chan]/16 ){ ptr->suppress = 1; }
      }
      // Germanium addback -
      //    earliest fragment has the sum energy, others are marked -1
      // Remember the other crystal channel number in ab_alt_chan for use in Compton Polarimetry
      if( (dt >= addback_window_min && dt <= addback_window_max) && alt->subsys == SUBSYS_HPGE_A ){
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
      if( (dt >= rcmp_fb_window_min && dt <= rcmp_fb_window_max) && alt->subsys == SUBSYS_RCMP && (ptr->ecal>0 && ptr->ecal<32768)){
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
      // In the TAC event we save the tile chan as ab_alt_chan, and the tile energy as e4cal.
      // So later in the main coincidence loop we only need to examine LBL and TAC.
      if( (dt >= art_tac_window_min && dt <= art_tac_window_max) && alt->subsys == SUBSYS_ARIES_A && crystal_table[ptr->chan] == 8){
        ptr->ab_alt_chan = alt->chan; ptr->e4cal = alt->ecal;
      }
      // For TAC01-07 we have a LBL-LBL coincidence
      // Here save the LBL Id number and the LBL energy in the TAC event
      // Save LBL channel number into ptr->q2 or q3 or q4
      // Save LBL energy ecal into TAC ptr-ecal2 or ecal3 or ecal4
      if( (dt >= lbl_tac_window_min && dt <= lbl_tac_window_max) && alt->subsys == SUBSYS_LABR_L && crystal_table[ptr->chan] < 8){
        if(ptr->e2cal<1){
          ptr->q2 = alt->chan; ptr->e2cal = alt->ecal;
        }else if(ptr->e3cal<1){
          ptr->q3 = alt->chan; ptr->e3cal = alt->ecal;
        }else{
          ptr->q4 = alt->chan; ptr->e4cal = alt->ecal; // If this is set then we have LBL multiplicity >2 for this TAC
        }
      }
      break;
      case SUBSYS_ZDS_B: // CAEN Zds
      if(alt->subsys == SUBSYS_DESWALL){
        if(dt >= desw_beta_window_min && dt <= desw_beta_window_max){
          // Calculate time-of-flight and correct it for this DESCANT detector distance
          tof = (spread(abs(ptr->cfd - alt->cfd))*2.0) + 100; //if( tof < 0 ){ tof = -1*tof; }
          //  fprintf(stdout,"tof: %d - %d = %f\n",ptr->cfd, alt->cfd, tof);
          alt->energy4 = (int)(tof); // Time of flight
          alt->e4cal = (int)(spread(tof) * DSW_tof_corr_factor[crystal_table[alt->chan]-1]); // Corrected Time of Flight
        }
      }
      break;
      default: break; // Unrecognized or unprocessed subsys type
    }// end of switch
  }// end of while
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
  // {(void **)&labr_tac_xtal,"TAC_LBL_ART_vs_LBL_Num", "TAC_ART_LBL_LBL_Xtal",       SUBSYS_LABR_T,   16,  E_2D_SPECLEN},
   {(void **)&art_tac_xtal, "TAC_LBL_ART_vs_ART_Num", "TAC_ART_LBL_ART_Xtal",       SUBSYS_ARIES_A,  80,  E_2D_SPECLEN},
   {(void **)&desw_e_xtal,  "DESWall_En_DetNum",    "DSW_En_Xtal",                 SUBSYS_DESWALL, 64,  E_2D_SPECLEN},
   {(void **)&desw_tof_xtal,"DESWall_TOF_DetNum",   "DSW_TOF_Xtal",                SUBSYS_DESWALL, 64,  E_2D_SPECLEN},
   {(void **)&desw_psd_e,   "DESWall_PSD_En",       "DES_Wall_PSD_En",             SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
   {(void **)&desw_psd_tof, "DES_Wall_PSD_TOF",      "DES_Wall_PSD_TOF",            SUBSYS_DESWALL, E_2D_SPECLEN, E_2D_SPECLEN},
   {(void **) rcmp_strips,  "RCS%02d_E_strips",         "",                         SUBSYS_RCMP,   2*N_RCMP_STRIPS, E_2D_RCMP_SPECLEN, N_RCMP_POS},
   {NULL,                   "Hits_and_Sums/Pileup",     "",                         },
   {(void **)&ge_pu_class,  "Pile-up_class",           "",                          SUBSYS_HPGE_A,         64},
   {(void **) ge_sum_class,    "", ge_pu_class_sum_titles[0],                       SUBSYS_HPGE_A,  E_SPECLEN,   0, N_PU_CLASSES},
   {(void **) ge_e_vs_k_class, "", ge_pu_class_2d_titles[0],                        SUBSYS_HPGE_A,       2048, 512, N_PU_CLASSES},
   {(void **)&ge_pu_type,   "Pile-up_type",            "",                          SUBSYS_HPGE_A,         64},
   {(void **)&ge_nhits_type,"nhits_type",              "",                          SUBSYS_HPGE_A,         64},
   {(void **) ge_1hit,      "Ge%02d_Single_hit",       "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
   {(void **) ge_2hit,      "Ge%02d_2_hit_pileup",     "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
   {(void **) ge_3hit,      "Ge%02d_3_hit_pileup",     "",                          SUBSYS_HPGE_A,  E_SPECLEN, 0, 64},
/*
//  These definitions of pileup spectra need to be converted to the new method.
// Filling of these spectra is currently disabled.

  sprintf(title,  "Pile-up_3Hits_detlaT_1_2"); sprintf(handle, "PU_dt12");
  ge_pu_dt12 = H1_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title,  "Pile-up_3Hits_detlaT_1_3"); sprintf(handle, "PU_dt13");
  ge_pu_dt13 = H1_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);

close_folder(cfg);
*/
   {NULL,                   "Hits_and_Sums/Pileup-corrections",     "",                         },
   {(void **) ge_e_vs_k_2hit_first,      "Ge%02d_E_vs_k_1st_of_2hit",              "", SUBSYS_HPGE_A,  2048, 512, 64},
   {(void **) ge_e_vs_k_2hit_second,     "Ge%02d_E_vs_k_2nd_of_2hit",              "", SUBSYS_HPGE_A,  2048, 512, 64},
   {(void **) ge_PU2_e2_v_k_gatedxrays,  "Ge%02d_PU2_E2_vs_k2_E1gated_on_Xrays",   "", SUBSYS_HPGE_A,  2048, 512, 64},
   {(void **) ge_PU2_e2_v_k_gated1408,   "Ge%02d_PU2_E2_vs_k2_E1gated_on_1408keV", "", SUBSYS_HPGE_A,   512, 512, 64},
   // Coinc
   {NULL,                  "Hits_and_Sums/Delta_t"," "},
   {(void **) dt_hist,     "",                   dt_handles[0], SUBSYS_HPGE_A,  DT_SPEC_LENGTH, 0, N_DT }, // leave subsys as GE -> all always defined
   {(void **) dt_tacs_hist,"dt_labr_tac%d",      "",  SUBSYS_LABR_T,  DT_SPEC_LENGTH, 0, N_TACS },
   {NULL,                  "Coinc/Coinc",        ""},
   {(void **)&gg_ab,       "Addback_GG",         "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
   {(void **)&gg,          "GG",                 "",  SUBSYS_HPGE_A,  E_2D_SPECLEN, SYMMETERIZE},
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
   {NULL,                  "Coinc/Hits",         ""},
   {(void **)&gg_hit,      "GeGeHit",            "",  SUBSYS_HPGE_A,  64,  64},
   {(void **)&bgobgo_hit,  "BgoBgoHit",          "",  SUBSYS_BGO,    512, 512},
   {(void **)&gea_hit,     "GeAriesHit",         "",  SUBSYS_HPGE_A,  80,  64},
   {(void **)&lba_hit,     "LaBrAriesHit",       "",  SUBSYS_LABR_L,  16,  80},
   {(void **)&aa_hit,      "AriesAriesHit",      "",  SUBSYS_ARIES_A, 80,  80},
   {(void **)&dsw_hit,     "DSWDSWHit",          "",  SUBSYS_DESWALL,64,  64},
   {(void **) rcmp_hit,    "RCS%d_PN_hit",       "",  SUBSYS_RCMP, N_RCMP_STRIPS,     N_RCMP_STRIPS,     N_RCMP_POS},
   {(void **) rcmp_fb,     "RCS%d__Front_Back",  "",  SUBSYS_RCMP, E_2D_RCMP_SPECLEN, E_2D_RCMP_SPECLEN, N_RCMP_POS},
   {NULL,                  "Ang_Corr/GG_Ang_Corr"     , ""},
   {(void **) gg_angcor_110,"Ge-Ge_110mm_angular_bin%d","", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  SYMMETERIZE, N_GE_ANG_CORR},
   {(void **) gg_angcor_145,"Ge-Ge_145mm_angular_bin%d","", SUBSYS_HPGE_A,  GE_ANGCOR_SPECLEN,  SYMMETERIZE, N_GE_ANG_CORR},
   {NULL,                   "Ang_Corr/DSW_DSW_Ang_Corr",""},
   {(void **) dsw_angcor,   "DSW-DSW_angular_bin%d",    "", SUBSYS_DESWALL,DSW_ANGCOR_SPECLEN, SYMMETERIZE, N_DSW_DSW_ANG_CORR },
   {NULL,                   "Ang_Corr/GG_ART_Ang_Corr", ""},
   {(void **) ge_art_angcor,"Ge-ART_angular_bin%d",     "", SUBSYS_ARIES_A,  GE_ANGCOR_SPECLEN, GE_ANGCOR_SPECLEN, N_GRG_ART_ANG_CORR },
   {NULL,                   "Fast-Timing/LBL_Walk",     ""},
   {(void **) tac_labr_CompWalk,"TAC_%d_CompWalk",      "", SUBSYS_LABR_T, 2048, 2048},
   {NULL,                   "Fast-Timing/ART_TACs",     ""},
   {(void **)&tac_aries_lbl_sum, "TAC_ART_LBL_LBLSUM",     "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
   {(void **)&tac_aries_art_sum, "TAC_ART_LBL_ARTSUM",     "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
   {(void **)&aries_tac,         "TAC-ARIES-LaBr3-1275keV","", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
   {(void **)&aries_tac_Egate,   "TAC-ARTE-LaBr3-1275keV", "", SUBSYS_ARIES_A, E_TAC_SPECLEN  },
   {(void **)&aries_tac_artEn,   "TAC-ARIES-Energy",       "", SUBSYS_ARIES_A, E_SPECLEN      },
   {(void **) tac_aries_lbl,    "TAC_ART_LBL%d",           "", SUBSYS_ARIES_A, E_TAC_SPECLEN, E_TAC_SPECLEN, N_LABR  },
   {(void **) tac_aries_art,    "TAC_LBL_ART%d",           "", SUBSYS_ARIES_A, E_TAC_SPECLEN, E_TAC_SPECLEN, N_ARIES }
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
   if( subsystem == SUBSYS_LABR_T ){ // TAC coincidence pair spectra
      open_folder(cfg, "Fast-Timing");
      open_folder(cfg, "LBL_TACs");
      k=0; memset(tac_labr_hist_index, -1, N_LABR*N_LABR*sizeof(int));
      for(i=1; i<=N_LABR; i++){
         for(j=(i+1); j<=N_LABR; j++){
            tac_labr_hist_index[i-1][j-1] = k++;
            sprintf(tmp,"TAC_%d_%d", i, j);
            tac_labr_hist[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
         }
      }
      sprintf(tmp,"TAC_%d_%d", 2, 1); // Add additional histogram (2_1) needed for Compton Walk corrections
      tac_labr_hist[++k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPECLEN, 0, E_TAC_SPECLEN);
      tac_labr_hist_index[1][0] = k;
      close_folder(cfg);
      close_folder(cfg);
   }
      // fill in subsys EvsE and Dt table pointers [** [X][<=Y]
      subsys_e_vs_e[SUBSYS_HPGE_A ][SUBSYS_HPGE_A ] = gg;
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
      subsys_dt[SUBSYS_LABR_L ][SUBSYS_LABR_T  ] = dt_hist[16];
      subsys_dt[SUBSYS_RCMP   ][SUBSYS_RCMP    ] = dt_hist[ 9];
      subsys_dt[SUBSYS_ARIES_A][SUBSYS_ARIES_A ] = dt_hist[13];
      subsys_dt[SUBSYS_ARIES_A][SUBSYS_LABR_T  ] = dt_hist[14];
      subsys_dt[SUBSYS_ZDS_A  ][SUBSYS_LABR_T  ] = dt_hist[15];
      subsys_dt[SUBSYS_ZDS_A  ][SUBSYS_ZDS_B   ] = dt_hist[22];
      subsys_dt[SUBSYS_DESWALL][SUBSYS_DESWALL ] = dt_hist[18];
      subsys_dt[SUBSYS_DESWALL][SUBSYS_ARIES_A ] = dt_hist[20];
      subsys_dt[SUBSYS_DESWALL][SUBSYS_ZDS_A   ] = dt_hist[21];
   return(0);
}

int fill_singles_histos(Grif_event *ptr)
{
  int i, j, dt, pu, nhits, chan, pos, sys, elem, clover, c1,c2, index, ge_addback_gate = 25, ge_sum_gate = 25;
  char *name, c;
  long ts;

  chan = ptr->chan;
  // Check for invalid channel numbers, prossibly due to data corruption
  if( chan < 0 || chan > odb_daqsize ){
    fprintf(stderr,"Invalid channel number in fill_singles_histos(), %d\n",chan);
    return(-1);
  }
  sys = ptr->subsys;
  // Check this is a valid susbsytem type
  if( sys <0 || sys > MAX_SUBSYS ){
    if( mult_hist[sys] != NULL && sys>=0 && sys<MAX_SUBSYS){ mult_hist[sys]->Fill(mult_hist[sys], ptr->fold, 1);  }
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

         // The PU class is assigned in the presort, use ptr->psd for pileup class
         ge_pu_class->Fill(ge_pu_class, ptr->psd, 1);  // pileup class value
         ge_sum_class[ptr->psd]->Fill(ge_sum_class[ptr->psd], (int)ptr->ecal, 1);  // energy spectrum of pileup class value
         ge_e_vs_k_class[ptr->psd]->Fill(ge_e_vs_k_class[ptr->psd], (int)ptr->ecal, (int)ptr->integ, 1);  // energy spectrum of pileup class value

         if(ptr->psd == 1){  // single hit
           ge_1hit[pos]->Fill(ge_1hit[pos], (int)ptr->ecal, 1);
           ge_xtal_1hit->Fill(ge_xtal_1hit, pos, (int)ptr->ecal, 1);
          }
         if(ptr->psd > 9 && ptr->psd < 13){ // 3-hit pileup
           ge_3hit[pos]->Fill(ge_3hit[pos], (int)ptr->ecal, 1);
           ge_xtal_3hit->Fill(ge_xtal_3hit, pos, (int)ptr->ecal, 1);
          }

         if(ptr->psd == 3 || ptr->psd == 5 || ptr->psd == 7){ // Select first Hit of two pileup events
           ge_2hit[pos]->Fill(ge_2hit[pos], (int)ptr->ecal, 1); // 2-hit pileup
           ge_xtal_2hit->Fill(ge_xtal_2hit, pos, (int)ptr->ecal, 1);
         // The following used for mapping the k2 dependant correction
           ge_e_vs_k_2hit_first[pos]->Fill(ge_e_vs_k_2hit_first[pos], (int)ptr->ecal, (int)ptr->integ, 1);  // energy1 vs k1 spectrum of 1st Hit of 2hit pileup events
         }
         if(ptr->psd == 4 || ptr->psd == 6 || ptr->psd == 8){ // Select second Hit of two pileup events
           ge_2hit[pos]->Fill(ge_2hit[pos], (int)ptr->ecal, 1); // 2-hit pileup
           ge_xtal_2hit->Fill(ge_xtal_2hit, pos, (int)ptr->ecal, 1);
           // The following is not used for mapping the corrections but is a useful diagnostic
           ge_e_vs_k_2hit_second[pos]->Fill(ge_e_vs_k_2hit_second[pos], (int)ptr->ecal, (int)ptr->integ, 1);  // energy2 vs k2 spectrum of 2nd Hit of 2hit pileup events

           // The following 1408keV matrix used for mapping the E1 offset correction
           if(ptr->e4cal > 1380 && ptr->e4cal < 1420){ // Require 152Eu 1408keV as E1 for mapping the E1 offset
             ge_PU2_e2_v_k_gated1408[pos]->Fill(ge_PU2_e2_v_k_gated1408[pos], (int)ptr->ecal, (int)ptr->integ, 1);  // e2 vs k2 for fixed e1
           }

           // The following x-rays matrix used for mapping the k2 dependant correction for E2
           // E2 vs k2 gated on fixed x-ray energies
           // The E2 energy has 1272keV subtracted from it to put it around 136keV to allow a smaller matrix side and easier processing in the app
           if(ptr->e4cal > 5 && ptr->e4cal < 50){ // Require 152Eu x rays as E1 for mapping the k2 dependance
             ge_PU2_e2_v_k_gatedxrays[pos]->Fill(ge_PU2_e2_v_k_gatedxrays[pos], (int)(ptr->ecal - 1272), (int)ptr->integ, 1);  // e2 vs k2 for fixed e1
           }
         }

      //   if(ptr->psd==7){ ge_pu_dt12->Fill(ge_pu_dt12, (int)(ptr->ts_int+DT_SPEC_LENGTH/2), 1); }
      //   if(ptr->psd==8){ ge_pu_dt13->Fill(ge_pu_dt13, (int)(ptr->ts_int+DT_SPEC_LENGTH/2), 1); }


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
       }else {
         fprintf(stderr,"bad ge crystal[%d] for chan %d\n", pos, ptr->chan);
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
   case SUBSYS_LABR_T:

   // Save LBL channel number into ptr->q2 or q3 or q4
   // Save LBL energy ecal into TAC ptr-ecal2 or ecal3 or ecal4
   if(ptr->e4cal>0){ break; } // more than two LaBr3 in coincidence with this TAC event so reject
   if( ptr->q2 >=          0 && ptr->q3 >=           0 &&
       ptr->q2 < MAX_DAQSIZE && ptr->q3  < MAX_DAQSIZE ){
      c1 = crystal_table[ptr->q2]-1; // c1 runs from 0 to 7
      c2 = crystal_table[ptr->q3]-1; // c2 runs from 0 to 7
   } else { c1 = c2 = -1; }
   if(c1>=0 && c1<N_LABR && c2>=0 && c2<N_LABR){
     index = tac_labr_hist_index[c1][c2];
     if(index>=0 && index<(int)((N_LABR)*(N_LABR-1)/2)+1){
       tac_labr_hist[index]->Fill(tac_labr_hist[index], (int)(ptr->ecal), 1);
     }
     // Compton Walk matrix for calibrations
     // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
     if(crystal_table[ptr->chan] == 1){ // Use TAC01
       if(c1 == 0 && c2>0 && ptr->e2cal>1252 && ptr->e2cal<1412){ // LBL01 gated on 1332keV
         tac_labr_CompWalk[c2]->Fill(tac_labr_CompWalk[c2], (int)(ptr->ecal/4), (int)ptr->e3cal, 1); // TAC01 and other LBL energy
       }
       if(c1 == 0 && c2==1 && ptr->e3cal>1252 && ptr->e3cal<1412){ // LBL02 gated on 1332keV
         tac_labr_CompWalk[c1]->Fill(tac_labr_CompWalk[c1], (int)(ptr->ecal/4), (int)ptr->e2cal, 1); // TAC01 and other LBL energy
       }
     }
   } break;
   case SUBSYS_DESCANT: break;
   case SUBSYS_DESWALL: // DESCANT Wall
    pos  = crystal_table[chan];
    if(pos>0 && pos<=N_DES_WALL){
       desw_psd[pos]       -> Fill(desw_psd[pos],   (int)ptr->psd,       1);
       if(ptr->energy4>0){
          desw_tof[pos]       -> Fill(desw_tof[pos],   (int)ptr->energy4,       1);
          desw_tof_corr[pos]  -> Fill(desw_tof_corr[pos],   (int)ptr->e4cal,       1);
          if(ptr->psd>10 && ptr->psd<710){
             desw_tof_psd[pos]  -> Fill(desw_tof_psd[pos],   (int)ptr->e4cal,       1);
          }
       }
    }
    desw_sum_e->Fill(desw_sum_e, (int)ptr->ecal, 1);
    if(ptr->psd>0){ desw_sum_psd->Fill(desw_sum_psd, (int)ptr->psd, 1); }
    if(ptr->e4cal>0){ desw_sum_tof->Fill(desw_sum_tof, (int)ptr->e4cal, 1); } // e4cal = corrected time-of-flight
    pos  = crystal_table[ptr->chan];
    if( pos < 1 || pos > 60 ){
      fprintf(stderr,"bad descant wall detector[%d] for chan %d\n", pos, ptr->chan);
    } else {
        desw_e_xtal->Fill(desw_e_xtal, pos, (int)ptr->ecal, 1);
        desw_tof_xtal->Fill(desw_tof_xtal, pos, (int)ptr->e4cal, 1);
        desw_psd_e->Fill(desw_psd_e, (int)ptr->psd, (int)ptr->ecal, 1);
        desw_psd_tof->Fill(desw_psd_tof, (int)ptr->psd, (int)ptr->e4cal, 1);
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
  case SUBSYS_ZDS_A: gc_hist->Fill(gc_hist, 2, 1); break;
  case SUBSYS_ZDS_B: gc_hist->Fill(gc_hist, 1, 1); break;
   case SUBSYS_RCMP:
       rcmp_sum->Fill(rcmp_sum, (int)ptr->ecal, 1);
       if(ptr->esum>0){
         rcmp_fb_sum->Fill(rcmp_fb_sum, (int)ptr->esum, 1);
       }
       pos  = crystal_table[ptr->chan];
       elem = (int)(elem + (int)(polarity_table[ptr->chan]*N_RCMP_STRIPS)); // polarity_table value is 0 or 1
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
   int g_aries_upper_gate=25;
   int g_aries_tac_gate=60;
   int g_rcmp_upper_gate=75;
   int gg_gate=25, c1, c2, angle_idx;
   switch(alt->subsys){
   case SUBSYS_HPGE_A:
      if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) ){
         if( ptr->esum >= 0 &&  alt->esum >= 0 ){ // addback energies
            gg_ab->Fill(gg_ab, (int)ptr->esum, (int)alt->esum, 1);
         }
         c1 = crystal_table[ptr->chan];
         c2 = crystal_table[alt->chan];
         if( c1 >= 0 && c1 < 64 && c2 >= 0 && c2 < 64 ){
            gg_hit->Fill(gg_hit, c1, c2, 1); // 2d crystal hitpattern
            if( c2 == grif_opposite[c1] ){
               // 180 degree coinc matrix for summing corrections
               gg_opp->Fill(gg_opp, (int)ptr->ecal, (int)alt->ecal, 1);
            }
            // Ge-Ge angular correlations
            // Fill the appropriate angular bin spectrum
            // c1 and c2 run from 0 to 63 for ge_angles_145mm.
            angle_idx = ge_angles_110mm[c1][c2];
            gg_angcor_110[angle_idx]->Fill(gg_angcor_110[angle_idx], (int)ptr->ecal, (int)alt->ecal, 1);
            angle_idx = ge_angles_145mm[c1][c2];
            gg_angcor_145[angle_idx]->Fill(gg_angcor_145[angle_idx], (int)ptr->ecal, (int)alt->ecal, 1);
         }
      }
      break;
   case SUBSYS_SCEPTAR:
      if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) ){
         ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
         ge_sum_b_sep->Fill(ge_sum_b_sep, (int)ptr->ecal, 1); // Sceptar-gated Ge sum energy spectrum
      }
      break;
      case SUBSYS_ARIES_A:
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
      }
      break;
   case SUBSYS_ZDS_A:
      if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ZDS_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ZDS_A]) ){
         ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1);         // beta-gated Ge sum energy spectrum
         ge_sum_b_zds->Fill(ge_sum_b_zds, (int)ptr->ecal, 1); // Zds-gated Ge sum energy spectrum
      }
      break;
   case SUBSYS_DESWALL: // ge-DSW
      ge_dsw->Fill(ge_dsw, (int)ptr->ecal, (int)alt->e4cal, 1); // e4cal = DSW corrected time-of-flight
      break;
   default: break;
   }
   return(0);
}

int fill_labr_coinc_histos(Grif_event *ptr, Grif_event *alt, int abs_dt)
{
   int lbl_tac_gate=15, tac_offset[8] = {-7300,-5585,-6804,0,-6488,-5682,-5416,0};
   int g_aries_upper_gate=25, c1, c2, corrected_tac_value;
   switch(alt->subsys){
   case SUBSYS_ARIES_A:
      if( (abs_dt >= time_diff_gate_min[SUBSYS_LABR_L][SUBSYS_ARIES_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_LABR_L][SUBSYS_ARIES_A]) ){
         c1 = crystal_table[ptr->chan];
         c2 = crystal_table[alt->chan];
         if(c1 >= 1 && c1 <=8 && c2 >= 1 && c2 <=76 ){
            lba_hit->Fill(lba_hit, c1, c2, 1);
         }
      } break;
   case SUBSYS_LABR_T:
      c1=crystal_table[alt->chan]-1;  // assign c1 as TAC number
      if(c1 >= 0 && c1 < 8 ){ // 8 LBL
         dt_tacs_hist[c1]->Fill(dt_tacs_hist[c1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
         if(((abs_dt >= time_diff_gate_min[SUBSYS_LABR_L][SUBSYS_LABR_T]) && (abs_dt <= time_diff_gate_max[SUBSYS_LABR_L][SUBSYS_LABR_T])) && c1 == 8){ // labr-tac with the ARIES TAC
            c2 = crystal_table[ptr->chan] - 1; // assign c2 as LBL number
            if(c2 >= 0 && c2 < 8 ){ // 8 LBL detectors
               corrected_tac_value = (int)alt->ecal + tac_offset[c2];
               tac_aries_lbl[c2]->Fill(tac_aries_lbl[c2], corrected_tac_value, 1); // tac spectrum per LBL
               tac_aries_lbl_sum->Fill(tac_aries_lbl_sum, corrected_tac_value, 1); // sum tac spectrum including all LBL
               if(ptr->ecal >1225 && ptr->ecal <1315){ // gate on LaBr3 energy 1275keV
                  aries_tac->Fill(aries_tac, (int)corrected_tac_value, 1); // tac spectrum gated on 1275keV
                  aries_tac_artEn->Fill(aries_tac_artEn, alt->e4cal, 1); // ARIES energy spectrum in coincidence with TAC
                  if(alt->e4cal >24 && alt->e4cal <36){ // gate on ARIES energy
                     aries_tac_Egate->Fill(aries_tac_Egate, corrected_tac_value, 1); // tac spectrum gated on 1275keV
                  }
               }
            }
         }
      } break;
   }
   return(0);
}

int reorder_rcmp_US_strips[32] = { // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25966
    1, 3, 5, 7, 9,11,13,15,31,29,27,25,23,21,19,17,
   16,18,20,22,24,26,28,30,14,12,10, 8, 6, 4, 2, 0
};
int reorder_rcmp_DS_strips[32] = { // Per GRIFFIN elog, https://grsilog.triumf.ca/GRIFFIN/25968
     0, 2, 4, 6, 8,10,12,14,30,28,26,24,22,20,18,16,
    17,19,21,23,25,27,29,31,15,13,11, 9, 7, 5, 3, 1
};

int frag_hist[MAX_COINC_EVENTS];
int fill_coinc_histos(int win_idx, int frag_idx)
{
   int global_window_size = 100; // size in grif-replay should be double this
   Grif_event *alt, *ptr = &grif_event[win_idx], *tmp;
   int dt, abs_dt,  pos, c1, c2, index, ptr_swap;
   int gg_gate=25, g_aries_upper_gate=25;
   TH2I *hist_ee; TH1I *hist_dt;

   // histogram of coincwin-size
   dt = (frag_idx - win_idx + 2*MAX_COINC_EVENTS) %  MAX_COINC_EVENTS; ++frag_hist[dt];

   while( win_idx != frag_idx ){ // check all conicidences in window
      if( ++win_idx == MAX_COINC_EVENTS ){ win_idx = 0; } // wrap
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
         }} break;
      case SUBSYS_RCMP:
         if(alt->subsys == SUBSYS_RCMP){
           if( (abs_dt >= time_diff_gate_min[SUBSYS_RCMP][SUBSYS_RCMP]) && (abs_dt <= time_diff_gate_max[SUBSYS_RCMP][SUBSYS_RCMP]) ){
            if((pos = crystal_table[ptr->chan]) == crystal_table[alt->chan] &&
                  polarity_table[ptr->chan] != polarity_table[alt->chan] ){ // front and back of same DSSD
            c1 = element_table[ptr->chan];
            c2 = element_table[alt->chan];
            rcmp_hit[(pos-1)]->Fill(rcmp_hit[(pos-1)], c1, c2, 1);
            rcmp_fb[(pos-1)]->Fill(rcmp_fb[(pos-1)], (int)ptr->ecal, (int)alt->ecal, 1);
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
         if( alt->subsys == SUBSYS_LABR_T && crystal_table[alt->chan] == 8 ){ // ARIES TAC
            dt_hist[14]->Fill(dt_hist[14], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
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
            dt_hist[20]->Fill(dt_hist[20], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
            art_dsw->Fill(art_dsw, (int)ptr->ecal, (int)alt->e4cal, 1);
            desw_sum_e_b->Fill(desw_sum_e_b, (int)alt->ecal, 1);
            desw_sum_tof_b->Fill(desw_sum_tof_b, (int)alt->e4cal, 1); // e4cal = corrected time-of-flight
         } break;
      case SUBSYS_ZDS_A: // grif16 zds
         if(alt->subsys == SUBSYS_ZDS_B ){ // ZDS GRIF-CAEN coincidence
            gc_hist->Fill(gc_hist, 3, 1);
            gc_hist->Fill(gc_hist, 5, 1);
            dt_hist[22]->Fill(dt_hist[22], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
            dt_hist[23]->Fill(dt_hist[23], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
         } break;
      case SUBSYS_ZDS_B: // CAEN zds
         if( alt->subsys == SUBSYS_DESWALL ){ // ZDS-DSW
            dt_hist[21]->Fill(dt_hist[21], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
            dt_hist[25]->Fill(dt_hist[25], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
            desw_sum_e_b->Fill(desw_sum_e_b, (int)alt->ecal, 1);
            desw_sum_tof_b->Fill(desw_sum_tof_b, (int)alt->e4cal, 1); // e4cal = corrected time-of-flight
         } break;
      case SUBSYS_DESWALL:
         if( alt->subsys == SUBSYS_DESWALL ){
            dt_hist[18]->Fill(dt_hist[18], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
            dt_hist[24]->Fill(dt_hist[24], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
            c1 = crystal_table[ptr->chan]-1; c2 = crystal_table[alt->chan]-1;
            if( c1 >= 0 && c1 < 60 && c2 >= 0 && c2 < 60 ){
               dsw_hit->Fill(dsw_hit, c1, c2, 1);
               dsw_dsw->Fill(dsw_dsw, (int)ptr->e4cal, (int)alt->e4cal, 1);
               // Fold 2 sum spectra
               desw_sum_e_nn->Fill(desw_sum_e_nn, (int)ptr->ecal, 1);
               desw_sum_tof_nn->Fill(desw_sum_tof_nn, (int)ptr->e4cal, 1);// e4cal = corrected time-of-flight
               // DSW-DSW angular correlations
               // Fill the appropriate angular bin spectrum with the corrected time-of-flight value
               index = DSW_DSW_angles[c1][c2];
               dsw_angcor[index]->Fill(dsw_angcor
                 [index],(int)ptr->e4cal,(int)alt->e4cal, 1);
               // Fold 2, angle greater than 60 degrees, sum spectra
               // index 13 = 58.555, index 14 = 61.535
               if( index > 13 ){                                         // e4cal = corrected time-of-flight
                  desw_sum_e_nn_a->Fill(desw_sum_e_nn_a, (int)ptr->ecal, 1);
                  desw_sum_tof_nn_a->Fill(desw_sum_tof_nn_a, (int)ptr->e4cal, 1);
               }
            }
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
#define ODBHANDLE_UNK  23
static char odb_handle[MAX_ODB_SUBSYS][8] = {
   "GRG", "GRS", "SEP", "PAC",  //  0- 3
   "LBS", "LBT", "LBL", "DSC",  //  4- 7
   "ART", "ZDS", "RCS", "XXX",  //  8-11
   "DSW", "DSG",    "",    "",
   "",    "",    "",    "",
   "",    "",    "",    "UNK"
};

static char   path[256];
static char dirname[64],value[32],type[32];
extern char midas_runtitle[SYS_PATH_LENGTH];

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
         i = (str-ptr) > SYS_PATH_LENGTH-1 ? SYS_PATH_LENGTH-1 : (str-ptr);
         memcpy( midas_runtitle, ptr, i ); midas_runtitle[i] = 0;
         ptr += i+1;
         if( (str = strchr(ptr,d)) == NULL ){ break; }
      } else if( strncmp(ptr,"</keyarray>",10) == 0 ){ active = 0; arrayptr = (void *)('\0');
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
      case ODBHANDLE_PAC: case ODBHANDLE_ZDS: case ODBHANDLE_DSW:
         crystal_table[i] = pos;
         if(        crystal == 'A' ){ element_table[i] = 1;
         } else if( crystal == 'B' ){ element_table[i] = 2;
         } else if( crystal == 'C' ){ element_table[i] = 3;
         } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3, ZDS, ART
         } else {
            fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
         } break;
      case ODBHANDLE_RCS:
         crystal_table[i] = pos;
         element_table[i] = element;
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
      case ODBHANDLE_LBT: subsys_table[i] = SUBSYS_LABR_T;    break;
      case ODBHANDLE_LBL: subsys_table[i] = SUBSYS_LABR_L;    break;
      case ODBHANDLE_DSC: subsys_table[i] = SUBSYS_DESCANT;   break;
      case ODBHANDLE_RCS: subsys_table[i] = SUBSYS_RCMP;      break;
      case ODBHANDLE_DSW: subsys_table[i] = SUBSYS_DESWALL;  break;
      case ODBHANDLE_DSG: subsys_table[i] = SUBSYS_DSG;  break;
      case ODBHANDLE_GRG: subsys_table[i] = (output_type == 1) ? SUBSYS_HPGE_A :SUBSYS_HPGE_B; break;
      case ODBHANDLE_ZDS: subsys_table[i] = (output_type == 1) ? SUBSYS_ZDS_A  :SUBSYS_ZDS_B;  break;
      case ODBHANDLE_ART: subsys_table[i] = (polarity_table[i] == 1) ? SUBSYS_ARIES_A:SUBSYS_ARIES_B;break;
      case ODBHANDLE_XXX: subsys_table[i] = SUBSYS_IGNORE;    break;
      case ODBHANDLE_UNK: subsys_table[i] = SUBSYS_UNKNOWN;   break;
      }
   }
   memset(subsys_initialized, 0, sizeof(int)*MAX_SUBSYS );

   return(0);
}
