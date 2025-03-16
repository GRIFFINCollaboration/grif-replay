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
int     subsys_dtype_mat[MAX_SUBSYS][MAX_SUBSYS]; // map using names and dtypes
int         subsys_dtype[MAX_SUBSYS];             // subsys->dtype, usage: subsys_dtype[dtype] = subsys_handle
int         dtype_subsys[MAX_SUBSYS];             // dtype->subsys
int     psc_dtype_subsys[MAX_SUBSYS];             // PSC psc_dtype->subsys, usage: psc_dtype_subsys[subsys_handle] = dtype from PSC table
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
float  pileupk1[MAX_DAQSIZE][7];
float  pileupk2[MAX_DAQSIZE][7];
float  pileupE1[MAX_DAQSIZE][7];
static short *chan_address = addr_table;
extern Grif_event grif_event[MAX_COINC_EVENTS];

// Default sort function declarations
extern int init_time_diff_gates(Config *cfg);
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
   init_singles_histos(cfg);
   init_coinc_histos(cfg);

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

  // Initialize all time differences between subsystems to be the default 250ns
  for(i=0; i<MAX_SUBSYS; i++){
    for(j=0; j<MAX_SUBSYS; j++){
      time_diff_gate_min[i][j] = 0;  // default is 0 nanoseconds
      time_diff_gate_max[i][j] = 25; // default is 250 nanoseconds
    }
  }

  // Search the globals for time difference settings and overwrite their values
  for(i=0; i<cfg->nglobal; i++){
    global = cfg->globals[i];
    sprintf(tmp,"%s",global->name);
    if(strncmp(tmp,"time_diff_",10) == 0){
      //fprintf(stdout,"Process %s\n",global->name);
      // This global is a time difference value
      // Identify the subsystem types and then save the value in the correct place
      for(j=0; j<MAX_SUBSYS; j++){
        if(strlen(subsys_handle[j])<2){ continue; }
        if(strstr(tmp,subsys_handle[j]) > 0){
          // Identiy the second subsystem type
          //  fprintf(stdout,"Found %s in %s\n",subsys_handle[j],global->name);
          for(k=0; k<MAX_SUBSYS; k++){
            if(strlen(subsys_handle[k])<2){ continue; }
            if(strstr(tmp,subsys_handle[k]) > 0 && (strstr(tmp,subsys_handle[k]) != strstr(tmp,subsys_handle[j]))){
              // save the value in the correct place
              fprintf(stdout,"time_diff_%s_%s set to %d,%d\n",subsys_handle[j],subsys_handle[k],global->min,global->max);
              time_diff_gate_min[j][k] = global->min;
              time_diff_gate_max[j][k] = global->max;
              time_diff_gate_min[k][j] = global->min;
              time_diff_gate_max[k][j] = global->max;
            }
          }
        }
      }
    }
  }

  return(0);
}

int apply_gains(Grif_event *ptr)
{
  //int tac_ts_offset[8] = {50,58,405,73,73,404,110,154};
  int tac_ts_offset[12] = {60,60,60,60,60,60,60,60,60,60,60,60}; // From Dec 2024
  //int tac_ts_offset[8] = {134,48,74,59,48,400,395,0}; // Rashmi S1723
  int caen_ts_offset = -60; // this value (-60) aligns the timestamps of HPGe with ZDS(CAEN)
   float energy,psd;
   int chan;
   if( (chan=ptr->chan) >= odb_daqsize ){
      fprintf(stderr,"unpack_event: ignored event in chan:%d [0x%04x]\n",
            	                                  chan, ptr->address );
      return(-1);
   }
   if( chan<0 ){
      fprintf(stderr,"unpack_event: ignored event with negative chan:%d\n",
            	                                  chan );
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

   // Assign the Sub System index based on dtype
   // The dtype to subsys mapping was determined from the PSC table in the function gen_derived_odb_tables()
    if( ptr->dtype >= 0 && ptr->dtype < MAX_SUBSYS ){
        ptr->subsys = dtype_subsys[ptr->dtype];
      if( debug ){ printf("--SET EVT[%4ld]=%d\n", ptr - grif_event, ptr->subsys ); }
      if( ptr->subsys != subsys_dtype[dtype_table[ptr->chan]] ){
         // Hack for non-DESCANT things in CAEN electronics
         // All CAEN channels are set to dtype 6 by default.
         // Reassign subsys value for these three channels that are not DSW
         if(ptr->subsys == SUBSYS_DES_WALL){ // These are all CAEN electronics channels
          ptr->subsys = subsys_dtype[dtype_table[ptr->chan]];
          //ptr->ts -= caen_ts_offset; // Subtract from CAEN timestamps to align coincidences

        }else{
                // non-CAEN electronics channel so there is an error here
        	     printf("--ERROR... Channel %d: There is a mismatch in the subsys type assigned [%d, %s] in comparison to the assignment in the PSC table [%d, %s]\n",ptr->chan,ptr->subsys,subsys_name[ptr->subsys],subsys_dtype[dtype_table[ptr->chan]],subsys_name[subsys_dtype[dtype_table[ptr->chan]]]);
        }
      }
   } else { ptr->subsys = MAX_SUBSYS-1; }

   // fprintf(stdout,"apply_gains %s chan%d: %d/%d=%d",subsys_handle[ptr->subsys],chan,ptr->q,ptr->integ,ptr->energy);
   // fprintf(stdout,", [%f,%f,%f] -> %d\n",quads[chan],gains[chan],offsets[chan],ptr->ecal);


   // HPGe pileup development
   if( ptr->subsys == SUBSYS_HPGE){
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
   //if( ptr->subsys == SUBSYS_DESCANT || ptr->subsys == SUBSYS_DES_WALL){
   if( ptr->subsys == SUBSYS_DES_WALL){
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
  int bgo_window = 20, addback_window = 20;
  int rcmp_fb_window = 10;
  int lbl_tac_window = 25;
  int art_tac_window = 25;
  int desw_beta_window = 80;
  float desw_median_distance = 1681.8328; // descant wall median source-to-detector distance in mm
  int i, j, dt, dt13, tof;
  float q1,q2,q12,k1,k2,k12,e1,e2,e12,m,c;
  int chan,found,pos;
  float energy,correction;
  float correction12, correction23;

      // Protect yourself
  if( ptr->chan<0 || ptr->chan >= odb_daqsize ){
     fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
     return(-1);
  }
  //printf("\n");

  //if( ptr->dtype ==  6 ){
  // printf("Dsc\n");;
  // }

  //if( ptr->dtype !=  0 ){ return(0); } // not Ge event
  //if( ptr->dtype == 15 ){ return(0); } // scalar
  i = frag_idx; ptr->fold = 1;
  while( i != end_idx ){ // need at least two events in window
    if( ++i >=  MAX_COINC_EVENTS ){ i=0; } alt = &grif_event[i]; // WRAP

    // Protect yourself
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
      case SUBSYS_HPGE:

      // HPGe pile-up corrections
      // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
      // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
      if(alt->subsys == SUBSYS_HPGE && alt->chan == ptr->chan){
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
            // Two-Hit pileup case ...
            //
            //      |    |   K1   | /|      K12      |\
            //      |    *________|/_|_______________| \
            //      |   /                             \ \ |    K2   |   .
            //      |  /                               \ \|_________|   .
            //      | /                                 \            \  .
            //    __*/                                   \_____       \___
            //
            //

            // The (ptr) fragement is the first Hit of a two Hit pile-up event.
            // It is identified as having (ptr->pileup==1 && ptr->nhit==2)

            // Assign the pileup class numbers to the two hits and calculate the time difference between them
            ptr->psd = 3; // Pileup class, 1st of 2Hits
            alt->psd = 4; // Pileup class, 2nd of 2Hits
            ptr->ts_int = alt->ts_int = dt;


            // Apply the k1 dependant correction to the energy of the first hit
            chan  = ptr->chan; // chan for ptr and alt are the same
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

          }else{
            // 2Hit error events - q12 is zero
            // 2 Hit, Type B

            // Assign the pileup class numbers to the two hits and calculate the time difference between them
            ptr->psd = 7; // Pileup class, 1st of 2Hits where Hits treated separately with no correction
            alt->psd = 8; // Pileup class, 2nd of 2Hits where Hits treated separately with no correction
            ptr->ts_int = alt->ts_int = dt;

            // Apply the k1 dependant correction to the energy of the first hit
            pos  = crystal_table[ptr->chan];
            k1 = ptr->integ;
            ptr->ecal=ptr->esum = ptr->ecal*( pileupk1[chan][0]+(k1*pileupk1[chan][1])+(k1*k1*pileupk1[chan][2])+(k1*k1*k1*pileupk1[chan][3])
            +(k1*k1*k1*k1*pileupk1[chan][4])+(k1*k1*k1*k1*k1*pileupk1[chan][5])+(k1*k1*k1*k1*k1*k1*pileupk1[chan][6]));
            alt->e4cal=ptr->ecal; // Remember the ecal of the first Hit in this second Hit

            // Apply the k2 dependant correction to the energy of the second hit
            k2 = alt->integ;
            correction = ptr->ecal*( pileupE1[chan][0]+(k2*pileupE1[chan][1])+(k2*k2*pileupE1[chan][2])+(k2*k2*k2*pileupE1[chan][3])
            +(k2*k2*k2*k2*pileupE1[chan][4])+(k2*k2*k2*k2*k2*pileupE1[chan][5])+(k2*k2*k2*k2*k2*k2*pileupE1[chan][6]));
            alt->ecal=alt->esum = (alt->ecal*( pileupk2[chan][0]+(k2*pileupk2[chan][1])+(k2*k2*pileupk2[chan][2])+(k2*k2*k2*pileupk2[chan][3])
            +(k2*k2*k2*k2*pileupk2[chan][4])+(k2*k2*k2*k2*k2*pileupk2[chan][5])+(k2*k2*k2*k2*k2*k2*pileupk2[chan][6])))+correction;

          }
        }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==1 && alt->nhit==1)){
          // 2 Hit, Type C
          // 2Hit pileup where 2nd Hit integration region starts after 1st Hit integration ends
          // k12<0, q12 is zero
          // -> Treat as separate hits
          // Correct second Hit for effect of first based on time between hits

          // Assign the pileup class numbers to the two hits and calculate the time difference between them
          ptr->psd = 5; // Pileup class, 1st of 2Hits where Hits treated separately with no correction
          alt->psd = 6; // Pileup class, 2nd of 2Hits where Hits treated separately with no correction
          ptr->ts_int = alt->ts_int = dt;

          // Apply the k1 dependant correction to the energy of the first hit
          pos  = crystal_table[ptr->chan];
          k1 = ptr->integ;
          ptr->ecal=ptr->esum = ptr->ecal*( pileupk1[chan][0]+(k1*pileupk1[chan][1])+(k1*k1*pileupk1[chan][2])+(k1*k1*k1*pileupk1[chan][3])
          +(k1*k1*k1*k1*pileupk1[chan][4])+(k1*k1*k1*k1*k1*pileupk1[chan][5])+(k1*k1*k1*k1*k1*k1*pileupk1[chan][6]));
          alt->e4cal=ptr->ecal; // Remember the ecal of the first Hit in this second Hit

          // Apply the k2 dependant correction to the energy of the second hit
          k2 = alt->integ;
          correction = ptr->ecal*( pileupE1[chan][0]+(k2*pileupE1[chan][1])+(k2*k2*pileupE1[chan][2])+(k2*k2*k2*pileupE1[chan][3])
          +(k2*k2*k2*k2*pileupE1[chan][4])+(k2*k2*k2*k2*k2*pileupE1[chan][5])+(k2*k2*k2*k2*k2*k2*pileupE1[chan][6]));
          alt->ecal=alt->esum = (alt->ecal*( pileupk2[chan][0]+(k2*pileupk2[chan][1])+(k2*k2*pileupk2[chan][2])+(k2*k2*k2*pileupk2[chan][3])
          +(k2*k2*k2*k2*pileupk2[chan][4])+(k2*k2*k2*k2*k2*pileupk2[chan][5])+(k2*k2*k2*k2*k2*k2*pileupk2[chan][6])))+correction;

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

    //    fprintf(stdout,"\n");
      }
    }
  }

      // BGO suppression of HPGe
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
      // In the TAC event we save the tile chan as ab_alt_chan, and the tile energy as e4cal.
      // So later in the main coincidence loop we only need to examine LBL and TAC.
      if(dt < art_tac_window && alt->subsys == SUBSYS_ARIES && crystal_table[ptr->chan] == 8){
      ptr->ab_alt_chan = alt->chan; ptr->e4cal = alt->ecal;
      }
      // For TAC01-07 we have a LBL-LBL coincidence
      // Here save the LBL Id number and the LBL energy in the TAC event
      // Save LBL channel number into ptr->q2 or q3 or q4
      // Save LBL energy ecal into TAC ptr-ecal2 or ecal3 or ecal4
      if(dt < lbl_tac_window && alt->subsys == SUBSYS_LABR_L && crystal_table[ptr->chan] < 8){
        if(ptr->e2cal<1){
          ptr->q2 = alt->chan; ptr->e2cal = alt->ecal;
        }else if(ptr->e3cal<1){
          ptr->q3 = alt->chan; ptr->e3cal = alt->ecal;
        }else{
          ptr->q4 = alt->chan; ptr->e4cal = alt->ecal; // If this is set then we have LBL multiplicity >2 for this TAC
        }
      }
      break;
      case SUBSYS_ZDS:
      if(output_table[ptr->chan]==0){ // CAEN Zds
        if(alt->subsys == SUBSYS_DES_WALL){
          if(dt < desw_beta_window){
          // Calculate time-of-flight and correct it for this DESCANT detector distance
          tof = (spread(abs(ptr->cfd - alt->cfd))*2.0) + 100; //if( tof < 0 ){ tof = -1*tof; }
          //  fprintf(stdout,"tof: %d - %d = %f\n",ptr->cfd, alt->cfd, tof);
          alt->energy4 = (int)(tof); // Time of flight
          alt->e4cal = (int)(spread(tof) * DSW_tof_corr_factor[crystal_table[alt->chan]-1]); // Corrected Time of Flight
        }
      }

      }
      break;
      /*
      //case SUBSYS_DESCANT:
      case SUBSYS_DES_WALL:
      // DESCANT detectors
      // use e4cal for Time-Of-Flight which is derived from the time difference between a beta hit and the DESCANT hit - equivalent to neutron energy
      if(dt < desw_beta_window){
	//  if( ((alt->subsys == SUBSYS_ARIES && polarity_table[alt->chan] == 0) || (alt->subsys == SUBSYS_ZDS && output_table[alt->chan]==0)) && alt->ecal > 5){
  // Use ARIES Fast output (polarity_table[alt->chan] == 0) in CAEN electronics
  // Use ZDS B output (output_table[alt->chan] == 0) in CAEN electronics
        //if( ((alt->subsys == SUBSYS_ARIES) && (polarity_table[alt->chan] == 0))
        if(alt->subsys == SUBSYS_ZDS && output_table[alt->chan]==0){
        // Calculate time-of-flight and correct it for this DESCANT detector distance
        tof = (ptr->cfd - alt->cfd) + 1000; //if( tof < 0 ){ tof = -1*tof; }
	  //  fprintf(stdout,"tof: %d - %d = %f\n",ptr->cfd, alt->cfd, tof);
        ptr->energy4 = (int)(tof); // Time of flight
        ptr->e4cal = (int)(tof * DSW_tof_corr_factor[crystal_table[ptr->chan]-1]); // Corrected Time of Flight
      }
    }

      break;
*/

      default: // Unrecognized or unprocessed subsys type
      break;
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
   int test_matrix_data[32][32] = {
     {2008,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,0,0,3, 3,0,3,3, 3,3,0,3, 3,0,0,3},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,0,0,3, 3,0,3,3, 3,3,0,3, 3,0,0,3},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,0,0,3, 3,0,3,3, 3,3,0,3, 3,0,0,3},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3},

     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,0,0,3, 3,0,0,3, 3,0,0,3, 3,0,0,3},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  3,0,0,3, 3,0,0,3, 3,0,0,3, 3,0,0,3},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,2008,  3,0,0,3, 3,0,0,3, 3,0,0,3, 3,0,0,3},

     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {0,0,8,0, 0,0,0,8, 0,0,0,0, 8,8,0,0,  0,0,8,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {0,8,2008,8, 0,0,8,0, 8,0,0,8, 0,0,0,0,  0,8,8,8, 0,0,0,0, 0,0,0,0, 0,0,0,0},

     {0,0,8,0, 0,0,8,8, 8,0,0,8, 8,0,0,0,  0,0,8,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {0,0,8,0, 0,0,8,0, 0,0,0,0, 0,8,0,0,  0,0,8,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {0,0,8,8, 0,0,0,8, 8,0,0,8, 8,8,0,0,  0,0,8,8, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},


     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,0,0,3, 3,0,3,3, 3,3,0,3, 3,0,0,3},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,0,0,3, 3,0,3,3, 3,3,0,3, 3,0,0,3},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,0,0,3, 3,0,3,3, 3,3,0,3, 3,0,0,3},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3},

     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,0,0,3, 3,0,0,3, 3,0,0,3, 3,0,0,3},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,0,0,3, 3,0,0,3, 3,0,0,3, 3,0,0,3},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  3,0,0,3, 3,0,0,3, 3,0,0,3, 3,0,0,3},

     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {10,10,2008,10, 10,10,10,2008, 10,10,10,10, 2008,2008,0,0,  0,0,8,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {10,2008,2008,2008, 10,10,2008,10, 2008,10,10,2008, 0,0,0,0,  0,8,8,8, 0,0,0,0, 0,0,0,0, 0,0,0,0},

     {10,10,2008,10, 10,10,2008,2008, 2008,10,10,2008, 2008,10,10,10,  0,0,8,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {10,10,2008,10, 10,10,2008,10, 10,10,10,10, 10,2008,10,10,  0,0,8,0, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {10,10,2008,2008, 10,10,10,2008, 2008,10,10,2008, 2008,2008,10,10,  0,0,8,8, 0,0,0,0, 0,0,0,0, 0,0,0,0},
     {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,  0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0}
   };

   open_folder(cfg, "Hits_and_Sums");
   open_folder(cfg, "Hits");
   for(i=0; i<N_HITPAT; i++){ // Create Hitpattern spectra
      sprintf(title,  "Hitpattern_%s",    hit_names[i] );
      hit_hist[i] = H1_BOOK(cfg, hit_handles[i], title, MAX_DAQSIZE, 0, MAX_DAQSIZE);
   }
   ts_hist = H1_BOOK(cfg, "ts", "Timestamp", 163840, 0, 163840);
   gc_hist = H1_BOOK(cfg, "gc", "ZDS GRIF-CAEN", 16, 0, 16);
   test_histogram = H2_BOOK(cfg, "Test", "Test matrix", 32, 0, 32, 32, 0, 32); // non symmeterized
  // test_histogram = H2_BOOK(cfg, "Test", "Test matrix", 32, 0, 32, 0, 0, 32); // symmeterized
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
      if( strcmp(subsys_handle[j],"XXX") == 0 ){     // suppress XXX channels
         ph_hist[i]->suppress = e_hist[i]->suppress =
            /*  cfd_hist[i]->suppress = wave_hist[i]->suppress = */ 1;
      }
   }

   // Fill the test matrix here after initialization
   for(i=0; i<32; i++){
     for(j=0; j<32; j++){
       for(k=0; k<test_matrix_data[i][j]; k++){
         test_histogram->Fill(test_histogram, j, i, 1);
       }
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
   if( ptr->dtype == 6 ){
      ++sys;
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
   if( ptr->subsys == SUBSYS_DES_WALL){
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
   }

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

int init_singles_histos(Config *cfg)
{
  char title[STRING_LEN], handle[STRING_LEN];
  int i,j;
  open_folder(cfg, "Hits_and_Sums");
  open_folder(cfg, "Sums");
  sprintf(title,  "Addback_Sum_Energy"); sprintf(handle, "Addback_Sum_E");
  ge_sum_ab = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_Energy"); sprintf(handle, "Ge_Sum_E");
  ge_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_En_betaTagged"); sprintf(handle, "Ge_Sum_E_B");
  ge_sum_b = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_En_SceptarTagged"); sprintf(handle, "Ge_Sum_E_B_SEP");
  ge_sum_b_sep = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_En_ZdsTagged"); sprintf(handle, "Ge_Sum_E_B_ZDS");
  ge_sum_b_zds = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Ge_Sum_En_AriesTagged"); sprintf(handle, "Ge_Sum_E_B_ART");
  ge_sum_b_art = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
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
  sprintf(title,  "RCMP_Sum_Energy"); sprintf(handle, "RCMP_Sum_E");
  rcmp_sum = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);

  sprintf(title,  "DES_Wall_Sum_Energy"); sprintf(handle, "DESW_Sum_E");
  desw_sum_e = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_TOF"); sprintf(handle, "DESW_Sum_TOF");
  desw_sum_tof = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_PSD"); sprintf(handle, "DESW_Sum_PSD");
  desw_sum_psd = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_En_betaTagged"); sprintf(handle, "DESW_Sum_E_B");
  desw_sum_e_b = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_TOF_betaTagged"); sprintf(handle, "DESW_Sum_TOF_B");
  desw_sum_tof_b = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_En_fold2"); sprintf(handle, "DESW_Sum_E_nn");
  desw_sum_e_nn = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_TOF_fold2"); sprintf(handle, "DESW_Sum_TOF_nn");
  desw_sum_tof_nn = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_En_fold2_ang60"); sprintf(handle, "DESW_Sum_E_nn_a");
  desw_sum_e_nn_a = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "DES_Wall_Sum_TOF_fold2_ang60"); sprintf(handle, "DESW_Sum_TOF_nn_a");
  desw_sum_tof_nn_a = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);

  sprintf(title,  "Upstream_Ge_Sum_Energy"); sprintf(handle, "US_Ge_Sum_E");
  ge_sum_us = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Downstream_Ge_Sum_Energy"); sprintf(handle, "DS_Ge_Sum_E");
  ge_sum_ds = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Upstream_AB_Sum_Energy"); sprintf(handle, "US_AB_Sum_E");
  ge_sum_ab_us = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  sprintf(title,  "Downstream_AB_Sum_Energy"); sprintf(handle, "DS_AB_Sum_E");
  ge_sum_ab_ds = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  for(i=0; i<N_CLOVER; i++){
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
  sprintf(title, "TAC_LBL_ART_vs_LBL_Num"); sprintf(handle, "TAC_ART_LBL_LBL_Xtal");
  labr_tac_xtal     = H2_BOOK(cfg, handle, title, 16, 0, 16,
			                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "TAC_LBL_ART_vs_ART_Num"); sprintf(handle, "TAC_ART_LBL_ART_Xtal");
  art_tac_xtal     = H2_BOOK(cfg, handle, title, 80, 0, 80,
			                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "DES_Wall_En_DetNum"); sprintf(handle, "DSW_En_Xtal");
  desw_e_xtal     = H2_BOOK(cfg, handle, title, 64, 0, 64,
                                            E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "DES_Wall_TOF_DetNum"); sprintf(handle, "DSW_TOF_Xtal");
  desw_tof_xtal     = H2_BOOK(cfg, handle, title, 64, 0, 64,
                                            E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "DES_Wall_PSD_En"); sprintf(handle, "DES_Wall_PSD_En");
  desw_psd_e     = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                            E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "DES_Wall_PSD_TOF"); sprintf(handle, "DES_Wall_PSD_TOF");
  desw_psd_tof     = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                            E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
    for(i=0; i<N_RCMP_POS; i++){ // Create RCMP DSSD strip spectra
      rcmp_strips[i] = H2_BOOK(cfg, rcmp_strips_handles[i], rcmp_strips_handles[i], 2*N_RCMP_STRIPS, 0, 2*N_RCMP_STRIPS,
                                            E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
    }
  close_folder(cfg);
  open_folder(cfg, "Pile-up");
  sprintf(title,  "Pile-up_type"); sprintf(handle, "PU_type");
  ge_pu_type = H1_BOOK(cfg, handle, title, 64, 0, 64);
  sprintf(title,  "nhits_type"); sprintf(handle, "nhits_type");
  ge_nhits_type = H1_BOOK(cfg, handle, title, 64, 0, 64);
  sprintf(title,  "Pile-up_class"); sprintf(handle, "PU_class");
  ge_pu_class = H1_BOOK(cfg, handle, title, 64, 0, 64);
  for(i=0; i<NUM_PILEUP_CLASSES; i++){
    sprintf( title, "%s", ge_pu_class_titles[i]);
    sprintf(handle, "%s", ge_pu_class_handles[i]);
    ge_sum_class[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
  }
  for(i=0; i<NUM_PILEUP_CLASSES; i++){
    sprintf( title, "E_vs_k_%s", ge_pu_class_titles[i]);
    sprintf(handle, "E_vs_k_%s", ge_pu_class_handles[i]);
    ge_e_vs_k_class[i] = H2_BOOK(cfg, handle, title, 2048, 0, 2048,
                                                      512, 0,  512);
  }
  sprintf(title, "GeEnergy_CrystalNum_single_hit"); sprintf(handle, "GeEnergy_Xtal_single_hit");
  ge_xtal_1hit     = H2_BOOK(cfg, handle, title, 64, 0, 64,
			                                    E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GeEnergy_CrystalNum_2_hit"); sprintf(handle, "GeEnergy_Xtal_2_hit");
  ge_xtal_2hit     = H2_BOOK(cfg, handle, title, 64, 0, 64,
			                                    E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GeEnergy_CrystalNum_3_hit"); sprintf(handle, "GeEnergy_Xtal_3_hit");
  ge_xtal_3hit     = H2_BOOK(cfg, handle, title, 64, 0, 64,
			                                    E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);

for(i=0; i<N_HPGE; i++){
    sprintf( title, "Ge%02d_single_hit", i); sprintf(handle, "Ge%02d_single_hit", i);
    ge_1hit[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
    sprintf( title, "Ge%02d_2hit_pileup", i); sprintf(handle, "Ge%02d_2hit_pileup", i);
    ge_2hit[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
    sprintf( title, "Ge%02d_3hit_pileup", i); sprintf(handle, "Ge%02d_3hit_pileup", i);
    ge_3hit[i] = H1_BOOK(cfg, handle, title, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
}

    sprintf(title,  "Pile-up_3Hits_detlaT_1_2"); sprintf(handle, "PU_dt12");
    ge_pu_dt12 = H1_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
    sprintf(title,  "Pile-up_3Hits_detlaT_1_3"); sprintf(handle, "PU_dt13");
    ge_pu_dt13 = H1_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);

  close_folder(cfg);
  open_folder(cfg, "Pile-up Corrections");
                  // 2d energy1 vs k1 matrix for pileup corrections
                  // Use for mapping k1 dependence of E1
  for(i=0; i<N_HPGE; i++){
    sprintf( title, "Ge%02d_E_vs_k_1st_of_2hit", i);
    sprintf(handle, "Ge%02d_E_vs_k_1st_of_2hit", i);
    ge_e_vs_k_2hit_first[i] = H2_BOOK(cfg, handle, title, 2048, 0, 2048,
                                                    512, 0,  512);
  }

                  // 2d energy2 vs k2 matrix for pileup corrections
                  // Not actually used for mapping corrections, but useful diagnostic
  for(i=0; i<N_HPGE; i++){
    sprintf( title, "Ge%02d_E_vs_k_2nd_of_2hit", i);
    sprintf(handle, "Ge%02d_E_vs_k_2nd_of_2hit", i);
    ge_e_vs_k_2hit_second[i] = H2_BOOK(cfg, handle, title, 2048, 0, 2048,
                                                    512, 0,  512);
  }

                // 2d energy2 vs k2 matrix for pileup corrections (will be gated on E1 = x rays)
                // Use for mapping k2 dependence of E2
        for(i=0; i<N_HPGE; i++){
                    sprintf( title, "Ge%02d_PU2_E2_vs_k2_E1gated_on_Xrays", i);
                    sprintf( handle, "Ge%02d_PU2_E2_vs_k2_E1gated_on_Xrays", i);
                    ge_PU2_e2_v_k_gatedxrays[i] = H2_BOOK(cfg, handle, title, 2048, 0, 2048,
                                                                            512, 0, 512);
          }

          // 2d energy1 vs energy2 matrix for pileup corrections
          // Use for mapping E1 dependence of E2
          for(i=0; i<N_HPGE; i++){
            sprintf( title, "Ge%02d_PU2_E2_vs_k2_E1gated_on_1408keV", i);
            sprintf( handle, "Ge%02d_PU2_E2_vs_k2_E1gated_on_1408keV", i);
            ge_PU2_e2_v_k_gated1408[i] = H2_BOOK(cfg, handle, title, 512, 0, 512,
                                                                     512, 0, 512);
            }

/*
            // To be removed...
            for(i=0; i<N_HPGE; i++){
              for(j=0; j<N_K; j++){
                sprintf( title, "Ge%02d_PU2_E2_vs_E1_k%d", i,((j*20)+10));
                sprintf(handle, "Ge%02d_PU2_E2_vs_E1_k%d", i,((j*20)+10));
                ge_PU2_e2_v_e1_k[i][j] = H2_BOOK(cfg, handle, title, 2048, 0, 2048,
                                                                      512, 0,  512);
                }
              }
              */

  close_folder(cfg);
  close_folder(cfg);
  return(0);
}

int init_coinc_histos(Config *cfg)
{
  char title[STRING_LEN], handle[STRING_LEN];
  char tmp[STRING_LEN];
  int i,j,k;

  // Set the ybins as SYMMETERIZE to create a symmeterized 2d histogram.

  open_folder(cfg, "Hits_and_Sums");
  open_folder(cfg, "Delta_t");
  for(i=0; i<N_DT; i++){ // Create delta-t spectra
    dt_hist[i] = H1_BOOK(cfg, dt_handles[i], dt_handles[i], DT_SPEC_LENGTH, 0, DT_SPEC_LENGTH);
  }
  for(i=0; i<N_TACS; i++){ // Create delta-t spectra separately for each TAC
  sprintf(title, "dt_labr_tac%d",i+1);
    dt_tacs_hist[i] = H1_BOOK(cfg, title, title, DT_SPEC_LENGTH, 0, DT_SPEC_LENGTH);
  }
  close_folder(cfg);
  close_folder(cfg);
  open_folder(cfg, "Coinc");
  open_folder(cfg, "Coinc");
  sprintf(title, "Addback_GG"); sprintf(handle, "Addback_GG");
  gg_ab     = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                      SYMMETERIZE, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GG"); sprintf(handle, "GG");
  gg        = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                      SYMMETERIZE, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GePaces"); sprintf(handle, "GePaces");
  ge_paces  = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "GeLabr"); sprintf(handle, "GeLabr");
   ge_labr   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                       E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "GeZds"); sprintf(handle, "GeZds");
    ge_zds   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
 		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
     sprintf(title, "GeAries"); sprintf(handle, "GeAries");
      ge_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
   		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "LaBrAries"); sprintf(handle, "LaBrAries");
    labr_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
   		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "PacesAries"); sprintf(handle, "PacesAries");
    paces_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
     		                                         E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "AriesAries"); sprintf(handle, "AriesAries");
    art_art   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
     		                                        SYMMETERIZE, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "LaBrLabr"); sprintf(handle, "LaBrLabr");
   labr_labr   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
 		                                         SYMMETERIZE, 0, E_2D_SPEC_LENGTH);
sprintf(title, "LaBrZds"); sprintf(handle, "LaBrZds");
labr_zds   = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                          E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
     sprintf(title, "DSWDSW"); sprintf(handle, "DSWDSW");
     dsw_dsw        = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
   		                                      E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
    sprintf(title, "GeDSW"); sprintf(handle, "GeDSW");
    ge_dsw        = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                           E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "ARTDSW"); sprintf(handle, "ARTDSW");
   art_dsw        = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                          E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH);
   sprintf(title, "GeRCMP"); sprintf(handle, "GeRCMP");
   ge_rcmp  = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
 		                                      E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
   sprintf(title, "LaBrRCMP"); sprintf(handle, "LaBrRCMP");
   labr_rcmp  = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                            E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
   sprintf(title, "GGoppo"); sprintf(handle, "GGoppo");
   gg_opp    = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
		                                       SYMMETERIZE, 0, E_2D_SPEC_LENGTH);
  sprintf(title, "Addback_GGoppo"); sprintf(handle, "AB_GGoppo");
  gg_ab_opp    = H2_BOOK(cfg, handle, title, E_2D_SPEC_LENGTH, 0, E_2D_SPEC_LENGTH,
                                          SYMMETERIZE, 0, E_2D_SPEC_LENGTH);
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
   sprintf(title, "DSWDSWHit"); sprintf(handle, "DSWDSWHit");
   dsw_hit    = H2_BOOK(cfg, handle, title, 64, 0, 64,
                                           64, 0, 64);
     for(i=0; i<N_RCMP_POS; i++){ // Create RCMP DSSD hitpatterns
       rcmp_hit[i] = H2_BOOK(cfg, rcmp_hit_handles[i], rcmp_hit_handles[i], N_RCMP_STRIPS, 0, N_RCMP_STRIPS,
                                                                            N_RCMP_STRIPS, 0, N_RCMP_STRIPS);
     }
for(i=0; i<N_RCMP_POS; i++){ // Create RCMP DSSD front energy vs back energy heatmaps
rcmp_fb[i] = H2_BOOK(cfg, rcmp_fb_handles[i], rcmp_fb_handles[i], E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH,
                                    E_2D_RCMP_SPEC_LENGTH, 0, E_2D_RCMP_SPEC_LENGTH);
}
   close_folder(cfg);
   close_folder(cfg);
   open_folder(cfg, "Ang_Corr");
   open_folder(cfg, "GG_Ang_Corr");
   for(i=0; i<N_GE_ANG_CORR; i++){ // Create Ge-Ge angular correlation spectra for HPGe at 110mm
     sprintf(tmp,"Ge-Ge_110mm_angular_bin%d",i);
     gg_ang_corr_110_hist[i] = H2_BOOK(cfg, tmp, tmp, GE_ANG_CORR_SPEC_LENGTH, 0, GE_ANG_CORR_SPEC_LENGTH,
                                                              SYMMETERIZE, 0, GE_ANG_CORR_SPEC_LENGTH);
   }
   for(i=0; i<N_GE_ANG_CORR; i++){ // Create Ge-Ge angular correlation spectra for HPGe at 145mm
     sprintf(tmp,"Ge-Ge_145mm_angular_bin%d",i);
     gg_ang_corr_145_hist[i] = H2_BOOK(cfg, tmp, tmp, GE_ANG_CORR_SPEC_LENGTH, 0, GE_ANG_CORR_SPEC_LENGTH,
                                                              SYMMETERIZE, 0, GE_ANG_CORR_SPEC_LENGTH);
   }
   close_folder(cfg);
   open_folder(cfg, "Ge_ART_Ang_Corr");
   for(i=0; i<N_GRG_ART_ANG_CORR; i++){ // Create Ge-ART angular correlation spectra
     sprintf(tmp,"Ge-ART_angular_bin%d",i);
     grg_art_ang_corr_hist[i] = H2_BOOK(cfg, tmp, tmp, GE_ANG_CORR_SPEC_LENGTH, 0, GE_ANG_CORR_SPEC_LENGTH,
                                                       GE_ANG_CORR_SPEC_LENGTH, 0, GE_ANG_CORR_SPEC_LENGTH);
   }
   close_folder(cfg);
   open_folder(cfg, "DSW_Ang_Corr");
   for(i=0; i<N_DSW_DSW_ANG_CORR; i++){ // Create DSW-DSW angular correlation spectra
     sprintf(tmp,"DSW-DSW_angular_bin%d",i);
     dsw_dsw_ang_corr_hist[i] = H2_BOOK(cfg, tmp, tmp, DSW_ANG_CORR_SPEC_LENGTH, 0, DSW_ANG_CORR_SPEC_LENGTH,
                                                                    SYMMETERIZE, 0, DSW_ANG_CORR_SPEC_LENGTH);
   }
   close_folder(cfg);
   close_folder(cfg);
   open_folder(cfg, "Fast-Timing");
   open_folder(cfg, "LBL_TACs");
   for(i=0; i<N_LABR; i++){ // intialize the index for the TAC coincidence pair spectra
     for(j=0; j<N_LABR; j++){
       tac_labr_hist_index[i][j] = -1;
     }
   }
   k=-1;
   for(i=1; i<=N_LABR; i++){ // Create TAC coincidence pair spectra
     for(j=(i+1); j<=N_LABR; j++){
       k++;
     sprintf(tmp,"TAC_%d_%d",i,j);
     tac_labr_hist[k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
     tac_labr_hist_index[i-1][j-1] = k;
     }
   }
 sprintf(tmp,"TAC_%d_%d",2,1); // Add additional histogram (2_1) needed for Compton Walk corrections
 tac_labr_hist[++k] = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
 tac_labr_hist_index[1][0] = k;

   close_folder(cfg);
      open_folder(cfg, "LBL_Walk");

      for(i=1; i<=N_LABR; i++){ // Create TAC Compton Walk matrices for calibrations
        sprintf(tmp,"TAC_%d_CompWalk",i);
        tac_labr_CompWalk[i-1] = H2_BOOK(cfg, tmp, tmp, 2048, 0, 2048,
      			                                            2048, 0, 2048);
      }
      close_folder(cfg);
   open_folder(cfg, "ART_TACs");
     sprintf(tmp,"TAC_ART_LBL_LBLSUM");
   tac_aries_lbl_sum = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
     sprintf(tmp,"TAC_ART_LBL_ARTSUM");
   tac_aries_art_sum = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
 sprintf(tmp,"TAC-ARIES-LaBr3-1275keV");
 aries_tac = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
sprintf(tmp,"TAC-ARTE-LaBr3-1275keV");
 aries_tac_Egate = H1_BOOK(cfg, tmp, tmp, E_TAC_SPEC_LENGTH, 0, E_TAC_SPEC_LENGTH);
sprintf(tmp,"TAC-ARIES-Energy");
 aries_tac_artEn = H1_BOOK(cfg, tmp, tmp, E_SPEC_LENGTH, 0, E_SPEC_LENGTH);
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

  // Check that this is the correctly assigned subsystem type, based on the datatype in the PSC table
  //if( sys != dtype_subsys[dtype_table[ptr->chan]] ){
  if( sys != subsys_dtype[dtype_table[ptr->chan]] ){
    fprintf(stderr,"Subsystem assigned [%s,%d] does not match expectation from PSC datatype [%d] which gives sys handle %d for channel %d\n",subsys_handle[sys],sys,dtype_table[ptr->chan],subsys_dtype[dtype_table[ptr->chan]],ptr->chan);
    return(-1);
  }

  // Get the position for this fragment
  pos  = crystal_table[ptr->chan];

  // Here we should not use dtype because the dtype can change between experiments.
  // Here we use the Subsytem name, which is obtained from dtype.
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
         if(pos<32){
           ge_sum_us->Fill(ge_sum_us, (int)ptr->ecal, 1);
         }else{
           ge_sum_ds->Fill(ge_sum_ds, (int)ptr->ecal, 1);
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

         if(ptr->psd==7){ ge_pu_dt12->Fill(ge_pu_dt12, (int)(ptr->ts_int+DT_SPEC_LENGTH/2), 1); }
         if(ptr->psd==8){ ge_pu_dt13->Fill(ge_pu_dt13, (int)(ptr->ts_int+DT_SPEC_LENGTH/2), 1); }

         // Ge-Addback
         clover = (int)(pos/16)+1;
         if( clover >= 0 && clover < N_CLOVER && ptr->esum >= 0 ){   // ge addback
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
      }
      break;
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
      }
      break;
   case SUBSYS_PACES: // PACES
    paces_sum->Fill(paces_sum, (int)ptr->ecal, 1);
    pos  = crystal_table[ptr->chan];
    if( pos < 1 || pos > 5 ){
      fprintf(stderr,"bad PACES crystal[%d] for chan %d\n", pos, ptr->chan);
    } else {
      paces_xtal->Fill(paces_xtal, pos, (int)ptr->ecal, 1);
    }
    break;
   case SUBSYS_LABR_L: // LaBr3 (LBL)
    labr_sum->Fill(labr_sum, (int)ptr->ecal, 1);
      pos  = crystal_table[ptr->chan];
      if( pos < 1 || pos > 8 ){
         fprintf(stderr,"bad LaBr3 crystal[%d] for chan %d\n", pos, ptr->chan);
      } else {
         labr_xtal->Fill(labr_xtal, pos, (int)ptr->ecal, 1);
      }
      break;
   case SUBSYS_SCEPTAR:
   break;
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
   }

   break;
   case SUBSYS_DESCANT:
   break;
   case SUBSYS_DES_WALL: // DESCANT Wall
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
   case SUBSYS_ARIES: // ARIES
   if(polarity_table[ptr->chan] == 1){ // Only use ARIES Standard Output
     aries_sum->Fill(aries_sum, (int)ptr->ecal, 1);
     pos  = crystal_table[ptr->chan];
     if( pos < 1 || pos > 76 ){
       fprintf(stderr,"bad aries tile[%d] for chan %d\n", pos, ptr->chan);
     } else {
       aries_xtal->Fill(aries_xtal, pos, (int)ptr->ecal, 1);
     }
   }
   break;
   case SUBSYS_ZDS:
    if(output_table[ptr->chan]==0){ // CAEN ZDS
     gc_hist->Fill(gc_hist, 2, 1);
    }
    if(output_table[ptr->chan]==1){ // GRIF16 ZDS
     gc_hist->Fill(gc_hist, 1, 1);
    }
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
  int lbl_tac_gate=15;
  int pos, c1, c2, index;
  int global_window_size = 100; // size in grif-replay should be double this
  int corrected_tac_value;
  int tac_offset[8] = {-7300,-5585,-6804,0,-6488,-5682,-5416,0};

        // Protect yourself
    if( ptr->chan<0 || ptr->chan >= odb_daqsize ){
       fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
       return(-1);
    }

  // histogram of coincwin-size
  dt = (frag_idx - win_idx + 2*MAX_COINC_EVENTS) %  MAX_COINC_EVENTS;
  ++frag_hist[dt];
  if( win_idx == frag_idx ){ return(0); } // window size = 1 - no coinc
  if( (i=win_idx+1) == MAX_COINC_EVENTS ){ i = 0; } // wrap
  // check all conicidences in window
  while( 1 ){
    alt = &grif_event[i];

          // Protect yourself
      if( alt->chan<0 || alt->chan >= odb_daqsize ){
         fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
         return(-1);
      }

    abs_dt = dt = ptr->ts - alt->ts; if( dt < 0 ){ abs_dt = -1*dt; }
    if( abs_dt > global_window_size ){ break; }

    switch(ptr->subsys){
      case SUBSYS_HPGE: // Ge matrices

      // Only use GRGa
      if(output_table[ptr->chan] == 1){

        switch(alt->subsys){
          case SUBSYS_HPGE:
          // Only use GRGa
          if(output_table[alt->chan] == 1){

            dt_hist[0]->Fill(dt_hist[0], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // ge-ge
            // Ge-Ge matrices
            if( abs_dt < gg_gate ){ // ge-ge (and addback)
              gg->Fill(gg, (int)ptr->ecal, (int)alt->ecal, 1);
              if( ptr->esum >= 0 &&  alt->esum >= 0 ){
                gg_ab->Fill(gg_ab, (int)ptr->esum, (int)alt->esum, 1);
              }

              c1 = crystal_table[ptr->chan];
              if( c1 >= 0 && c1 <= 63 ){
                c2 = crystal_table[alt->chan];
                if( c2 >= 0 && c2 <= 63 ){
                  gg_hit->Fill(gg_hit, c1, c2, 1);

                  // Ge-Ge with 180 degrees between Ge1 and Ge2 used for summing corrections
                  if( c2 == grif_opposite[c1] ){
                    gg_opp->Fill(gg_opp, (int)ptr->ecal, (int)alt->ecal, 1);
                    gg_ab_opp->Fill(gg_ab_opp, (int)ptr->esum, (int)alt->esum, 1);
                  }

                  // Ge-Ge angular correlations
                  // Fill the appropriate angular bin spectrum
                  // c1 and c2 run from 0 to 63 for ge_angles_145mm.
                  index = ge_angles_110mm[c1][c2];
                  gg_ang_corr_110_hist[index]->Fill(gg_ang_corr_110_hist[index], (int)ptr->ecal, (int)alt->ecal, 1);
                  index = ge_angles_145mm[c1][c2];
                  gg_ang_corr_145_hist[index]->Fill(gg_ang_corr_145_hist[index], (int)ptr->ecal, (int)alt->ecal, 1);
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
            ge_sum_b_sep->Fill(ge_sum_b_sep, (int)ptr->ecal, 1); // Sceptar-gated Ge sum energy spectrum
          }
          break;
          case SUBSYS_ARIES: // ge-aries
          if(polarity_table[alt->chan] == 1){ // Only use ARIES Standard Output
            dt_hist[10]->Fill(dt_hist[10], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
            if( abs_dt < g_aries_upper_gate ){
              ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
              ge_sum_b_art->Fill(ge_sum_b_art, (int)ptr->ecal, 1); // Aries-gated Ge sum energy spectrum
              ge_art->Fill(ge_art, (int)ptr->ecal, (int)alt->esum, 1);
              c1 = crystal_table[ptr->chan];
              if(c1 >= 1 && c1 <=64 ){
                c2 = crystal_table[alt->chan];
                if(c2 >= 1 && c2 <=76 ){
                  gea_hit->Fill(gea_hit, c1, c2, 1);
                  c1--; c2--;

                  // Ge-ARIES angular correlations
                  // Fill the appropriate angular bin spectrum
                  index = GRG_ART_angles_110mm[c1][c2];
                  grg_art_ang_corr_hist[index]->Fill(grg_art_ang_corr_hist[index], (int)ptr->ecal, (int)alt->ecal, 1);
                }
              }
            }
          }
          break;
          case SUBSYS_ZDS: // ge-zds
          if(output_table[alt->chan]==1){ // GRIF16 ZDS
            dt_hist[3]->Fill(dt_hist[3], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
            if( abs_dt < gg_gate ){
              ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
              ge_sum_b_zds->Fill(ge_sum_b_zds, (int)ptr->ecal, 1); // Zds-gated Ge sum energy spectrum
                ge_zds->Fill(ge_zds, (int)ptr->ecal, (int)alt->ecal, 1);
            }
          }
          break;
          case SUBSYS_RCMP: // ge-rcmp
          dt_hist[6]->Fill(dt_hist[6], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          if( abs_dt < g_rcmp_upper_gate && alt->esum>0){
            ge_rcmp->Fill(ge_rcmp, (int)ptr->ecal, (int)alt->esum, 1);
          }
          break;
          case SUBSYS_DES_WALL: // ge-DSW
          dt_hist[19]->Fill(dt_hist[19], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          ge_dsw->Fill(ge_dsw, (int)ptr->e4cal, (int)alt->ecal, 1);
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
            ge_sum_b_sep->Fill(ge_sum_b_sep, (int)alt->ecal, 1); // Sceptar-gated Ge sum energy spectrum
          }
        }
        break;
        default: break;
      }
      break;

      case SUBSYS_PACES: // paces matrices
      switch(alt->subsys){
        case SUBSYS_ZDS: // paces-zds
        if(output_table[alt->chan]==1){ // GRIF16 ZDS
          dt_hist[7]->Fill(dt_hist[7], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-zds
        }
        break;
        case SUBSYS_HPGE: // paces-ge
        // Only use GRGa
        if(output_table[alt->chan] == 1){
          dt_hist[4]->Fill(dt_hist[4], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-ge
        }
        break;
        case SUBSYS_LABR_L: // paces-LBL
        dt_hist[8]->Fill(dt_hist[8], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-labr
        break;
        case SUBSYS_ARIES: // paces-aries
        if(polarity_table[alt->chan] == 1){ // Only use ARIES Standard Output
          dt_hist[13]->Fill(dt_hist[13], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // pac-aries
          if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0){
            paces_art->Fill(paces_art, (int)ptr->ecal, (int)alt->ecal, 1);
          }
        }
        break;
        default: break; // unprocessed coincidence combinations
      } // end of inner switch(ALT)

      case SUBSYS_LABR_L: // Labr matrices
      switch(alt->subsys){
        case SUBSYS_HPGE: // labr-ge
        // Only use GRGa
        if(output_table[alt->chan] == 1){
          dt_hist[5]->Fill(dt_hist[5], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // labr-ge
          if( abs_dt < gg_gate ){
            ge_labr->Fill(ge_labr, (int)alt->ecal, (int)ptr->ecal, 1);
          }
        }
        break;
        case SUBSYS_LABR_L: // labr-labr
        dt_hist[8]->Fill(dt_hist[8], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // labr-labr
        labr_labr->Fill(labr_labr, (int)ptr->ecal, (int)alt->ecal, 1);
        break;
        case SUBSYS_ZDS: // labr-zds
        if(output_table[alt->chan]==1){ // GRIF16 ZDS
          dt_hist[17]->Fill(dt_hist[17], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
          labr_zds->Fill(labr_zds, (int)ptr->ecal, (int)alt->ecal, 1);
        }
        break;
        case SUBSYS_ARIES: // labr-aries
        if(polarity_table[alt->chan] == 1){ // Only use ARIES Standard Output
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
        }
        break;
        case SUBSYS_LABR_T: // labr-tac
        dt_hist[16]->Fill(dt_hist[16], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        c1=crystal_table[alt->chan];  // assign c1 as TAC number
        if(c1 >= 1 && c1 <=8 ){ // 8 LBL

          dt_tacs_hist[c1-1]->Fill(dt_tacs_hist[c1-1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);


          if(abs_dt < lbl_tac_gate && c1 == 8){ // labr-tac with the ARIES TAC
            c2 = crystal_table[ptr->chan] - 1; // assign c2 as LBL number

            if(c2 >= 0 && c2 <8 ){ // 8 LBL detectors
              corrected_tac_value = (int)alt->ecal + tac_offset[c2];
              tac_aries_lbl_hist[c2]->Fill(tac_aries_lbl_hist[c2], corrected_tac_value, 1); // tac spectrum per LBL
              tac_aries_lbl_sum->Fill(tac_aries_lbl_sum, corrected_tac_value, 1); // sum tac spectrum including all LBL
              if(ptr->ecal >1225 && ptr->ecal <1315){ // gate on LaBr3 energy 1275keV
                aries_tac->Fill(aries_tac, (int)corrected_tac_value, 1); // tac spectrum gated on 1275keV
                aries_tac_artEn->Fill(aries_tac_artEn, alt->e4cal, 1); // ARIES energy spectrum in coincidence with TAC
                if(alt->e4cal >24 && alt->e4cal <36){ // gate on ARIES energy
                  aries_tac_Egate->Fill(aries_tac_Egate, corrected_tac_value, 1); // tac spectrum gated on 1275keV
                }
              }
            }
          } // end of ARIES TAC

        } // end of all LaBr3 TAC

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
  if(polarity_table[ptr->chan] == 1){ // Only use ARIES Standard Output
    switch(alt->subsys){
      case SUBSYS_HPGE: // aries-ge
      // Only use GRGa
      if(output_table[alt->chan] == 1){

        dt_hist[10]->Fill(dt_hist[10], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        if( ( abs_dt < g_aries_upper_gate) && ptr->ecal>0 && alt->ecal>0){
          ge_sum_b->Fill(ge_sum_b, (int)alt->ecal, 1); // beta-gated Ge sum energy spectrum
          ge_sum_b_art->Fill(ge_sum_b_art, (int)alt->ecal, 1); // Aries-gated Ge sum energy spectrum
          ge_art->Fill(ge_art, (int)alt->ecal, (int)ptr->ecal, 1);
          c1 = crystal_table[ptr->chan];
          if(c1 >= 1 && c1 <=76 ){
            c2 = crystal_table[alt->chan];
            if(c2 >= 1 && c2 <=64 ){
              gea_hit->Fill(gea_hit, c2, c1, 1);
              c1--; c2--;

              // Ge-ARIES angular correlations
              // Fill the appropriate angular bin spectrum
              index = GRG_ART_angles_110mm[c2][c1];
              grg_art_ang_corr_hist[index]->Fill(grg_art_ang_corr_hist[index], (int)alt->ecal, (int)ptr->ecal, 1);
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
        tac_aries_art_sum->Fill(tac_aries_art_sum, (int)alt->ecal, 1); // sum tac spectrum including all art
        c2 = crystal_table[ptr->chan]-1;
        if(c2>=0 && c2<N_ARIES){
        tac_aries_art_hist[c2]->Fill(tac_aries_art_hist[c2], (int)alt->ecal, 1); // tac spectrum per ART tiles
      }
      }
      break;
      default: break;
    } // end of inner switch(ALT)
  }else{ // end of if output_table for ARIES SO
    // ARIES Fast Output in CAEN
    if(alt->subsys == SUBSYS_DES_WALL){ // aries-DSW
      dt_hist[20]->Fill(dt_hist[20], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      art_dsw->Fill(art_dsw, (int)ptr->ecal, (int)alt->e4cal, 1);

      desw_sum_e_b->Fill(desw_sum_e_b, (int)alt->ecal, 1);
      desw_sum_tof_b->Fill(desw_sum_tof_b, (int)alt->e4cal, 1); // e4cal = corrected time-of-flight
    }
  }
  break; // end of inner switch(ALT) for ptr=aries


  case SUBSYS_ZDS: // zds matrices
  if(output_table[ptr->chan]==1){ // GRIF16 ZDS
    switch(alt->subsys){
      case SUBSYS_ZDS: // ZDS GRIF-CAEN
      if(output_table[alt->chan]==0){ // ZDS GRIF-CAEN coincidence
        gc_hist->Fill(gc_hist, 3, 1);
        gc_hist->Fill(gc_hist, 5, 1);
        dt_hist[22]->Fill(dt_hist[22], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        dt_hist[23]->Fill(dt_hist[23], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
      }
      break;
      case SUBSYS_LABR_T: // zds-tac
      dt_hist[15]->Fill(dt_hist[15], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      break;
      case SUBSYS_HPGE: // zds-HPGe
      // Only use GRGa
      if(output_table[alt->chan] == 1){
        dt_hist[3]->Fill(dt_hist[3], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        if( abs_dt < gg_gate ){
          ge_sum_b->Fill(ge_sum_b, (int)alt->ecal, 1); // beta-gated Ge sum energy spectrum
          ge_sum_b_zds->Fill(ge_sum_b_zds, (int)alt->ecal, 1); // Zds-gated Ge sum energy spectrum
          ge_zds->Fill(ge_zds, (int)alt->ecal, (int)ptr->ecal, 1);
        }
      }
      break;
      case SUBSYS_LABR_L: // zds-labr
      dt_hist[17]->Fill(dt_hist[17], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      break;
      default: break;
    } // end of inner switch(ALT)
  }else if(output_table[ptr->chan]==0){ // CAEN ZDS
    if(alt->subsys == SUBSYS_ZDS && output_table[alt->chan]==1){ // ZDS CAEN-GRIF coincidence
      gc_hist->Fill(gc_hist, 4, 1);
      gc_hist->Fill(gc_hist, 5, 1);
      dt_hist[22]->Fill(dt_hist[22], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      dt_hist[23]->Fill(dt_hist[23], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);
    }

    if(alt->subsys == SUBSYS_DES_WALL){ // ZDS-DSW
      dt_hist[21]->Fill(dt_hist[21], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      dt_hist[25]->Fill(dt_hist[25], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);

      desw_sum_e_b->Fill(desw_sum_e_b, (int)alt->ecal, 1);
      desw_sum_tof_b->Fill(desw_sum_tof_b, (int)alt->e4cal, 1); // e4cal = corrected time-of-flight
    }

  } // end of CAEN ZDS
  break; // end of ptr ZDS

  case SUBSYS_LABR_T: // tac matrices
  switch(alt->subsys){
    case SUBSYS_ZDS: // tac-zds
    if(output_table[alt->chan]==1){ // GRIF16 ZDS
      dt_hist[15]->Fill(dt_hist[15], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    }
    break;
    case SUBSYS_ARIES: // tac-aries
    if(polarity_table[alt->chan] == 1){ // Only use ARIES Standard Output
      if(crystal_table[ptr->chan] == 8){ // ARIES TAC
        dt_hist[14]->Fill(dt_hist[14], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
        tac_aries_art_hist[crystal_table[alt->chan]-1]->Fill(tac_aries_art_hist[crystal_table[alt->chan]-1], (int)ptr->ecal, 1); // tac spectrum per ART tiles
        tac_aries_art_sum->Fill(tac_aries_art_sum, (int)ptr->ecal, 1); // sum tac spectrum including all art
      }
    }
    break;
    case SUBSYS_LABR_L: // tac-labr
    dt_hist[16]->Fill(dt_hist[16], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    c1 = crystal_table[ptr->chan]-1;
    if(c1>=0 && c1<N_LABR){
      dt_tacs_hist[c1]->Fill(dt_tacs_hist[c1], (int)(abs_dt+DT_SPEC_LENGTH/2), 1); // LBL TACs


      if(abs_dt < lbl_tac_gate && crystal_table[ptr->chan] == 8){ // ARIES TAC

        corrected_tac_value = (int)ptr->ecal + tac_offset[crystal_table[alt->chan]-1];
        tac_aries_lbl_hist[crystal_table[alt->chan]-1]->Fill(tac_aries_lbl_hist[crystal_table[alt->chan]-1], corrected_tac_value, 1); // tac spectrum per LBL
        tac_aries_lbl_sum->Fill(tac_aries_lbl_sum, corrected_tac_value, 1); // sum tac spectrum including all LBL
        if(alt->ecal >1225 && alt->ecal <1315){ // gate on LaBr3 energy 1275keV

          aries_tac->Fill(aries_tac, corrected_tac_value, 1); // tac spectrum gated on 1275keV
          aries_tac_artEn->Fill(aries_tac_artEn, ptr->e4cal, 1); // ARIES energy spectrum in coincidence with TAC

          if(ptr->e4cal >24 && ptr->e4cal <36){ // gate on ARIES energy
            aries_tac_Egate->Fill(aries_tac_Egate, corrected_tac_value, 1); // tac spectrum gated on 1275keV
          }
        }
      }
    }
    break;
    default: break;
  } // end of inner switch(ALT)
  break; // end of ptr tac

  case SUBSYS_DES_WALL: // DESCANT Wall matrices
  switch(alt->subsys){
    case SUBSYS_DES_WALL: // DSW-DSW
    dt_hist[18]->Fill(dt_hist[18], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
    dt_hist[24]->Fill(dt_hist[24], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);

    c1 = crystal_table[ptr->chan];
    if( c1 >= 1 && c1 <= 60 ){
      c2 = crystal_table[alt->chan];
      if( c2 >= 1 && c2 <= 60 ){
        dsw_hit->Fill(dsw_hit, c1, c2, 1);
        dsw_dsw->Fill(dsw_dsw, (int)ptr->e4cal, (int)alt->e4cal, 1);
        c1--; c2--;

        // Fold 2 sum spectra
        desw_sum_e_nn->Fill(desw_sum_e_nn, (int)ptr->ecal, 1);
        desw_sum_tof_nn->Fill(desw_sum_tof_nn, (int)ptr->e4cal, 1); // e4cal = corrected time-of-flight

        // DSW-DSW angular correlations
        // Fill the appropriate angular bin spectrum with the corrected time-of-flight value
        index = DSW_DSW_angles[c1][c2];
        dsw_dsw_ang_corr_hist[index]->Fill(dsw_dsw_ang_corr_hist[index], (int)ptr->e4cal, (int)alt->e4cal, 1);

        // Fold 2, angle greater than 60 degrees, sum spectra
        // index 13 = 58.555, index 14 = 61.535
        if(index>13){
          desw_sum_e_nn_a->Fill(desw_sum_e_nn_a, (int)ptr->ecal, 1);
          desw_sum_tof_nn_a->Fill(desw_sum_tof_nn_a, (int)ptr->e4cal, 1); // e4cal = corrected time-of-flight
        }
      }
    }
    break;
    case SUBSYS_HPGE: // DSW-HPGe
    // Only use GRGa
    if(output_table[alt->chan] == 1){
      dt_hist[19]->Fill(dt_hist[19], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      ge_dsw->Fill(ge_dsw, (int)alt->ecal, (int)ptr->e4cal, 1);
    }
    break;
    case SUBSYS_ARIES: // DSW-aries
    if(polarity_table[alt->chan] == 0){ // Only use ARIES Fast Output in CAEN
      dt_hist[20]->Fill(dt_hist[20], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      art_dsw->Fill(art_dsw, (int)alt->ecal, (int)ptr->e4cal, 1);

      desw_sum_e_b->Fill(desw_sum_e_b, (int)ptr->ecal, 1);
      desw_sum_tof_b->Fill(desw_sum_tof_b, (int)ptr->e4cal, 1); // e4cal = corrected time-of-flight
    }
    break;
    case SUBSYS_ZDS: // DSW-ZDS
    if(output_table[alt->chan]==0){ // CAEN ZDS
      dt_hist[21]->Fill(dt_hist[21], (int)(abs_dt+DT_SPEC_LENGTH/2), 1);
      dt_hist[25]->Fill(dt_hist[25], (int)(abs(ptr->cfd - alt->cfd)+DT_SPEC_LENGTH/2), 1);

      desw_sum_e_b->Fill(desw_sum_e_b, (int)ptr->ecal, 1);
      desw_sum_tof_b->Fill(desw_sum_tof_b, (int)ptr->e4cal, 1); // e4cal = corrected time-of-flight
    }
    break;
    default: break;
  } // end of inner switch(ALT)
  break; // end of ptr DESCANT Wall matrices

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

extern int read_caen_odb_addresses(int odb_daqsize, unsigned short *addr_table);
int gen_derived_odb_tables()
{
  char sys_name[64], crystal, polarity, type;
  int i, j, tmp, pos, element;

  read_caen_odb_addresses(odb_daqsize, (unsigned short *)addrs);

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
  memset(dtype_subsys,   0xff,  MAX_SUBSYS*sizeof(int));
  memset(subsys_dtype,   0xff,  MAX_SUBSYS*sizeof(int));
  for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
    if( (tmp=sscanf(chan_name[i], "%3c%d%c%c%d%c", sys_name, &pos, &crystal, &polarity, &element, &type)) != 6 ){
      fprintf(stderr,"can't decode name[%s] decoded %d of 6 items\n", chan_name[i], tmp );
      continue;
    }

    // Determine Polarity
    // 1 is N, 0 is P or T or S, -1 is anything else
    if(        polarity == 'N' ){ polarity_table[i] = 1;
    } else if( polarity == 'P' ){ polarity_table[i] = 0;
    } else if( polarity == 'T' ){ polarity_table[i] = 0; // TAC signal
    } else if( polarity == 'S' ){ polarity_table[i] = 1; // ARIES Standard Ouput signal
    } else if( polarity == 'F' ){ polarity_table[i] = 0; // ARIES Fast Output signal
    } else if( polarity == 'X' ){ polarity_table[i] = 0; // XXX type
    } else { fprintf(stderr,"unknown polarity[=%c] in %s\n", polarity, chan_name[i]); }

    // Determine Output
    // Some detector elements have more than one output (HPGe A and B)
    // 1 is A, 0 is B, -1 is X or unknown
    output_table[i] = type=='A' ? 1 : (type=='B' ? 0 : -1);

    // Determine crystal and element position numbers for each Subsystem
    if((strncmp(sys_name,"ART",3) == 0) || (strncmp(sys_name,"DSW",3) == 0) || (strncmp(sys_name,"ZDS",3) == 0) || (strncmp(sys_name,"PAC",3) == 0)
    || (strncmp(sys_name,"LBL",3) == 0) || (strncmp(sys_name,"LBS",3) == 0) || (strncmp(sys_name,"LBT",3) == 0)){ // LBL and LBS, LaBr3 and ancillary BGOs, PAC paces
      crystal_table[i] = pos;
      if(        crystal == 'A' ){ element_table[i] = 1;
      } else if( crystal == 'B' ){ element_table[i] = 2;
      } else if( crystal == 'C' ){ element_table[i] = 3;
      } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3, ZDS, ART
      } else {
        fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
      }
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
    address_chan[ (unsigned short)chan_address[i] ] = i;
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

// Redefine dtype_subsys so we can use it to convert dtype from PSC table to subsys handle index
memset(psc_dtype_subsys,    0xFF,  MAX_SUBSYS*sizeof(int));
for(i=0; i<MAX_SUBSYS; i++){
  psc_dtype_subsys[subsys_dtype[i]] = i;
}

/*
// Print out all the unpacked PSC table information for checking/debugging
fprintf(stdout,"chan\tname\t\taddr\tchan\tdtype\tsubsys\n");
for(i=0; i<odb_daqsize; i++){
fprintf(stdout,"%d\t%s\t0x%04X (%d)\t%d\t%d\t%s\n",i,chan_name[i],chan_address[i],chan_address[i],address_chan[(unsigned short)chan_address[i]],dtype_table[i],subsys_handle[dtype_subsys[dtype_table[i]]]);
}
*/

/*
// Print out all the dtype-subsys tables for checking/debugging
fprintf(stdout,"\nsubsys_dtype\n");
for(i=0; i<MAX_SUBSYS; i++){
fprintf(stdout,"%d=%d\n",i,subsys_dtype[i]);
}
fprintf(stdout,"\ndtype_subsys\n");
for(i=0; i<MAX_SUBSYS; i++){
fprintf(stdout,"%d=%d\n",i,dtype_subsys[i]);
}
fprintf(stdout,"\npsc_dtype_subsys\n");
for(i=0; i<MAX_SUBSYS; i++){
fprintf(stdout,"%d=%d\n",i,psc_dtype_subsys[i]);
}
*/

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
