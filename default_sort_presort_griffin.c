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
#include "default_sort_griffin.h"

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

      // Protect against invalid channel numbers
      if( (unsigned int)chan >= (unsigned int)odb_daqsize ){
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
        init_histos(NULL, ptr->subsys);
      }

      // Check this timestamp against the current cycle bin and fill deadtime histograms if it is a new bin
      if(ppg_cycles_active==1 && ppg_last_ptr_ts>ppg_bin_end && ptr->ts>ppg_last_ptr_ts){
        // Fill deadtime histograms and reset the deadtime count
        bin = (int)((ppg_last_ptr_ts-ppg_cycle_start)/ppg_cycles_binning_factor);  // convert 10ns to binning size set as Global
        if(ppg_cycle_number<MAX_CYCLES){
          //  gea_cycle_num_dt[ppg_cycle_number]->Fill(gea_cycle_num_dt[ppg_cycle_number], bin, subsys_deadtime_count[SUBSYS_HPGE_A]);
          //  geb_cycle_num_dt[ppg_cycle_number]->Fill(geb_cycle_num_dt[ppg_cycle_number], bin, subsys_deadtime_count[SUBSYS_HPGE_B]);
          //  cycle_num_vs_ge_dt->Fill(cycle_num_vs_ge_dt, ppg_cycle_number, bin, subsys_deadtime_count[SUBSYS_HPGE_A]);
          //  cycle_num_vs_ge_b_dt->Fill(cycle_num_vs_ge_b_dt, ppg_cycle_number, bin, subsys_deadtime_count[SUBSYS_HPGE_B]);
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
        //        ppg_cycle_number, ppg_cycle_start, ppg_cycle_end, ppg_cycle_step, ppg_handles[ppg_current_pattern], ppg_pattern_start, ppg_pattern_end);
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
          chan2 = alt->chan;
          if( (unsigned int)chan2 >= (unsigned int)odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
            continue;
          }
          if((dt=ptr->ts - alt->ts)>479 || alt->subsys != SUBSYS_HPGE_A){ continue; }

          if(chan2 != chan ){
            if((clover=(int)(crystal_table[chan2]>>2)) == (int)(crystal_table[chan]>>2)){
              // HPGe Clover time-dependant crosstalk corrections within same clover
              // dt is always positive here
              // The original hit (ptr) came after the crosstalk-inducing hit (alt)
              // Make correction to ptr hit based on energy of alt.
              bin = (int)((1940+dt)/160);
              if(bin<0 || bin>15){ fprintf(stderr,"pre_sort_enter bin [%d] out of bounds for dt %d\n",bin,dt); continue; }
              ge1 = crystal_table[chan];
              c1 = ge1%4;
              c2 = ct_index[c1][(int)(crystal_table[chan2]&3)];
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
          chan2 = alt->chan;
          if( (unsigned int)chan2 >= (unsigned int)odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
            continue;
          }
          if((dt=ptr->ts - alt->ts)>1000 || alt->subsys != SUBSYS_HPGE_A){ continue; }

          if(alt->ecal>1327 && alt->ecal<1337){ // Crosstalk inducing hit was 1332keV
            if(chan2 != chan ){
              if((clover=(int)(crystal_table[chan2]>>2)) == (int)(crystal_table[chan]>>2)){

                // (c2%4) = Crystal Color [B, G, R, W]
                // Hits with ptr arriving after alt
                switch(crystal_table[chan2]&3){
                  case 0:  ct_e_vs_dt_B[crystal_table[chan]]->Fill(ct_e_vs_dt_B[crystal_table[chan]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)ptr->ecal-1100, 1); break;
                  case 1:  ct_e_vs_dt_G[crystal_table[chan]]->Fill(ct_e_vs_dt_G[crystal_table[chan]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)ptr->ecal-1100, 1); break;
                  case 2:  ct_e_vs_dt_R[crystal_table[chan]]->Fill(ct_e_vs_dt_R[crystal_table[chan]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)ptr->ecal-1100, 1); break;
                  case 3:  ct_e_vs_dt_W[crystal_table[chan]]->Fill(ct_e_vs_dt_W[crystal_table[chan]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)ptr->ecal-1100, 1); break;
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
      //    all other events are later than frag_idx
      //    check other frags in window for possible suppression and/or summing
      //  also calculate multiplicities[store in frag_idx only]
      int pre_sort_exit(int frag_idx, int end_idx)
      {
        Grif_event *alt2, *alt, *ptr = &grif_event[frag_idx];
        float desw_median_distance = 1681.8328; // descant wall median source-to-detector distance in mm
        int i, j, dt, dt13, tof;
        float q1,integ2,q12,k1,k2,k12,e1,e2,e12,m,c;
        int chan,chan2,found,pos;
        int clover, ge1, c1,c2,bin, p_strip, n_strip;
        float energy,ecal,correction,angle;

        // Assign chan local variable and check it is a valid channel number
        chan = ptr->chan;
        if( (unsigned int)chan >= (unsigned int)odb_daqsize ){
          fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
          return(-1);
        }
        i = frag_idx; ptr->multiplicity = 1;
        while( i != end_idx ){ // need at least two events in window
          if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
          chan2 = alt->chan;
          if( (unsigned int)chan2 >= (unsigned int)odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
            continue;
          }

          // Determine absolute time difference between timestamps
          dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }

          // Restrict to 2 microseconds presort window for everything except HPGe crosstalk to maintain speed
          if(alt->subsys != SUBSYS_HPGE_A && dt>250){ continue; }

          // Determine multiplicity
          if( alt->subsys == ptr->subsys ){
            if(ptr->subsys == SUBSYS_HPGE_A){
              // Prompt coincidence window for HPGe multiplicity
              if( (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) & (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) ){
                ++ptr->multiplicity;
              }
            }else if(ptr->subsys == SUBSYS_QED_STRIP){
              // Same DSSD and Prompt coincidence window for QED multiplicity
              if( (crystal_table[ptr->chan] == crystal_table[alt->chan]) && (dt >= qed_fb_window_min) && (dt <= qed_fb_window_max) ){
                ++ptr->multiplicity;
              }
            }else{ // 2 microseconds presort window for multiplicity of all other subsystem types
              ++ptr->multiplicity;
            }

          }

          // SubSystem-specific pre-processing
          switch(ptr->subsys){
            case SUBSYS_HPGE_A:

            // HPGe Clover time-dependant crosstalk corrections
            if(alt->subsys == SUBSYS_HPGE_A && chan2 != chan){
              if((clover=(int)(crystal_table[chan2]>>2)) == (int)(crystal_table[chan]>>2)){
                dt = ptr->ts - alt->ts; // Use relative time difference. This is always negative.

                if(dt>-1940 && dt <= 0){
                  // The original hit (ptr) came earlier then the crosstalk-inducing hit (alt)
                  // Make correction to ptr hit based on energy of alt.
                  //  dt -= 1920; // Correct the timestamp difference so that the right correction is calculated
                  bin = (int)((1940+dt)/160);
                  if(bin<0 || bin>15){ fprintf(stderr,"pre_sort_exit bin [%d] out of bounds for dt %d\n",bin,dt); continue; }
                  ge1 = crystal_table[chan];
                  c1 = ge1&3;
                  c2 = ct_index[c1][(int)(crystal_table[chan2]&3)];
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
              //  gea_self_dt->Fill(gea_self_dt, dt, crystal_table[chan], 1);
              perform_pileup_correction(ptr, alt, dt, chan, chan2, i, end_idx);
            }
            // BGO suppression of HPGe
            if(alt->subsys == SUBSYS_BGO){
              if( (dt >= bgo_window_min && dt <= bgo_window_max) ){
                // could alternatively use crystal numbers rather than clover#
                //    (don't currently have this for BGO)
                if( (int)(crystal_table[ptr->chan]>>2) == (int)(crystal_table[alt->chan]>>2) ){ ptr->suppress = 1; }
              }
            }
            // Germanium addback -
            //    earliest fragment has the sum energy, others are marked -1
            // Remember the other crystal channel number in alt_chan for use in Compton Polarimetry
            if( (dt >= addback_window_min && dt <= addback_window_max) && alt->subsys == SUBSYS_HPGE_A ){
              if( alt->esum >= 0 && (int)(crystal_table[alt->chan]>>2) == (int)(crystal_table[ptr->chan]>>2) ){
                ptr->esum += alt->esum; alt->esum = -1; ptr->alt_chan = alt->chan; if(alt->suppress==1){ ptr->suppress = 1; }
              }
            }
            // Look to build Compton events (a QED-Ge coincidence with energy of 511keV)
            // Just tag the Ge information into the strip as this Ge leaves the presort window
            if( alt->subsys == SUBSYS_QED_STRIP ){
              if(alt->ecal>QED_PIXEL_THRESHOLD && ((ptr->esum+alt->ecal) > QED_GAMMA_ENERGY-QED_GAMMA_ENERGY_WINDOW) && ((ptr->esum+alt->ecal) < QED_GAMMA_ENERGY+QED_GAMMA_ENERGY_WINDOW) ){
                if(dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL] && dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL]){
                  // Here ptr is Ge and alt is QED_PIXEL
                  // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL, esum will be the full energy
                  // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
                  alt->alt_ecal = ptr->esum;
                  alt->esum += ptr->esum; // Add this Ge energy to the pixel sum energy
                  alt->net_id = crystal_table[ptr->chan];
                  alt->delta_t = ptr->ts-alt->ts;
                  if(ptr->esum>ptr->ecal){
                    alt->alt2_ecal = ptr->esum-ptr->ecal;
                    alt->alt2_chan = crystal_table[ptr->alt_chan];
                  }else{
                    alt->alt2_ecal = 0;
                    alt->alt2_chan = -1;
                  }
                }
              }
            }
            // Tag HPGe events which have a beta-coincidence
            // Use the ptr->tof
            // ptr->tof = 0 = No beta coincidence
            // ptr->tof > 0 = Beta coincidence with any beta detector
            // ptr->tof = 1 = Beta coincidence with SCEPTAR only
            // ptr->tof = 2 = Beta coincidence with ZDS only
            // ptr->tof = 4 = Beta coincidence with ARIES only
            // The bit assignments are cumulative eg. pt->tof = 5 = Beta coincidence with both SCEPTAR and ARIES.
            if( alt->subsys == SUBSYS_SCEPTAR && (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) && (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) ){ ptr->tof = ptr->tof | 1; break; }
            if( alt->subsys == SUBSYS_ZDS_A   && (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ZDS_A])   && (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ZDS_A])   ){ ptr->tof = ptr->tof | 2; break; }
            if( alt->subsys == SUBSYS_ARIES_A && (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) && (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) ){
              ptr->tof = ptr->tof | 4;
              pos = crystal_table[alt->chan];
              if((pos>12 && pos<17) || (pos>64 && pos<69)){
                // Triangles
                ptr->tof = ptr->tof | 8;
              }else if(pos<25 || pos>56){
                // Squares
                ptr->tof = ptr->tof | 16;
              }else{
                // Rectangles
                ptr->tof = ptr->tof | 32;
              }break;
            }
            break;
            case SUBSYS_SCEPTAR:
            if(alt->subsys == SUBSYS_HPGE_A && (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) && (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) ){ alt->tof = alt->tof | 1; }
            if(alt->subsys == SUBSYS_PACES  && (dt >= time_diff_gate_min[SUBSYS_PACES][SUBSYS_SCEPTAR])  && (dt <= time_diff_gate_max[SUBSYS_PACES][SUBSYS_SCEPTAR])  ){ alt->tof = alt->tof | 1; }
            break;
            case SUBSYS_ZDS_A:
            if(alt->subsys == SUBSYS_HPGE_A && (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ZDS_A])   && (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ZDS_A])   ){ alt->tof = alt->tof | 2; }
            if(alt->subsys == SUBSYS_PACES  && (dt >= time_diff_gate_min[SUBSYS_PACES][SUBSYS_ZDS_A])    && (dt <= time_diff_gate_max[SUBSYS_PACES][SUBSYS_ZDS_A])    ){ ptr->tof = alt->tof | 2; }
            break;
            case SUBSYS_ARIES_A:
            if(alt->subsys == SUBSYS_PACES  && (dt >= time_diff_gate_min[SUBSYS_PACES][SUBSYS_ARIES_A])  && (dt <= time_diff_gate_max[SUBSYS_PACES][SUBSYS_ARIES_A])  ){ alt->tof = alt->tof | 4; }
            if(alt->subsys == SUBSYS_HPGE_A && (dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) && (dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) ){
              alt->tof = alt->tof | 4;
              pos = crystal_table[alt->chan];
              if((pos>12 && pos<17) || (pos>64 && pos<69)){
                // Triangles
                alt->tof = alt->tof | 8;
              }else if(pos<25 || pos>56){
                // Squares
                alt->tof = alt->tof | 16;
              }else{
                // Rectangles
                alt->tof = alt->tof | 32;
              }
            }
            break;
            case SUBSYS_HPGE_B:
            // HPGe B pile-up corrections
            if(alt->subsys == SUBSYS_HPGE_B && chan2 == chan){
              //  geb_self_dt->Fill(geb_self_dt, dt, crystal_table[chan], 1);
              perform_pileup_correction(ptr, alt, dt, chan, chan2, i, end_idx);
            }
            break;
            case SUBSYS_PACES:
            // Tag PACES events which have a beta-coincidence
            // Use the ptr->tof
            // ptr->tof = 0 = No beta coincidence
            // ptr->tof > 0 = Beta coincidence with any beta detector
            // ptr->tof = 1 = Beta coincidence with SCEPTAR only
            // ptr->tof = 2 = Beta coincidence with ZDS only
            // ptr->tof = 4 = Beta coincidence with ARIES only
            // The bit assignments are cumulative eg. pt->tof = 5 = Beta coincidence with both SCEPTAR and ARIES.
            if( alt->subsys == SUBSYS_SCEPTAR && (dt >= time_diff_gate_min[SUBSYS_PACES][SUBSYS_SCEPTAR]) && (dt <= time_diff_gate_max[SUBSYS_PACES][SUBSYS_SCEPTAR]) ){ ptr->tof = ptr->tof | 1; break; }
            if( alt->subsys == SUBSYS_ZDS_A   && (dt >= time_diff_gate_min[SUBSYS_PACES][SUBSYS_ZDS_A])   && (dt <= time_diff_gate_max[SUBSYS_PACES][SUBSYS_ZDS_A])   ){ ptr->tof = ptr->tof | 2; break; }
            if( alt->subsys == SUBSYS_ARIES_A && (dt >= time_diff_gate_min[SUBSYS_PACES][SUBSYS_ARIES_A]) && (dt <= time_diff_gate_max[SUBSYS_PACES][SUBSYS_ARIES_A]) ){ ptr->tof = ptr->tof | 4; break; }
            break;
            case SUBSYS_RCMP:
            // RCMP Front-Back coincidence
            // Ensure its the same DSSD and that the two events are front and back
            // The charged particles enter the P side and this has superior energy resolution
            // Ensure the energy collected in the front and back is similar
            ptr->esum = -1; // Need to exclude any noise and random coincidences.
            if( alt->subsys == SUBSYS_RCMP && (dt >= rcmp_fb_window_min && dt <= rcmp_fb_window_max) && (ptr->ecal>0 && ptr->ecal<32768)){
              if((crystal_table[ptr->chan] == crystal_table[alt->chan]) && (polarity_table[ptr->chan] != polarity_table[alt->chan]) && (alt->ecal > 0  && alt->ecal<32768)){
                if( ((ptr->ecal / alt->ecal)<=1.1 && (ptr->ecal / alt->ecal)>=0.9)){
                  // Ensure esum comes from P side, but use this timestamp
                  ptr->esum = polarity_table[ptr->chan]==0 ? ptr->ecal : (polarity_table[ptr->chan]==1 ? alt->ecal : -1);
                  ptr->suppress = alt->suppress = 1;
                }
              }
            }
            break;
            case SUBSYS_QED_STRIP:
            // QED Front-Back coincidence of two QED strips will define a QED pixel
            // Ensure its the same DSSD and that the two events are front and back strips
            // The charged particles enter the P side and this has superior energy resolution
            // Ensure the energy collected in the front and back is similar
            if( alt->subsys == SUBSYS_QED_STRIP && (dt >= qed_fb_window_min && dt <= qed_fb_window_max) && (ptr->ecal>QED_STRIP_THRESHOLD && ptr->ecal<32768)){
              if((crystal_table[ptr->chan] == crystal_table[alt->chan]) && (polarity_table[ptr->chan] != polarity_table[alt->chan]) && (alt->ecal > QED_STRIP_THRESHOLD && alt->ecal<32768)){
                //  if( ((ptr->ecal / alt->ecal)<=1.1 && (ptr->ecal / alt->ecal)>=0.9)){ // Energy-sharing only works if both strips are calibrated!
                // The ptr strip now changes to a PIXEL
                // Ensure ecal comes from P side, alt_ecal will be N side
                // Use the timestamp of ptr which came earlier
                // 32*32 = 1024 pixels in each DSSD.
                // crystal_table[ptr->chan] is DSSD number [1-6]
                // ptr->alt_chan is the pixel number within this DSSD [0-1023]
                // P strip number is Floor(pixel number / N_QED_STRIPS) [0-31]
                // N strip number is pixel number % N_QED_STRIPS [0-31]
                if(polarity_table[ptr->chan]==0){ // ptr is P strip, alt is N strip

                  if((unsigned int)ptr->net_id < 64){
                    angle = scattering_angle_QEDGe(crystal_table[ptr->chan],(((element_table[ptr->chan] * N_QED_STRIPS) + element_table[alt->chan])%1024),ptr->net_id);
                  }else{ angle=0; }
                  pos = crystal_table[ptr->chan]-1;
                  p_strip = element_table[ptr->chan];
                  n_strip = element_table[alt->chan];

                  if(ptr->alt_ecal>0 && (angle>=compton_angle(ptr->alt_ecal,QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle<=compton_angle(ptr->alt_ecal,QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){ // Already have a Ge tag to this strip
                    // esum and alt_ecal already set to P strip + Ge values
                    ptr->subsys = SUBSYS_COMPTON;

                    if(ptr->alt2_ecal>0){
                      ptr->subsys = SUBSYS_DCOMPTONA;
                    }
                    // Fill P and N strip energy calibration histograms here before N strip (and energy) is discarded
                    // Single strip vs theta needed for strip-energy calibration
                    qedp_ge_theta[p_strip + pos*N_QED_STRIPS]->Fill(qedp_ge_theta[p_strip + pos*N_QED_STRIPS], ptr->ecal, angle, 1);
                    qedn_ge_theta[n_strip + pos*N_QED_STRIPS]->Fill(qedn_ge_theta[n_strip + pos*N_QED_STRIPS], alt->ecal, angle, 1);
                  }else{

                    // PIXEL identified
                    ptr->esum = ptr->ecal;
                    ptr->alt_ecal = alt->ecal;
                    ptr->subsys = SUBSYS_QED_PIXEL;
                  }
                  ptr->alt_chan = (element_table[ptr->chan] * N_QED_STRIPS) + element_table[alt->chan];
                }else{ // alt is P strip, ptr is N strip

                  if((unsigned int)ptr->net_id < 64){
                    angle = scattering_angle_QEDGe(crystal_table[ptr->chan],(((element_table[alt->chan] * N_QED_STRIPS) + element_table[ptr->chan])%1024),ptr->net_id);
                  }else{ angle=0; }
                  pos = crystal_table[alt->chan]-1;
                  p_strip = element_table[alt->chan];
                  n_strip = element_table[ptr->chan];
                  if(ptr->alt_ecal>0 && (angle>=compton_angle(ptr->alt_ecal,QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle<=compton_angle(ptr->alt_ecal,QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){ // Already have a Ge tag to this strip but it is the N
                    // alt_ecal already set Ge value. Change esum to use P strip energy
                    ptr->esum = ptr->alt_ecal + alt->ecal;
                    ptr->subsys = SUBSYS_COMPTON;
                    if(ptr->alt2_ecal>0){
                      ptr->subsys = SUBSYS_DCOMPTONA;
                    }
                    // Fill P and N strip energy calibration histograms here before N strip (and energy) is discarded
                    // Single strip vs theta needed for strip-energy calibration
                    qedp_ge_theta[p_strip + pos*N_QED_STRIPS]->Fill(qedp_ge_theta[p_strip + pos*N_QED_STRIPS], alt->ecal, angle, 1);
                    qedn_ge_theta[n_strip + pos*N_QED_STRIPS]->Fill(qedn_ge_theta[n_strip + pos*N_QED_STRIPS], ptr->ecal, angle, 1);
                  }else{
                    ptr->alt_ecal = ptr->ecal;
                    ptr->esum = ptr->ecal = alt->ecal;
                    ptr->subsys = SUBSYS_QED_PIXEL;
                  }
                  ptr->alt_chan = (element_table[alt->chan] * N_QED_STRIPS) + element_table[ptr->chan];
                }
                --ptr->multiplicity; // Adjust multiplicity because we are combining two strips into one pixel.
                //  }
              }
            }
            break;
            case SUBSYS_QED_PIXEL:
            // Check for any additional neighbouring strips in this DSSD to use in addback
            // If a neighbour strip has higher energy then change the pixel to that one, otherwise just add the energy
            if( alt->subsys == SUBSYS_QED_STRIP && (dt >= qed_fb_window_min && dt <= qed_fb_window_max) && (alt->ecal>QED_STRIP_THRESHOLD && alt->ecal<32768)){
              if(crystal_table[ptr->chan] == crystal_table[alt->chan]){ // same DSSD


                if(polarity_table[alt->chan]==0){ // alt is P strip
                  if( abs((ptr->alt_chan / N_QED_STRIPS) - element_table[alt->chan]) == 1 ){ // neighbour strip
                    if(alt->ecal > ptr->ecal){
                      // Change pixel number to use this alt strip
                      ptr->alt_chan = (element_table[ptr->chan] * N_QED_STRIPS) + element_table[ptr->alt_chan% N_QED_STRIPS];
                      ptr->ecal += alt->ecal;
                    }else{ // just add the energies
                      ptr->ecal += alt->ecal;
                    }
                  }
                }else{ // alt is N strip
                  // This could change the pixel number, but pixel energy came from p strips
                  if(alt->ecal > ptr->alt_ecal){
                    // Change pixel number to use this alt strip
                    ptr->alt_chan = (element_table[ptr->alt_chan/N_QED_STRIPS] * N_QED_STRIPS) + element_table[ptr->chan];
                    ptr->alt_ecal += alt->ecal;
                  }else{
                    ptr->alt_ecal += alt->ecal;
                  }
                }
              }
            }

            // Look to build Compton events (a QED-Ge coincidence with energy of 511keV)
            if( alt->subsys == SUBSYS_HPGE_A ){
              if(ptr->ecal>QED_PIXEL_THRESHOLD && ptr->ecal+alt->esum > QED_GAMMA_ENERGY-QED_GAMMA_ENERGY_WINDOW && ptr->ecal+alt->esum < QED_GAMMA_ENERGY+QED_GAMMA_ENERGY_WINDOW){
                if(dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL] && dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL]){
                  angle = scattering_angle_QEDGe(crystal_table[ptr->chan],(ptr->alt_chan&1023),crystal_table[alt->chan]);
                  if( (angle>=compton_angle(alt->esum,QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle<=compton_angle(alt->esum,QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW) ){

                    pos = crystal_table[ptr->chan]-1;
                    p_strip = (int)((ptr->alt_chan&1023)/N_QED_STRIPS);
                    n_strip = (ptr->alt_chan&1023)%N_QED_STRIPS;
                    // Fill P and N strip energy calibration histograms here before N strip (and energy) is discarded
                    // Single strip vs theta needed for strip-energy calibration
                    qedp_ge_theta[p_strip + pos*N_QED_STRIPS]->Fill(qedp_ge_theta[p_strip + pos*N_QED_STRIPS], ptr->ecal, angle, 1);
                    qedn_ge_theta[n_strip + pos*N_QED_STRIPS]->Fill(qedn_ge_theta[n_strip + pos*N_QED_STRIPS], ptr->alt_ecal, angle, 1);

                    // Here alt is Ge and ptr is QED_PIXEL
                    // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL
                    // In COMPTON event, chan will be Ge, alt_chan will be packed with DSSD and PIXEL numbers
                    // OR
                    // 13 bits needed for array pixel number [0-6143] (Can just be treated as a decimal)
                    ptr->alt_ecal = alt->esum;
                    ptr->esum += alt->esum; // Add this Ge energy to the pixel energy
                    ptr->net_id = crystal_table[alt->chan];
                    ptr->delta_t = alt->ts-ptr->ts;
                    if(alt->esum>alt->ecal){ // Ge with addback is a Double Compton scatter (DSSD-Ge-Ge)
                      ptr->subsys = SUBSYS_DCOMPTONA;
                      ptr->alt2_ecal = (alt->esum-alt->ecal);
                      ptr->alt2_chan = crystal_table[alt->alt_chan];
                    }else{ // Single Ge crsytal is a Compton scatter (DSSD-Ge)
                      ptr->subsys = SUBSYS_COMPTON;
                      ptr->alt2_ecal = 0;
                      ptr->alt2_chan = -1;
                    }
                  }
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
                //  desw_psd_zdse->Fill(desw_psd_zdse, (int)alt->psd, (int)ptr->ecal, 1); // Fill desw_psd_zdse
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
                if((clover=(int)(crystal_table[chan2]>>2)) == (int)(crystal_table[chan]>>2)){

                  // Fill crosstalk histograms
                  // (c2%4) = Crystal Color [B, G, R, W]
                  // Hits with ptr arriving after alt
                  // Fill higher-channels of the x axis on these plots
                  switch(crystal_table[chan]%4){
                    case 0:  ct_e_vs_dt_B[crystal_table[chan2]]->Fill(ct_e_vs_dt_B[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)alt->ecal-1100, 1); break;
                    case 1:  ct_e_vs_dt_G[crystal_table[chan2]]->Fill(ct_e_vs_dt_G[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)alt->ecal-1100, 1); break;
                    case 2:  ct_e_vs_dt_R[crystal_table[chan2]]->Fill(ct_e_vs_dt_R[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)alt->ecal-1100, 1); break;
                    case 3:  ct_e_vs_dt_W[crystal_table[chan2]]->Fill(ct_e_vs_dt_W[crystal_table[chan2]], ((int)((alt->ts - ptr->ts)>>2)+300), (int)alt->ecal-1100, 1); break;
                  }
                }
              }
            }
          }// end of while
        }
        return(0);
      }


      #define TRIPLES_GE 0
      #define TRIPLES_QED 1
      // Presort - The final presort to identify triple coincidence events
      //  - frag_idx is about to leave coinc window (which ends at end_idx)
      //    all other events are later than frag_idx
      int pre_sort_triples(int frag_idx, int end_idx)
      {
        Grif_event *tmp, *alt2, *alt, *ptr = &grif_event[frag_idx];
        int i, j, dt, dt13, tof;
        int pos1, pos2, qed1, qed2;
        float q1,integ2,q12,k1,k2,k12,e1,e2,e12,m,c;
        int chan,chan2,chan3,found,pos, ptr_swap=0, alt_swap=0;
        int clover, ge1, c1,c2,bin, p_strip, n_strip;
        float energy,ecal,correction,angle1,angle2;
        int triples_multiplicity[2] = {0,0}; // subsys: Ge is 0, PIXEL is 22
        int debug_triples=0;

        // Only interested in HPGe - PIXEL coincidences in this presort
        if(ptr->subsys != SUBSYS_HPGE_A && ptr->subsys != SUBSYS_QED_PIXEL){ return(0); }
        triples_multiplicity[(int)(ptr->subsys/22)]++;

        // Assign chan local variable and check it is a valid channel number
        chan = ptr->chan;
        if( (unsigned int)chan >= (unsigned int)odb_daqsize ){
          fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
          return(-1);
        }
        i = frag_idx;
        while( i != end_idx ){ // need at least two events in window
          if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
          chan2 = alt->chan;
          if( (unsigned int)chan2 >= (unsigned int)odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan2:%d\n",chan2 );
            continue;
          }
          // Determine absolute time difference between timestamps
          dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }
          if(dt>40){ continue; } // Limit to +/- 400 nanoseconds
          //  fprintf(stdout,"Scanning Doubles: %d %d: %d: %d %d, %.1f %.1f (%.1f)\n",frag_idx,i,(ptr->ts - alt->ts),ptr->subsys,alt->subsys,ptr->ecal,alt->ecal,(ptr->ecal+alt->ecal));

          // Only interested in HPGe - PIXEL coincidences in this presort
          if(alt->subsys != SUBSYS_HPGE_A && alt->subsys != SUBSYS_QED_PIXEL){ continue; }
          if(alt->ecal == ptr->ecal){ continue; }
          triples_multiplicity[(int)(alt->subsys/22)]++;
          //  fprintf(stdout,"Accepted Doubles: %d %d: %d: %d %d, %.1f %.1f (%.1f)\n",frag_idx,i,(ptr->ts - alt->ts),ptr->subsys,alt->subsys,ptr->ecal,alt->ecal,(ptr->ecal+alt->ecal));

          while( i != end_idx ){ // need at least three events in window
            if( ++i >=  PTR_BUFSIZE ){ i=0; } alt2 = &grif_event[i]; // WRAP
            chan3 = alt2->chan;
            if( (unsigned int)chan3 >= (unsigned int)odb_daqsize ){
              fprintf(stderr,"presort error: ignored event in chan3:%d\n",chan3 );
              continue;
            }
            // Determine absolute time difference between timestamps
            dt13 = ptr->ts - alt2->ts; if( dt13 < 0 ){ dt13 = -1*dt13; }
            if(dt13>40){ continue; } // Limit to +/- 400 nanoseconds
            //    fprintf(stdout,"Triples ALL %d %d: %d %d: %d %d %d, %.1f %.1f %.1f (%.1f)\n",frag_idx,i,dt,dt13,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal));

            // Only interested in HPGe - PIXEL coincidences in this presort
            if(alt2->subsys != SUBSYS_HPGE_A && alt2->subsys != SUBSYS_QED_PIXEL){ continue; }
            if(alt2->ecal == ptr->ecal){ continue; }
            if(alt2->ecal == alt->ecal){ continue; }
            triples_multiplicity[(int)(alt2->subsys/22)]++;

            //    fprintf(stdout,"Triples presort ONLY Ge+Si %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13);
            //  if(ptr->subsys == SUBSYS_QED_PIXEL && alt->subsys == SUBSYS_QED_PIXEL){
            //    fprintf(stdout,"Triples presort Ge+Si ONLY %d %d: %d %d: %d %d %d, %.1f %.1f %.1f (%.1f)\n",frag_idx,i,dt,dt13,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal));
            //  }

            // Sort fragments by subsys not ts. ptr now can only be HPGE_A
            //  if( alt->subsys > alt2->subsys ){ tmp = alt; alt = alt2; alt2 = tmp; alt_swap = 1; }
            //  if( ptr->subsys > alt->subsys ){ tmp = ptr; ptr = alt; alt = tmp; ptr_swap = 1; }

            //  fprintf(stdout,"Triples presort REORDERED %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13);

            // Only interested in HPGe - PIXEL coincidences in this presort
            //  if(ptr->subsys != SUBSYS_HPGE_A || alt->subsys != SUBSYS_QED_PIXEL || alt2->subsys != SUBSYS_QED_PIXEL){ continue; }

            //fprintf(stdout,"Triples Multiplicity (GE,QED): [%d,%d]\n",triples_multiplicity[0],triples_multiplicity[1]);
            if(triples_multiplicity[TRIPLES_GE] == 1 && triples_multiplicity[TRIPLES_QED] == 2){
              if((ptr->ecal+alt->ecal+alt2->ecal)>495 && (ptr->ecal+alt->ecal+alt2->ecal)<527){
                if(debug_triples){fprintf(stdout,"================Triples presort IDENTIFIED %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13);}
                // QED DCOMPTONA EVENTS
                // Ge with addback is a Double Compton scatter (DSSD-Ge-Ge)
                // Identified as subsys==SUBSYS_DCOMPTONA
                // In DCOMPTONA event, ecal will be QED_PIXEL and alt_ecal will be Ge addback sum energy
                // In DCOMPTONA event, crystal_table[chan]=pos will be QED DSSD number [1-6], alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
                // In DCOMPTONA event, net_id will be first HPGe crystal number [1-64], alt2_chan will be second HPGe crystal number [1-64]

                // QED DCOMPTONB EVENTS
                // Ge with addback is a Double Compton scatter (DSSD-DSSD-Ge)
                // Identified as subsys==SUBSYS_DCOMPTONB
                // In DCOMPTONB event, ecal will be the first QED_PIXEL, alt_ecal will be the second QED_PIXEL, alt2_ecal will be Ge energy
                // In DCOMPTONB event, crystal_table[chan]=pos will be  first QED DSSD number [1-6],  alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
                // In DCOMPTONB event, crystal_table[ tof]=pos will be second QED DSSD number [1-6], alt2_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
                // In DCOMPTONB event, net_id will be the HPGe crystal number [1-64]

                // Need to identify the correct order
                // First, which is the HPGe: use subsys. This is the final secondary photon.
                //
                if(ptr->subsys == 0){
                  // ptr is SUBSYS_HPGE_A
                  // initial energy is 511keV

                  pos1  = crystal_table[alt->chan]; // QED DSSD number [1-6]
                  qed1 = alt->alt_chan; // QED pixel number [0-1023]
                  pos2  = crystal_table[alt2->chan]; // QED DSSD number [1-6]
                  qed2 = alt2->alt_chan; // QED pixel number [0-1023]
                  angle1 = scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
                  angle2 = scattering_angle_QEDQED(pos2, qed2, pos1, qed1);
                  if((angle1>=compton_angle((ptr->ecal+alt->ecal),QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle1<=compton_angle((ptr->ecal+alt->ecal),QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
                    // alt is first Si hit, ptr is Ge.
                    if(debug_triples){fprintf(stdout,"alt2 is first Si hit, alt is second Si hit, ptr is Ge\n");}
                    alt2->subsys = SUBSYS_DCOMPTONB;                 // First Si hit becomes DCOMPTONB event
                    alt2->alt_ecal  = alt->ecal;                     // Second QED pixel energy
                    alt2->alt2_ecal = ptr->ecal;                     // Ge energy
                    alt2->esum = ptr->ecal + alt->ecal + alt2->ecal; // Add this Ge energy to the pixel energy
                    alt2->tof = alt->chan;                           // Second QED DSSD number [1-6]
                    alt2->alt2_chan = alt->alt_chan;                 // Second QED pixel number [DSSD*PIXELnumber],[0-5*0-1023]
                    alt2->net_id = crystal_table[ptr->chan];         // Ge crystal number [1-64]
                    alt2->delta_t = alt2->ts-ptr->ts;                // Time difference covering all three Hits
                  }else{
                    if((angle2>=compton_angle((ptr->ecal+alt2->ecal),QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle2<=compton_angle((ptr->ecal+alt2->ecal),QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
                      // alt2 is first Si hit, ptr is Ge.
                      if(debug_triples){fprintf(stdout,"alt is first Si hit, alt2 is second Si hit, ptr is Ge\n");}
                      alt->subsys = SUBSYS_DCOMPTONB;                 // First Si hit becomes DCOMPTONB event
                      alt->alt_ecal  = alt2->ecal;                    // Second QED pixel energy
                      alt->alt2_ecal = ptr->ecal;                     // Ge energy
                      alt->esum = ptr->ecal + alt->ecal + alt2->ecal; // Add this Ge energy to the pixel energy
                      alt->tof = alt2->chan;                          // Second QED DSSD number [1-6]
                      alt->alt2_chan = alt2->alt_chan;                // Second QED pixel number [DSSD*PIXELnumber],[0-5*0-1023]
                      alt->net_id = crystal_table[ptr->chan];         // Ge crystal number [1-64]
                      alt->delta_t = alt2->ts-ptr->ts;                // Time difference covering all three Hits
                    }else{
                      // No solution
                      if(debug_triples){fprintf(stdout,"No solution for this with ptr=Ge: %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d, [%.1f %.1f], [%.1f %.1f]\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13,angle1,(ptr->ecal+alt->ecal),angle2,(ptr->ecal+alt2->ecal));}
                    }
                  }

                }else if(alt->subsys == 0){
                  // alt is SUBSYS_HPGE_A

                  pos1  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
                  qed1 = ptr->alt_chan; // QED pixel number [0-1023]
                  pos2  = crystal_table[alt2->chan]; // QED DSSD number [1-6]
                  qed2 = alt2->alt_chan; // QED pixel number [0-1023]
                  angle1 = scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
                  angle2 = scattering_angle_QEDQED(pos2, qed2, pos1, qed1);
                  if((angle1>=compton_angle((alt->ecal+ptr->ecal),QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle1<=compton_angle((alt->ecal+ptr->ecal),QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
                    // alt2 is first Si hit, alt is Ge.
                    if(debug_triples){fprintf(stdout,"alt2 is first Si hit, ptr is second Si hit, alt is Ge\n");}
                    alt2->subsys = SUBSYS_DCOMPTONB;                 // First Si hit becomes DCOMPTONB event
                    alt2->alt_ecal  = ptr->ecal;                     // Second QED pixel energy
                    alt2->alt2_ecal = alt->ecal;                     // Ge energy
                    alt2->esum = ptr->ecal + alt->ecal + alt2->ecal; // Add this Ge energy to the pixel energy
                    alt2->tof = ptr->chan;                           // Second QED DSSD number [1-6]
                    alt2->alt2_chan = ptr->alt_chan;                 // Second QED pixel number [DSSD*PIXELnumber],[0-5*0-1023]
                    alt2->net_id = crystal_table[alt->chan];         // Ge crystal number [1-64]
                    alt2->delta_t = alt2->ts-ptr->ts;                // Time difference covering all three Hits
                  }else{
                    if((angle2>=compton_angle((alt->ecal+alt2->ecal),QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle2<=compton_angle((alt->ecal+alt2->ecal),QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
                      // ptr is first Si hit, alt is Ge.
                      if(debug_triples){fprintf(stdout,"ptr is first Si hit, alt2 is second Si hit, alt is Ge\n");}
                      ptr->subsys = SUBSYS_DCOMPTONB;                 // First Si hit becomes DCOMPTONB event
                      ptr->alt_ecal  = alt2->ecal;                    // Second QED pixel energy
                      ptr->alt2_ecal = alt->ecal;                     // Ge energy
                      ptr->esum = ptr->ecal + alt->ecal + alt2->ecal; // Add this Ge energy to the pixel energy
                      ptr->tof = alt2->chan;                          // Second QED DSSD number [1-6]
                      ptr->alt2_chan = alt2->alt_chan;                // Second QED pixel number [DSSD*PIXELnumber],[0-5*0-1023]
                      ptr->net_id = crystal_table[alt->chan];         // Ge crystal number [1-64]
                      ptr->delta_t = alt2->ts-ptr->ts;                // Time difference covering all three Hits
                    }else{
                      // No solution
                      if(debug_triples){fprintf(stdout,"No solution for this with alt=Ge: %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d, [%.1f %.1f], [%.1f %.1f]\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13,angle1,(alt->ecal+ptr->ecal),angle2,(alt->ecal+alt2->ecal));}
                    }
                  }
                }else{
                  // alt2 is SUBSYS_HPGE_A
                  // initial energy is 511keV
                  // Secondary photon energy is 511keV minus the first DSSD energy
                  // Secondary photon energy is either 511.0-ptr->ecal OR 511.0-alt->ecal
                  // Secondary photon energy is either ptr->ecal+alt2->ecal OR alt->ecal+alt2->ecal
                  // Check condition for these secondary energies for each DSSD pixel angle

                  pos1  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
                  qed1 = ptr->alt_chan; // QED pixel number [0-1023]
                  pos2  = crystal_table[alt->chan]; // QED DSSD number [1-6]
                  qed2 = alt->alt_chan; // QED pixel number [0-1023]
                  angle1 = scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
                  angle2 = scattering_angle_QEDQED(pos2, qed2, pos1, qed1);
                  if((angle1>=compton_angle((alt2->ecal+alt->ecal),QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle1<=compton_angle((alt2->ecal+alt->ecal),QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
                    // ptr is first Si hit, alt2 is Ge.
                    if(debug_triples){fprintf(stdout,"ptr is first Si hit, alt is second Si hit, alt2 is Ge\n");}
                    ptr->subsys = SUBSYS_DCOMPTONB;                 // First Si hit becomes DCOMPTONB event
                    ptr->alt_ecal  = alt->ecal;                     // Second QED pixel energy
                    ptr->alt2_ecal = alt2->ecal;                    // Ge energy
                    ptr->esum = ptr->ecal + alt->ecal + alt2->ecal; // Add this Ge energy to the pixel energy
                    ptr->tof = alt->chan;                           // Second QED DSSD number [1-6]
                    ptr->alt2_chan = alt->alt_chan;                 // Second QED pixel number [DSSD*PIXELnumber],[0-5*0-1023]
                    ptr->net_id = crystal_table[alt2->chan];        // Ge crystal number [1-64]
                    ptr->delta_t = alt2->ts-ptr->ts;                // Time difference covering all three Hits
                  }else{
                    if((angle2>=compton_angle((alt2->ecal+ptr->ecal),QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle2<=compton_angle((alt2->ecal+ptr->ecal),QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
                      // alt is first Si hit, alt2 is Ge.
                      if(debug_triples){fprintf(stdout,"alt is first Si hit, ptr is second Si hit, alt2 is Ge\n");}
                      alt->subsys = SUBSYS_DCOMPTONB;                 // First Si hit becomes DCOMPTONB event
                      alt->alt_ecal  = ptr->ecal;                     // Second QED pixel energy
                      alt->alt2_ecal = alt2->ecal;                    // Ge energy
                      alt->esum = ptr->ecal + alt->ecal + alt2->ecal; // Add this Ge energy to the pixel energy
                      alt->tof = ptr->chan;                           // Second QED DSSD number [1-6]
                      alt->alt2_chan = ptr->alt_chan;                 // Second QED pixel number [DSSD*PIXELnumber],[0-5*0-1023]
                      alt->net_id = crystal_table[alt2->chan];        // Ge crystal number [1-64]
                      alt->delta_t = alt2->ts-ptr->ts;                // Time difference covering all three Hits
                    }else{
                      // No solution
                      if(debug_triples){fprintf(stdout,"No solution for this with alt2=Ge: %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d, [%.1f %.1f], [%.1f %.1f]\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13,angle1,(alt2->ecal+alt->ecal),angle2,(alt2->ecal+ptr->ecal));}
                    }
                  }
                }

              }
            }

            // SubSystem-specific pre-processing
            //  fprintf(stdout,"Triples presort %d: %d %d %d, %.1f %.1f %.1f (%.1f), %d %d\n",frag_idx,ptr->subsys,alt->subsys,alt2->subsys,ptr->ecal,alt->ecal,alt2->ecal,(ptr->ecal+alt->ecal+alt2->ecal),dt,dt13);
            if(DEBUG_OUTPUT){ fprintf(stdout,"\nTriples: %d: %d %d %ld | %.1f %.1f %.1f | %d vs %d %d %ld | %.1f %.1f %.1f | %d, dt=%d, sumE=%.1f",frag_idx,ptr->subsys,ptr->chan,ptr->ts,ptr->ecal,ptr->alt_ecal,ptr->esum,ptr->net_id,alt->subsys,alt->chan,alt->ts,alt->ecal,alt->alt_ecal,alt->esum,alt->net_id,dt,(ptr->ecal+alt->ecal)); }

            triples_multiplicity[(int)(alt2->subsys/22)]--; // Remove this event from consideration
          }// end of while
          triples_multiplicity[(int)(alt->subsys/22)]--; // Remove this event from consideration
        }// end of while


        if(DEBUG_OUTPUT){ fprintf(stdout,"\nEnd of pre_sort_triples\n"); }
        return(0);
      }


/*
      // Presort - The final presort to construct weighting factors based on event mixing
      //  - frag_idx is about to leave coinc window (which ends at end_idx)
      //    all other events are later than frag_idx
      int pre_sort_qed_weights(int frag_idx, int end_idx)
      {
        Grif_event *tmp, *alt2, *alt, *ptr = &grif_event[frag_idx];
        int i, j, dt, dt13, tof;
        float q1,integ2,q12,k1,k2,k12,e1,e2,e12,m,c;
        int chan1,chan2,chan3,found,pos1,pos2,qed1,qed2, ptr_swap=0, alt_swap=0;
        int clover, ge1,ge2, c1,c2,bin, p_strip, n_strip, omega, theta1,theta2,azimuthal,azimuthal2,energy_derived_theta1,energy_derived_theta2;
        float energy,ecal,correction,angle;

        // Only interested in COMPTON coincidences in this presort
        if(ptr->subsys != SUBSYS_COMPTON){ return(0); }

        // Assign chan local variable and check it is a valid channel number
        chan1 = ptr->chan;
        if( (unsigned int)chan1 >= (unsigned int)odb_daqsize ){
          fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
          return(-1);
        }
        i = frag_idx;
        while( i != end_idx ){ // need at least two events in window
          if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
          chan2 = alt->chan;
          if( (unsigned int)chan2 >= (unsigned int)odb_daqsize ){
            fprintf(stderr,"presort error: ignored event in chan2:%d\n",chan2 );
            continue;
          }
          // Only interested in COMPTON coincidences in this presort
          if(alt->subsys != SUBSYS_COMPTON){ continue; }
          //  fprintf(stdout,"Accepted Doubles: %d %d: %d: %d %d, %.1f %.1f (%.1f)\n",frag_idx,i,(ptr->ts - alt->ts),ptr->subsys,alt->subsys,ptr->ecal,alt->ecal,(ptr->ecal+alt->ecal));

          // Determine absolute time difference between timestamps
          // As ptr is about to leave the window, this will always be positive here
          // The value of presort_window_width sets the maximum (19.4us). Here set minimum of 1.1us to ensure event mixing
          if((dt = alt->ts - ptr->ts) > 40){ // > 400ns
            //  fprintf(stdout,"Scanning Doubles: %d %d: %d: %d %d, %.1f %.1f (%.1f)\n",frag_idx,i,(ptr->ts - alt->ts),ptr->subsys,alt->subsys,ptr->ecal,alt->ecal,(ptr->ecal+alt->ecal));


            // QED COMPTON EVENTS
            // Identified as subsys==SUBSYS_COMPTON
            // pos is HPGE crystal number
            // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL
            // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
            //  pos  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
            //  c2 = (ptr->alt_chan%1024);        // Pixel number [0-1023]
            //  c1 = ptr->net_id; // HPGe crystal number
            pos1 = crystal_table[chan1];
            pos2 = crystal_table[chan2];
            if(pos1 != pos2){
              qed1 = (ptr->alt_chan&1023);
              ge1 = ptr->net_id;
              qed2 = (alt->alt_chan&1023);
              ge2 = alt->net_id;

              if(ge1 != ge2){
                omega = (int)angular_diff_QEDQED(pos1, qed1, pos2, qed2);
                qed_dcs_omega_dt_TRWF->Fill(qed_dcs_omega_dt_TRWF, (int)(dt+2048), omega, 1);
                if(omega>170){ // Require back-to-back coincidence
                  // Calculate the angles
                  theta1 = (int)scattering_angle_QEDGe(pos1, qed1, ge1);
                  theta2 = (int)scattering_angle_QEDGe(pos2, qed2, ge2);
                  azimuthal = (int)azimuthal_DCS(pos1, qed1, ge1, pos2, qed2, ge2);
                  azimuthal2 = (int)energy_corrected_azimuthal_DCS(pos1, qed1, ge1, ptr->alt_ecal,pos2, qed2, ge2, alt->alt_ecal);
                  energy_derived_theta1 = (int)compton_angle(ptr->alt_ecal, 511.0);
                  energy_derived_theta2 = (int)compton_angle(alt->alt_ecal, 511.0);

                  // Time-random for weighting factors
                  // Build weighting factors here
                  qed_dcs_azi_TRWF_t->Fill(qed_dcs_azi_TRWF_t, azimuthal, 1);
                  // Scattering angle 70 to 110
                  if(theta1>69 && theta1<111 && theta2>69 && theta2<111){
                    qed_dcs_azi_TRWF->Fill(qed_dcs_azi_TRWF, azimuthal, 1);
                    // Scattering angle 93 to 103
                    if(theta1>92 && theta1<104 && theta2>92 && theta2<104){
                      qed_dcs_azi_TRWF_tg->Fill(qed_dcs_azi_TRWF_tg, azimuthal, 1);
                    }
                  }
                  // I know this code is ugly. Sorry.
                  if(energy_derived_theta1>=0 && energy_derived_theta1<=180 && energy_derived_theta2>=0 && energy_derived_theta2<=180){
                    qed_dcs_azi_TRWF_bins1->Fill(qed_dcs_azi_TRWF_bins1, azimuthal2, 1);

                    if(energy_derived_theta1>=10 && energy_derived_theta1<=170 && energy_derived_theta2>=10 && energy_derived_theta2<=170){
                      qed_dcs_azi_TRWF_bins2->Fill(qed_dcs_azi_TRWF_bins2, azimuthal2, 1);

                      if(energy_derived_theta1>=20 && energy_derived_theta1<=160 && energy_derived_theta2>=20 && energy_derived_theta2<=160){
                        qed_dcs_azi_TRWF_bins3->Fill(qed_dcs_azi_TRWF_bins3, azimuthal2, 1);

                        if(energy_derived_theta1>=30 && energy_derived_theta1<=150 && energy_derived_theta2>=30 && energy_derived_theta2<=150){
                          qed_dcs_azi_TRWF_bins4->Fill(qed_dcs_azi_TRWF_bins4, azimuthal2, 1);

                          if(energy_derived_theta1>=40 && energy_derived_theta1<=140 && energy_derived_theta2>=40 && energy_derived_theta2<=140){
                            qed_dcs_azi_TRWF_bins5->Fill(qed_dcs_azi_TRWF_bins5, azimuthal2, 1);

                            if(energy_derived_theta1>=50 && energy_derived_theta1<=130 && energy_derived_theta2>=50 && energy_derived_theta2<=130){
                              qed_dcs_azi_TRWF_bins6->Fill(qed_dcs_azi_TRWF_bins6, azimuthal2, 1);

                              if(energy_derived_theta1>=60 && energy_derived_theta1<=120 && energy_derived_theta2>=60 && energy_derived_theta2<=120){
                                qed_dcs_azi_TRWF_bins7->Fill(qed_dcs_azi_TRWF_bins7, azimuthal2, 1);

                                if(energy_derived_theta1>=70 && energy_derived_theta1<=110 && energy_derived_theta2>=70 && energy_derived_theta2<=110){
                                  qed_dcs_azi_TRWF_bins8->Fill(qed_dcs_azi_TRWF_bins8, azimuthal2, 1);

                                  if(energy_derived_theta1>=93 && energy_derived_theta1<=103 && energy_derived_theta2>=93 && energy_derived_theta2<=103){
                                    qed_dcs_azi_TRWF_bins8a->Fill(qed_dcs_azi_TRWF_bins8a, azimuthal2, 1);
                                  }

                                  if(energy_derived_theta1>=80 && energy_derived_theta1<=100 && energy_derived_theta2>=80 && energy_derived_theta2<=100){
                                    qed_dcs_azi_TRWF_bins9->Fill(qed_dcs_azi_TRWF_bins9, azimuthal2, 1);

                                    if(energy_derived_theta1>=85 && energy_derived_theta1<=95 && energy_derived_theta2>=85 && energy_derived_theta2<=95){
                                      qed_dcs_azi_TRWF_bins10->Fill(qed_dcs_azi_TRWF_bins10, azimuthal2, 1);
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                } // end of omega condition
              }
            }
          } // End of dt condition

        }// end of while


        if(DEBUG_OUTPUT){ fprintf(stdout,"\nEnd of pre_sort_qed_weights\n"); }
        return(0);
      }
      */

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
