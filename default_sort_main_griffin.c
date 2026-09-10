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
        //printf("=================End of Event================\n");

        return(0);
      }

            //#######################################################################
            //########        Individual channel singles HISTOGRAMS        ##########
            //#######################################################################

      int fill_chan_histos(Grif_event *ptr)
      {
        static int event;
        int chan, sys, pos;

        // Check for invalid channel numbers, prossibly due to data corruption
        chan = ptr->chan;
        if( (unsigned int)chan >= (unsigned int)odb_daqsize ){
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
          if( (unsigned int)sys < MAX_SUBSYS ){
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

          int fill_singles_histos(Grif_event *ptr)
          {
            int i, j, dt, pu, nhits, chan, pos, sys, elem, clover, class, c1,c2, index, offset, bin, ecal, esum, alt_ecal, integ, psd;
            int angle, ge_corrected_angle;
            char *name, c;

            int quick_reorder_qed_strips[NUM_QED_REORDERS][32] = {
              // A B C D
              {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},
              {1, 0, 3, 2, 5, 4, 7, 6, 9, 8,11,10,13,12,15,14,17,16,19,18,21,20,23,22,25,24,27,26,29,28,31,30},

              // A B C D, A reversed
              {7, 6, 5, 4, 3, 2, 1, 0,   8, 9,10,11,12,13,14,15,  16,17,18,19,20,21,22,23,  24,25,26,27,28,29,30,31},
              {6, 7, 4, 5, 2, 1, 0, 1,   9, 8,11,10,13,12,15,14,  17,16,19,18,21,20,23,22,  25,24,27,26,29,28,31,30},

              // A B C D, B reversed
              {0, 1, 2, 3, 4, 5, 6, 7,  15,14,13,12,11,10, 9, 8,  16,17,18,19,20,21,22,23,  24,25,26,27,28,29,30,31},
              {1, 0, 3, 2, 5, 4, 7, 6,  14,15,12,13,10,11, 8, 9,  17,16,19,18,21,20,23,22,  25,24,27,26,29,28,31,30},

              // A B C D, C reversed
              {0, 1, 2, 3, 4, 5, 6, 7,   8, 9,10,11,12,13,14,15,  23,22,21,20,19,18,17,16,  24,25,26,27,28,29,30,31},
              {1, 0, 3, 2, 5, 4, 7, 6,   9, 8,11,10,13,12,15,14,  22,23,20,21,18,19,16,17,  25,24,27,26,29,28,31,30},

              // A B C D
              {0, 1, 2, 3, 4, 5, 6, 7,   8, 9,10,11,12,13,14,15,  16,17,18,19,20,21,22,23,  31,30,29,28,27,26,25,24},
              {1, 0, 3, 2, 5, 4, 7, 6,   9, 8,11,10,13,12,15,14,  17,16,19,18,21,20,23,22,  30,31,28,29,26,27,24,25},

            };

            chan = ptr->chan;
            // Check for invalid channel numbers, prossibly due to data corruption
            if( (unsigned int)chan >= (unsigned int)odb_daqsize ){
              fprintf(stderr,"Invalid channel number in fill_singles_histos(), %d\n",chan);
              return(-1);
            }
            sys = ptr->subsys;
            // Check this is a valid susbsytem type
            if( (unsigned int)sys > MAX_SUBSYS ){
              return(-1);
            }
            mult_hist[sys]->Fill(mult_hist[sys], ptr->multiplicity, 1); // Fill multiplicity histograms
            // Get the position for this fragment
            pos  = crystal_table[chan];
            ecal = (int)ptr->ecal;

            switch (sys){
              case SUBSYS_HPGE_A: // GRGa
              if( (unsigned int)pos < 64 ){
                ge_sum->Fill(ge_sum, ecal, 1);
                ge_xtal->Fill(ge_xtal, pos, ecal, 1);


                index = (pos<32);
                //  ge_sum_hem[index]->Fill(ge_sum_hem[index], ecal, 1);
                //    ge_sum_hem[0]->Fill(ge_sum_hem[0], ecal, 1);
                //  ge_sum_hem[1]->Fill(ge_sum_hem[1], ecal, 1);
                ge_sum_us->Fill(ge_sum_us, ecal, 1);
                //  TH1I_Fill(ge_sum_hem[0], ecal, 1);
                /*
                // Separate sum spectra for upstream and downstream
                if(pos<32){
                ge_sum_ds->Fill(ge_sum_ds, ecal, 1);
              }else{
              ge_sum_us->Fill(ge_sum_us, ecal, 1);
            }
            */

            /*
            // Use Function Pointer Array instead of conditional logic for execution speed.
            // The outcome is 50-50 so we want to avoid CPU branch misprediction and pipeline stalls
            typedef void (*func_ptr)(void);
            static func_ptr branch_table[2] = {TH1I_Fill(ge_sum_us, ecal, 1),TH1I_Fill(ge_sum_ds, ecal, 1)};
            branch_table[pos<32]();
            */

            // Beta-gated HPGe singles spectra
            if(ptr->tof>0){
              ge_sum_b->Fill(ge_sum_b, ecal, 1);       // beta-gated Ge sum energy spectrum
              ge_sum_b_ab->Fill(ge_sum_b_ab, (int)ptr->esum, 1); // beta-gated Ge addback spectrum
            }
            if(ptr->tof & 1){ ge_sum_b_sep->Fill(ge_sum_b_sep, ecal, 1); } // Sceptar-gated Ge sum energy spectrum
            if(ptr->tof & 2){ ge_sum_b_zds->Fill(ge_sum_b_zds, ecal, 1); } // Zds-gated Ge sum energy spectrum
            if(ptr->tof & 4){ ge_sum_b_art->Fill(ge_sum_b_art, ecal, 1); } // Aries-gated Ge sum energy spectrum
            if(ptr->tof & 8){ ge_sum_b_artT->Fill(ge_sum_b_artT, ecal, 1); } // Aries-gated Ge sum energy spectrum
            if(ptr->tof & 16){ ge_sum_b_artR->Fill(ge_sum_b_artR, ecal, 1); } // Aries-gated Ge sum energy spectrum
            if(ptr->tof & 32){ ge_sum_b_artS->Fill(ge_sum_b_artS, ecal, 1); } // Aries-gated Ge sum energy spectrum

            // Pile-up
            pu = ptr->pileup;
            ge_pu_type->Fill(ge_pu_type, pu, 1);
            nhits = ptr->nhit;
            ge_nhits_type->Fill(ge_nhits_type, nhits, 1);

            // The PU class is assigned in the presort, use ptr->pu_class for pileup class
            class = ptr->pu_class;
            integ = ptr->integ1;
            ge_pu_class->Fill(ge_pu_class, class, 1);  // pileup class value
            ge_sum_class[class]->Fill(ge_sum_class[class], ecal, 1);  // energy spectrum of pileup class value
            ge_e_vs_k_class[class]->Fill(ge_e_vs_k_class[class], ecal, integ, 1);  // energy spectrum of pileup class value

            // Fill pileup histograms
            switch(class){
              case PU_SINGLE_HIT:  // single hit
              ge_1hit[pos]->Fill(ge_1hit[pos], ecal, 1);
              ge_xtal_1hit->Fill(ge_xtal_1hit, pos, ecal, 1);
              break;
              case PU_3HIT_1ST: // 3-hit pileup
              case PU_3HIT_2ND:
              case PU_3HIT_3RD:
              ge_3hit[pos]->Fill(ge_3hit[pos], ecal, 1);
              ge_xtal_3hit->Fill(ge_xtal_3hit, pos, ecal, 1);
              break;
              case PU_2HIT_A1ST: // Select first Hit of two pileup events
              case PU_2HIT_B1ST:
              case PU_2HIT_C1ST:
              ge_2hit[pos]->Fill(ge_2hit[pos], ecal, 1); // 2-hit pileup
              ge_xtal_2hit->Fill(ge_xtal_2hit, pos, ecal, 1);
              // The following used for mapping the k2 dependant correction
              ge_e_vs_k_2hit_first[pos]->Fill(ge_e_vs_k_2hit_first[pos], ecal, integ, 1);  // energy1 vs k1 spectrum of 1st Hit of 2hit pileup events
              break;
              case PU_2HIT_A2ND: // Select second Hit of two pileup events
              case PU_2HIT_B2ND:
              case PU_2HIT_C2ND:
              ge_2hit[pos]->Fill(ge_2hit[pos], ecal, 1); // 2-hit pileup
              ge_xtal_2hit->Fill(ge_xtal_2hit, pos, ecal, 1);
              // The following is not used for mapping the corrections but is a useful diagnostic
              ge_e_vs_k_2hit_second[pos]->Fill(ge_e_vs_k_2hit_second[pos], ecal, integ, 1);  // energy2 vs k2 spectrum of 2nd Hit of 2hit pileup events
              alt_ecal = ptr->alt_ecal;
              // The following 1408keV matrix used for mapping the E1 offset correction
              if(alt_ecal > 1380 && alt_ecal < 1420){ // Require 152Eu 1408keV as E1 for mapping the E1 offset
                ge_PU2_e2_v_k_gated1408[pos]->Fill(ge_PU2_e2_v_k_gated1408[pos], ecal, integ, 1);  // e2 vs k2 for fixed e1
                // The following x-rays matrix used for mapping the k2 dependant correction for E2
                // E2 vs k2 gated on fixed x-ray energies
                // The E2 energy has 1272keV subtracted from it to put the 1408keV peak around 136keV to allow a smaller matrix side and easier processing in the app
                if(alt_ecal > 5 && alt_ecal < 130){ // Require 152Eu x rays (or 121keV because x rays are attenuated in some channels) as E1 for mapping the k2 dependance
                  ge_PU2_e2_v_k_gatedxrays[pos]->Fill(ge_PU2_e2_v_k_gatedxrays[pos], (ecal - 1272), integ, 1);  // e2 vs k2 for fixed e1
                }
              }
              break;
              default: break;
            }

            clover = (int)(pos>>2)+1;
            if( clover >= 0 && clover < N_CLOVER && esum >= 0 ){   // ge addback
              esum = (int)ptr->esum;
              ge_ab_e[clover]->Fill(ge_ab_e[clover], esum, 1);
              ge_sum_ab   ->Fill(ge_sum_ab, esum, 1);
              //  if(ptr->suppress != 1){ ge_sum_ab_sup->Fill(ge_sum_ab_sup,(int)ptr->esum, 1); } // Addback and BGO suppressed
              //  else{  ge_sum_ab_sup_rej->Fill(ge_sum_ab_sup_rej,(int)ptr->esum, 1); }          // Addback and What is rejected by BGO suppressed
              if(ptr->suppress == 1){ ge_sum_ab_sup_rej->Fill(ge_sum_ab_sup_rej,esum, 1); } // Addback and What is rejected by BGO suppressed
              else{
                ge_sum_ab_sup->Fill(ge_sum_ab_sup,esum, 1);               // Addback and BGO suppressed
                ge_ab_sup_e[clover]->Fill(ge_ab_sup_e[clover],esum, 1); // Addback and BGO suppressed per clover
              }

              if( clover < 9 ){ // Separate Addback sum for upstream and downstream
                ge_sum_ab_us->Fill(ge_sum_ab_us, ecal, 1);
              } else {
                ge_sum_ab_ds->Fill(ge_sum_ab_ds, ecal, 1);
              }
            }

            // PPG Cycles histograms
            if(ppg_cycles_active==1){
              ge_cycle_code[ppg_current_pattern]->Fill(ge_cycle_code[ppg_current_pattern], ecal, 1);
              bin = (int)((ptr->ts-ppg_cycle_start)/ppg_cycles_binning_factor);  // convert 10ns to binning size set as Global
              ge_cycle_activity->Fill(ge_cycle_activity, bin, 1);
              ge_e_vs_cycle_time->Fill(ge_e_vs_cycle_time, bin, ecal, 1);
              if(ppg_cycle_number<MAX_CYCLES){
                gea_cycle_num[ppg_cycle_number]->Fill(gea_cycle_num[ppg_cycle_number], bin, 1);
                cycle_num_vs_ge->Fill(cycle_num_vs_ge, ppg_cycle_number, bin, 1);
                cycle_num_vs_geEnergy[pos]->Fill(cycle_num_vs_geEnergy[pos], ppg_cycle_number, ecal, 1);
                if(class == PU_SINGLE_HIT){
                  gea_cycle_num_sh[ppg_cycle_number]->Fill(gea_cycle_num_sh[ppg_cycle_number], bin, 1);
                  cycle_num_vs_sh->Fill(cycle_num_vs_sh, ppg_cycle_number, bin, 1);
                }else{
                  gea_cycle_num_pu[ppg_cycle_number]->Fill(gea_cycle_num_pu[ppg_cycle_number], bin, 1);
                  cycle_num_vs_pu->Fill(cycle_num_vs_pu, ppg_cycle_number, bin, 1);
                }
                if((ecal>=ppg_cycles_gamma_gate_min) && (ecal<=ppg_cycles_gamma_gate_max)){
                  gea_cycle_num_g[ppg_cycle_number]->Fill(gea_cycle_num_g[ppg_cycle_number], bin, 1);
                  if(class == PU_SINGLE_HIT){
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
            geb_xtal->Fill(geb_xtal, pos, ecal, 1);
            // PPG Cycles histograms
            if(ppg_cycles_active==1){
              bin = (int)((ptr->ts-ppg_cycle_start)/ppg_cycles_binning_factor);  // convert 10ns to binning size set as Global
              if(ppg_cycle_number<MAX_CYCLES){
                geb_cycle_num[ppg_cycle_number]->Fill(geb_cycle_num[ppg_cycle_number], bin, 1);
                cycle_num_vs_ge_b->Fill(cycle_num_vs_ge_b, ppg_cycle_number, bin, 1);
                class = ptr->pu_class;
                if(class == PU_SINGLE_HIT){
                  geb_cycle_num_sh[ppg_cycle_number]->Fill(geb_cycle_num_sh[ppg_cycle_number], bin, 1);
                  cycle_num_vs_sh_b->Fill(cycle_num_vs_sh_b, ppg_cycle_number, bin, 1);
                }else{
                  geb_cycle_num_pu[ppg_cycle_number]->Fill(geb_cycle_num_pu[ppg_cycle_number], bin, 1);
                  cycle_num_vs_pu_b->Fill(cycle_num_vs_pu_b, ppg_cycle_number, bin, 1);
                }
                if((ecal>=ppg_cycles_gamma_gate_min) && (ecal<=ppg_cycles_gamma_gate_max)){
                  geb_cycle_num_g[ppg_cycle_number]->Fill(geb_cycle_num_g[ppg_cycle_number], bin, 1);
                  if(class == PU_SINGLE_HIT){
                    geb_cycle_num_sh_g[ppg_cycle_number]->Fill(geb_cycle_num_sh_g[ppg_cycle_number], bin, 1);
                    cycle_num_vs_ge_b_sh_g->Fill(cycle_num_vs_ge_b_sh_g, ppg_cycle_number, bin, 1);
                  }
                }
              }
            }
          } break;
          case SUBSYS_BGO: // BGOs
          pos  = crystal_table[chan];
          elem = element_table[chan];
          if( pos < 0 || pos > 63 ){
            fprintf(stderr,"bad bgo crystal[%d] for chan %d\n", pos, chan);
          } else if( elem < 1 || elem > 5 ){
            fprintf(stderr,"bad bgo element[%d] for chan %d, %s, subsys %d\n", elem, chan, chan_name[chan], sys);
          } else {
            pos *= 5; pos += (elem-1);
            bgo_xtal->Fill(bgo_xtal, pos, ecal, 1);
            if(elem <3){ // front
              pos  = crystal_table[chan];
              pos *= 2; pos += (elem-1);
              bgof_xtal->Fill(bgof_xtal, pos, ecal, 1);
            } else if(elem>4){ // back
              pos  = crystal_table[chan];
              bgob_xtal->Fill(bgob_xtal, pos, ecal, 1);
            } else{ // side
              pos  = crystal_table[chan];
              pos *= 2; pos += (elem-3);
              bgos_xtal->Fill(bgos_xtal, pos, ecal, 1);
            }
          }  break;
          case SUBSYS_LABR_BGO: // Ancillary BGOs
          pos  = crystal_table[chan];
          elem = element_table[chan];
          if( pos < 1 || pos > 8 ){
            fprintf(stderr,"bad ancillary bgo crystal[%d] for chan %d\n", pos, chan);
          } else if( elem < 1 || elem > 3 ){
            fprintf(stderr,"bad ancillary bgo element[%d] for chan %d\n", elem, chan);
          } else {
            pos *= 3; pos += (elem-1);
            bgoa_xtal->Fill(bgoa_xtal, pos, ecal, 1);
          } break;
          case SUBSYS_PACES: // PACES
          paces_sum->Fill(paces_sum, ecal, 1);
          if(ptr->tof>0){ paces_sum_b->Fill(paces_sum_b, ecal, 1); }      // beta-gated PACES sum energy spectrum
          pos  = crystal_table[chan];
          if( pos < 1 || pos > 5 ){
            fprintf(stderr,"bad PACES crystal[%d] for chan %d\n", pos, chan);
          } else {
            paces_xtal->Fill(paces_xtal, pos, ecal, 1);
          } break;
          case SUBSYS_LABR_L: // LaBr3 (LBL)
          labr_sum->Fill(labr_sum, ecal, 1);
          pos  = crystal_table[chan];
          if( pos < 1 || pos > 8 ){
            fprintf(stderr,"bad LaBr3 crystal[%d] for chan %d\n", pos, chan);
          } else {
            labr_xtal->Fill(labr_xtal, pos, ecal, 1);
          } break;
          case SUBSYS_SCEPTAR:
          sceptar_xtal->Fill(sceptar_xtal, crystal_table[chan], ecal, 1);
          break;
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
                if(index>=0 && index<(int)((N_LABR)*((N_LABR-1)>>1))+2){
                  offset = tac_lbl_combo_offset[index];
                  tac_labr_hist[index]->Fill(tac_labr_hist[index], ecal+offset, 1);
                  tac_labr_hist_uncal[index]->Fill(tac_labr_hist_uncal[index], (int)(( ptr->integ1 == 0 ) ? ptr->q1 : spread(ptr->q1)/ptr->integ1), 1);
                }
                // Calibrated TAC spectra
                final_tac[crystal_table[chan]-1]->Fill(final_tac[crystal_table[chan]-1], ecal+offset, 1);
                final_tac_sum->Fill(final_tac_sum, ecal+offset, 1);

                // A 3d histogram of first LBL energy vs second LBL energy vs TAC
                // lbl_lbl_tac->Fill(lbl_lbl_tac, (int)((ptr->e2cal/10)*(ptr->e3cal/10)), (int)(ptr->ecal)+offset, 1); // LBL energy vs LBL energy vs TAC
                if(ptr->ecal>5 && ptr->q2>5 && ptr->q3>5){
                  bin = (int)(((ptr->q2/10)*400)+(ptr->q3/10));
                  lbl_lbl_tac->Fill(lbl_lbl_tac, ecal+offset-250, bin, 1); // LBL energy vs LBL energy vs TAC
                }

                // Compton Walk matrix for calibrations
                // First LBL gated on 1332keV, this matrix is second LBL E vs TAC
                if(ecal>5 && crystal_table[chan] == 1){ // Use the First TAC (TAC01)
                  if(c1 == 0 && c2>0 && ptr->q2>1252 && ptr->q2<1412 && ptr->q3>5){ // LBL01 gated on 1332keV
                    tac_labr_CompWalk[c2]->Fill(tac_labr_CompWalk[c2], ecal+offset, (int)ptr->q3, 1); // TAC01 and other LBL energy
                  }
                }else if(ecal>5 && crystal_table[chan] == 2){
                  if(c1 == 1 && c2 == 2 && ptr->q2>1252 && ptr->q2<1412 && ptr->q3>5){ // LBL02 gated on 1332keV
                    tac_labr_CompWalk0->Fill(tac_labr_CompWalk0, ecal+offset, (int)ptr->q3, 1); // TAC02 to check LBL01
                  }
                }
              }
            } break;
            case SUBSYS_DESCANT: break;
            case SUBSYS_DESWALL: // DESCANT Wall
            pos  = crystal_table[chan];
            if( pos < 1 || pos > 60 ){
              fprintf(stderr,"bad descant wall detector[%d] for chan %d\n", pos, chan);
            }else{
              psd = ptr->psd;
              alt_ecal = ptr->alt_ecal;
              desw_psd[pos]       -> Fill(desw_psd[pos],   psd,       1);
              if(ptr->tof>0){
                desw_tof[pos]       -> Fill(desw_tof[pos],   (int)ptr->tof,       1);
                desw_tof_corr[pos]  -> Fill(desw_tof_corr[pos],   alt_ecal,       1);
                if(psd>10 && psd<710){
                  desw_tof_psd[pos]  -> Fill(desw_tof_psd[pos],   alt_ecal,       1);
                }
              }

              desw_sum_e->Fill(desw_sum_e, ecal, 1);
              if(psd>0){ desw_sum_psd->Fill(desw_sum_psd, psd, 1); }
              if(alt_ecal>0){ desw_sum_tof->Fill(desw_sum_tof, alt_ecal, 1); } // alt_ecal = corrected time-of-flight
              if(ecal>5){
                desw_e_xtal->Fill(desw_e_xtal, pos, ecal, 1);
                if(psd>5){
                  desw_psd_e->Fill(desw_psd_e, psd, ecal, 1);
                  desw_psd_q->Fill(desw_psd_q, psd, (int)ptr->q1, 1);
                  desw_psd_cc->Fill(desw_psd_cc, psd, (int)ptr->cc_short, 1);
                  desw_q_cc->Fill(desw_q_cc, (int)ptr->q1, (int)ptr->cc_short, 1);
                }
              }
              if(alt_ecal>5){ // DESCANT Wall ptr->alt_ecal=corrected-TOF (ptr->tof is TOF)
                desw_tof_xtal->Fill(desw_tof_xtal, pos, alt_ecal, 1);
                if(psd>5){ desw_psd_tof->Fill(desw_psd_tof, psd, alt_ecal, 1); }
              }
            }
            break;
            case SUBSYS_ARIES_A: // ARIES Standard Output
            aries_sum->Fill(aries_sum, ecal, 1);
            pos  = crystal_table[chan];
            if( pos < 1 || pos > 76 ){
              fprintf(stderr,"bad aries tile[%d] for chan %d\n", pos, chan);
            } else {
              aries_xtal->Fill(aries_xtal, pos, ecal, 1);
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
            rcmp_sum->Fill(rcmp_sum, ecal, 1);
            if(esum>0){
              rcmp_fb_sum->Fill(rcmp_fb_sum, esum, 1);
            }
            pos  = crystal_table[chan];
            elem = (int)(element_table[chan] + (int)(polarity_table[chan]*N_RCMP_STRIPS)); // polarity_table value is 0 or 1
            if( pos < 1 || pos > 6 ){
              fprintf(stderr,"bad RCMP DSSD[%d] for chan %d, elem%d, pol%d\n", pos, chan, elem, polarity_table[chan]);
            } else if( elem < 0 || elem > 63 ){
              fprintf(stderr,"bad RCMP strip[%d] for chan %d, pos%d, pol%d\n", elem, chan, pos, polarity_table[chan]);
            } else {
              rcmp_strips[(pos-1)]->Fill(rcmp_strips[(pos-1)], elem, ecal, 1);
            }
            break;
            case SUBSYS_QED_PIXEL: // QED pixel is a coincidence between a front and back strip
            qed_sum->Fill(qed_sum, ecal, 1);     // p strips
            qed_sum->Fill(qed_sum, (int)ptr->alt_ecal, 1); // n strips
            qed_fb_sum->Fill(qed_fb_sum, ecal, 1);     // fb coincidence
            pos  = crystal_table[ptr->chan]-1; // QED DSSD number [1-6]
            elem = ptr->alt_chan; // QED pixel number [0-1023]
            qed_fb[pos]->Fill(qed_fb[pos], ecal, (int)ptr->alt_ecal, 1); // front-back energy
            //  qed_psd_e[pos]->Fill(qed_psd_e[pos], ecal, ptr->psd, 1); // qed psd
            if(ecal>QED_STRIP_THRESHOLD){
              qed_hit[pos]->Fill(qed_hit[pos], (int)(elem/N_QED_STRIPS), (elem%N_QED_STRIPS), 1); // QED DSSD hitpattern
              qed_strips[pos]->Fill(qed_strips[pos], (int)(elem/N_QED_STRIPS), ecal, 1); // p strip energies
              qed_strips[pos]->Fill(qed_strips[pos], (elem%N_QED_STRIPS)+N_QED_STRIPS, (int)ptr->alt_ecal, 1); // n strip energies

              /*
              // Channel mapping, strip reorder trials.
              for(i=0; i<NUM_QED_REORDERS; i++){
              for(j=0; j<NUM_QED_REORDERS; j++){
              qed_hit_trials[pos]->Fill(qed_hit_trials[pos], (int)(quick_reorder_qed_strips[i][elem/N_QED_STRIPS])+(i*64), (quick_reorder_qed_strips[j][elem%N_QED_STRIPS])+(j*64), 1); // QED DSSD hitpattern
            }
          }
          */

          // PPG Cycles histograms
          if(ppg_cycles_active==1){
            if(ppg_cycle_number<MAX_CYCLES){
              if((int)(elem/N_QED_STRIPS) == 10){
                cycle_num_vs_qedEnergy[pos]->Fill(cycle_num_vs_qedEnergy[pos], ppg_cycle_number, ecal, 1);
              }
            }
          }

        }
        break;
        case SUBSYS_COMPTON: // COMPTON is a coincidence between a DSSD pixel and a HPGE with sum energy of 511keV
        // QED COMPTON EVENTS
        // Identified as subsys==SUBSYS_COMPTON
        // pos is HPGE crystal number
        // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL, esum is to total energy
        // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
        pos  = crystal_table[chan]; // QED DSSD number [1-6]
        c2 = (ptr->alt_chan&1023);        // Pixel number [0-1023]
        c1 = ptr->net_id; // HPGe crystal number
        angle = (int)scattering_angle_QEDGe(pos,c2,c1);
        alt_ecal = (int)ptr->alt_ecal;

        if(DEBUG_OUTPUT){ fprintf(stdout,"\nfill_singles_histos(): COMPTON FILL HISTOS: %d %d %ld | %.1f %.1f %.1f | %d dt=%d\n",ptr->subsys,ptr->chan,ptr->ts,ptr->ecal,ptr->alt_ecal,ptr->esum,ptr->net_id,ptr->delta_t); }
        //fprintf(stdout,"\nfill_singles_histos(): COMPTON FILL HISTOS: %d %d %ld | %.1f %.1f %.1f | %d dt=%d | %d %d %d %d\n",ptr->subsys,ptr->chan,ptr->ts,ptr->ecal,ptr->alt_ecal,ptr->esum,ptr->net_id,ptr->delta_t,(int)(scattering_angle_QEDGe(pos,c2,c1)),(int)(scattering_angle_GeQED(pos,c2,c1)),(int)compton_angle(ptr->alt_ecal,511.0),(int)compton_angle(ptr->ecal,511.0));

        qedE_ge_dt_c->Fill(qedE_ge_dt_c, ptr->delta_t+512, ecal, 1);
        qed_geE_dt_c->Fill(qed_geE_dt_c, ptr->delta_t+512, alt_ecal, 1);
        qed_theta_dt_c->Fill(qed_theta_dt_c, ptr->delta_t+512, angle, 1);
        ge_qed_c->Fill(ge_qed_c, alt_ecal, ecal, 1);
        qedE_ge_theta_sum_c->Fill(qedE_ge_theta_sum_c, ecal, angle, 1);
        qed_geE_theta_sum_c->Fill(qed_geE_theta_sum_c, alt_ecal, angle, 1);
        qed_theta->Fill(qed_theta, angle, 1);

        if((angle>=compton_angle(ptr->alt_ecal,QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && (angle<=compton_angle(ptr->alt_ecal,QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
          qedE_ge_theta_sum_c_s->Fill(qedE_ge_theta_sum_c_s, ecal, angle, 1);
          qed_geE_theta_sum_c_s->Fill(qed_geE_theta_sum_c_s, alt_ecal, angle, 1);

        }else if(((int)(scattering_angle_GeQED(pos,c2,c1))>=compton_angle(ecal,QED_GAMMA_ENERGY)-QED_ANGLE_WINDOW) && ((int)(scattering_angle_GeQED(pos,c2,c1))<=compton_angle(ecal,QED_GAMMA_ENERGY)+QED_ANGLE_WINDOW)){
          // Look for Ge-Si backwards scattering where Ge was first
          qedE_ge_theta_sum_c_g->Fill(qedE_ge_theta_sum_c_g, ecal, angle, 1);
          qed_geE_theta_sum_c_g->Fill(qed_geE_theta_sum_c_g, alt_ecal, angle, 1);
        }
        break;

        case SUBSYS_DCOMPTONA: // DCOMPTONA is a coincidence between a DSSD pixel and a HPGE addback with sum energy of 511keV
        // QED DCOMPTONA EVENTS
        // Ge with addback is a Double Compton scatter (DSSD-Ge-Ge)
        // Identified as subsys==SUBSYS_DCOMPTONA
        // In DCOMPTONA event, ecal will be QED_PIXEL and alt_ecal will be Ge addback sum energy
        // In DCOMPTONA event, crystal_table[chan]=pos will be QED DSSD number [1-6], alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
        // In DCOMPTONA event, net_id will be first HPGe crystal number [1-64], alt2_chan will be second HPGe crystal number [1-64]
        pos  = crystal_table[chan]; // QED DSSD number [1-6]
        c2 = (ptr->alt_chan&1023);  // Pixel number [0-1023] (faster than alt_chan%1024)
        c1 = ptr->net_id; // HPGe crystal number
        angle = (int)scattering_angle_QEDGe(pos,c2,c1); // This is the initial theta angle in DCompton

        dcsa_theta->Fill(dcsa_theta, (int)(angle), 1);
        dcsaE_ge_theta->Fill(dcsaE_ge_theta, ecal, (int)(angle), 1);
        dcsa_geE_theta->Fill(dcsa_geE_theta, (int)ptr->alt_ecal, (int)(angle), 1);
        break;
        default: break; // Unrecognized or unprocessed dtype
      }// end of switch
      return(0);
    }

    int fill_ge_coinc_histos(Grif_event *ptr, Grif_event *alt, int abs_dt)
    {
      int c1, c2, c3, c4, pos, bin, angle_idx, coinc_ecal, scatt_esum, totalEnergy, ge_corrected_angle;
      int pos1, qed1, pos2, qed2, ptr_ecal, alt_ecal, ptr_esum, alt_esum, p_strip, n_strip;
      double angle, initial_theta, omega, azimuthal;
      switch(alt->subsys){
        case SUBSYS_HPGE_A:
        ptr_ecal = (int)ptr->ecal; alt_ecal = (int)alt->ecal;
        gg_dt->Fill(gg_dt, (int)((ptr->ts - alt->ts)+(DT_SPEC_LENGTH>>1)), ptr_ecal, 1); // This dt result is always negative
        gg_dt->Fill(gg_dt, (int)((alt->ts - ptr->ts)+(DT_SPEC_LENGTH>>1)), alt_ecal, 1); // This dt result is always positive
        if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_HPGE_A]) ){

          // PPG Cycles histograms
          if(ppg_cycles_active==1){
            gg_cycle_code[ppg_current_pattern]->Fill(gg_cycle_code[ppg_current_pattern], ptr_ecal, alt_ecal, 1);
          }

          c1 = crystal_table[ptr->chan];
          c2 = crystal_table[alt->chan];
          if( c1 >= 0 && c1 < 64 && c2 >= 0 && c2 < 64 ){
            if(ptr_ecal < 5 || alt_ecal < 5){ break; } // Both crystals must have a valid energy
            ptr_esum = (int)ptr->esum; alt_esum = (int)alt->esum;
            gg_hit->Fill(gg_hit, c1, c2, 1); // 2d crystal hitpattern

            // Individual Ge crystal energy in coincidence, used for angular correlations weighting factors
            gg_energy[c1]->Fill(gg_energy[c1], ptr_ecal, 1);
            gg_energy[c2]->Fill(gg_energy[c2], alt_ecal, 1);

            if( c2 == grif_opposite[c1] ){
              // 180 degree coinc matrix for summing corrections
              gg_opp->Fill(gg_opp, ptr_ecal, alt_ecal, 1);
              gg_ab_opp->Fill(gg_ab_opp, ptr_esum, alt_esum, 1);
            }

            // Ge-Ge angular correlations
            // Fill the appropriate angular bin spectrum
            // c1 and c2 run from 0 to 63 for ge_angles_145mm.
            angle_idx = ge_angles_110mm[c1][c2];
            gg_angcor_110[angle_idx]->Fill(gg_angcor_110[angle_idx], ptr_ecal, alt_ecal, 1);
            //fprintf(stdout,"%d %d have angular difference of: calculate %lf [%f]\n",c1,c2,angular_diff_GeGe(c1,c2,110),angular_bins_110mm[angle_idx]);
            // double atan2(double y, double x);
            angle_idx = ge_angles_145mm[c1][c2];
            gg_angcor_145[angle_idx]->Fill(gg_angcor_145[angle_idx], ptr_ecal, alt_ecal, 1);

            if( ptr_esum > 5 || alt_esum > 5 ){ // addback energies
              gg_ab->Fill(gg_ab, ptr_esum, alt_esum, 1);
              // Fill Compton polarimetry matrices Here
              // c1 and c2 are the same as for angular correlations. These define the plane.
              // c3 is another crystal in one of the clovers.
              // The azimuthal scattering angle is between c3 and the plane of polarization.
              if(angle_idx>11 && angle_idx<41){ // Only look at Compton Polarimetry for clovers at 90 degrees to each other. (Always checked with 145mm)
                c3 = ptr->alt_chan;
                if(c3<0){ // c1 is coincident, c2+c3 are the scattering event
                  c3 = alt->alt_chan; scatt_esum = alt_esum; coinc_ecal = ptr_esum;
                }else{ // c2 is coincident, c1+c3 are the scattering event
                  scatt_esum = alt_esum; coinc_ecal = ptr_esum;
                }
                if(c3<0 || c3>63){ break; }
                // Determine angle_idx from azimuthal.
                angle = azimuthal_GeGeGe(c1,c2,c3,110);
                angle_idx = (int)(angle / 15);
                if(angle_idx==12){ angle_idx=11; } // Put 180 degree scatters into the last valid bin
                if(angle_idx>=0 && angle_idx<N_GE_COMP_POL){
                  gg_comp_pol_110[angle_idx]->Fill(gg_comp_pol_110[angle_idx], coinc_ecal, scatt_esum, 1); // Asymmetric matrix
                }
              }
            }
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
        gb_dt->Fill(gb_dt, (int)((ptr->ts - alt->ts)+(DT_SPEC_LENGTH>>1)), (int)ptr->esum, 1);
        if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_SCEPTAR]) ){
          //  ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1); // beta-gated Ge sum energy spectrum
          //  ge_sum_b_sep->Fill(ge_sum_b_sep, (int)ptr->ecal, 1); // Sceptar-gated Ge sum energy spectrum
        }else if((ptr->ts - alt->ts)<-25){
          ge_isomer_popu->Fill(ge_isomer_popu, (int)ptr->ecal, 1); // Early gamma rays appearing earlier in time than the prompt
        }else if((ptr->ts - alt->ts)>25){
          ge_isomer_depop->Fill(ge_isomer_depop, (int)ptr->ecal, 1); // Delayed gamma rays appearing later in time than the prompt
        }
        break;
        case SUBSYS_ARIES_A:
        gb_dt->Fill(gb_dt, (int)((ptr->ts - alt->ts)+(DT_SPEC_LENGTH>>1)), (int)ptr->esum, 1);
        if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ARIES_A]) ){
          if(ptr->ecal > 10 && alt->ecal > 10){
            //  ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1);         // beta-gated Ge sum energy spectrum
            //  ge_sum_b_art->Fill(ge_sum_b_art, (int)ptr->ecal, 1); // Aries-gated Ge sum energy spectrum
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
        gb_dt->Fill(gb_dt, (int)((ptr->ts - alt->ts)+(DT_SPEC_LENGTH>>1)), (int)ptr->esum, 1);
        if( (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_ZDS_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_ZDS_A]) ){
          //  ge_sum_b->Fill(ge_sum_b, (int)ptr->ecal, 1);         // beta-gated Ge sum energy spectrum
          //  ge_sum_b_zds->Fill(ge_sum_b_zds, (int)ptr->ecal, 1); // Zds-gated Ge sum energy spectrum
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
        case SUBSYS_QED_PIXEL: // ge-QED
        c1 = crystal_table[ptr->chan]; // HPGe crystal number
        pos  = crystal_table[alt->chan]; // QED DSSD number [1-6]
        c2 = alt->alt_chan; // QED pixel number [0-1023]
        //  (int)(c2/N_QED_STRIPS); // QED p strip number [0-31]
        //  (c2%N_QED_STRIPS); // QED n strip number [0-31]
        if( c1 >= 0 && c1 < 64 && c2 >= 0 && c2 < 1024 && ptr->ecal>5 && (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL])){
          qed_p_ge_hit->Fill(qed_p_ge_hit, (int)((int)(c2/N_QED_STRIPS) + (int)((pos-1)*N_QED_STRIPS)), c1, 1);
          qed_n_ge_hit->Fill(qed_n_ge_hit, (int)((c2%N_QED_STRIPS) + (int)((pos-1)*N_QED_STRIPS)), c1, 1);

          //  fprintf(stdout,"Theta Geometric,idealized: %d,%d\n",(int)(scattering_angle_QEDGe(pos,c2,c1)),compton_angle(ptr->ecal,662.0));
          angle = (int)scattering_angle_QEDGe(pos,c2,c1);
          ptr_ecal = (int)ptr->ecal;
          alt_ecal = (int)alt->ecal;
          qedE_ge_theta_sum->Fill(qedE_ge_theta_sum, alt_ecal, angle, 1);
          qed_geE_theta_sum->Fill(qed_geE_theta_sum, ptr_ecal, angle, 1);
          qed_geE_theta_clov[(int)(c1>>2)]->Fill(qed_geE_theta_clov[(int)(c1>>2)], ptr_ecal, angle, 1);
          qed_E_theta_dssd[pos-1]->Fill(qed_E_theta_dssd[pos-1], alt_ecal, angle, 1);
          qed_geE_theta_dssd[pos-1]->Fill(qed_geE_theta_dssd[pos-1], ptr_ecal, angle, 1);
          qed_geE_theta_clov_t[(int)(c1>>2)]->Fill(qed_geE_theta_clov_t[(int)(c1>>2)], ptr_ecal, angle, 1);
          totalEnergy = ptr_ecal+alt_ecal;
          //  qed_totE_theta_sum->Fill(qed_totE_theta_sum, totalEnergy, (int)(angle), 1);
          //  qed_totE_theta[pos-1]->Fill(qed_totE_theta[pos-1], totalEnergy, (int)(angle), 1);
          qed_E_totE_sum_t->Fill(qed_E_totE_sum_t, alt_ecal, totalEnergy, 1);
          qed_geE_totE_sum_t->Fill(qed_geE_totE_sum_t, ptr_ecal, totalEnergy, 1);
          qed_theta_dt->Fill(qed_theta_dt, (int)(alt->ts-ptr->ts)+512, angle, 1);
          qed_theta_dt_cfd->Fill(qed_theta_dt_cfd, (int)(alt->cfd-ptr->cfd)+512, angle, 1);

          if(ptr_ecal>=1264 && ptr_ecal<=1284 && alt->ecal>450 && alt->ecal<650){ // Gate on 1274keV in Ge - to show 511keV peak in Si
            // Single strip vs theta needed for strip-energy calibration
            pos--; p_strip = (int)(c2/N_QED_STRIPS); n_strip = c2%N_QED_STRIPS;
            qedp_ge_theta[p_strip + pos*N_QED_STRIPS]->Fill(qedp_ge_theta[p_strip + pos*N_QED_STRIPS], alt->ecal, angle, 1);
            qedn_ge_theta[n_strip + pos*N_QED_STRIPS]->Fill(qedn_ge_theta[n_strip + pos*N_QED_STRIPS], alt->alt_ecal, angle, 1);
          }
          // 60 to 120 degree scattering angles to use the Compton scattering of 1274keV for Si strip calibration
          if(angle>100 && angle<170 && ptr_ecal>=210 && ptr_ecal<=330 && alt->ecal>850 && alt->ecal<1300){ // Gate on scattered 1274keV in Ge - to show Comptons in Si
            //if((angle>=compton_angle(ptr->ecal,QED_GAMMA2_ENERGY)-QED_ANGLE_WINDOW) && (angle<=compton_angle(ptr->ecal,QED_GAMMA2_ENERGY)+QED_ANGLE_WINDOW)){
            //  if(ptr->ecal+alt->ecal > QED_GAMMA2_ENERGY-QED_GAMMA_ENERGY_WINDOW && ptr->ecal+alt->ecal < QED_GAMMA2_ENERGY+QED_GAMMA_ENERGY_WINDOW){
            if(ptr->ecal >= secondary_energy(angle, QED_GAMMA2_ENERGY)-5 && ptr->ecal <= secondary_energy(angle, QED_GAMMA2_ENERGY)+35){
              // Single strip vs theta needed for strip-energy calibration
              pos--; p_strip = (int)(c2/N_QED_STRIPS); n_strip = c2%N_QED_STRIPS;
              qedp_ge_theta[p_strip + pos*N_QED_STRIPS]->Fill(qedp_ge_theta[p_strip + pos*N_QED_STRIPS], alt->ecal, angle, 1);
              qedn_ge_theta[n_strip + pos*N_QED_STRIPS]->Fill(qedn_ge_theta[n_strip + pos*N_QED_STRIPS], alt->alt_ecal, angle, 1);
            }
          }

          if(totalEnergy>QED_GAMMA_ENERGY-QED_GAMMA_ENERGY_WINDOW && totalEnergy<QED_GAMMA_ENERGY+QED_GAMMA_ENERGY_WINDOW){

            //    fprintf(stdout,"ge-QED coinc. [%d][%d,%d,%d,%d], dt = %ld - %ld = %ld\n",crystal_table[ptr->chan],crystal_table[alt->chan],alt->alt_chan,(int)(alt->alt_chan/N_QED_STRIPS),(int)(alt->alt_chan%N_QED_STRIPS),ptr->ts,alt->ts,ptr->ts - alt->ts);
            /*
            if(c1==13){
            fprintf(stdout,"COINCS: %d %f %d %ld %d %ld\n",c1,ptr->ecal,ptr->trig_acc,ptr->ts,alt->trig_acc,alt->ts);
            fprintf(stdout,"COINC, %d %d scattering_angle_QEDGe(%d,%d,%d) %f\n",(int)ptr->ecal,(int)alt->ecal,pos,c2,c1,angle);
          }
          */
          if(DEBUG_OUTPUT){ fprintf(stdout,"\nCOINC PIXEL: %d %d %ld | %.1f %.1f %.1f | %d vs %d %d %ld | %.1f %.1f %.1f | %d, dt=%ld, sumE=%.1f\n",ptr->subsys,ptr->chan,ptr->ts,ptr->ecal,ptr->alt_ecal,ptr->esum,ptr->net_id,alt->subsys,alt->chan,alt->ts,alt->ecal,alt->alt_ecal,alt->esum,alt->net_id,(ptr->ts-alt->ts),(ptr->ecal+alt->ecal)); }
          qedE_ge_dt->Fill(qedE_ge_dt, (int)(alt->ts-ptr->ts)+512, alt_ecal, 1);
          qed_geE_dt->Fill(qed_geE_dt, (int)(alt->ts-ptr->ts)+512, ptr_ecal, 1);
          qedE_ge_theta_sum_t->Fill(qedE_ge_theta_sum_t, alt_ecal, angle, 1);
          qed_geE_theta_sum_t->Fill(qed_geE_theta_sum_t, ptr_ecal, angle, 1);

          ge_corrected_angle = compton_angle(ptr->ecal,QED_GAMMA_ENERGY);
          qedE_ge_thetaI_sum_t->Fill(qedE_ge_thetaI_sum_t, alt_ecal, ge_corrected_angle, 1);
          qed_geE_thetaI_sum_t->Fill(qed_geE_thetaI_sum_t, ptr_ecal, ge_corrected_angle, 1);
          qed_geE_thetaDiff_sum_t->Fill(qed_geE_thetaDiff_sum_t, ptr_ecal, (int)(angle - ge_corrected_angle)+90, 1);
        }
      }
      break;
      case SUBSYS_COMPTON: // ge-COMPTON where is a coincidence between a DSSD pixel and a HPGE with sum energy of 511keV
      // QED COMPTON EVENTS
      // Identified as subsys==SUBSYS_COMPTON
      // pos is HPGE crystal number
      // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL
      // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
      //  pos  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
      //  c2 = (ptr->alt_chan%1024);        // Pixel number [0-1023]
      //  c1 = ptr->net_id; // HPGe crystal number
      pos  = crystal_table[alt->chan]; // QED DSSD number [1-6]
      c2 = (alt->alt_chan&1023);        // Pixel number [0-1023] (faster than alt_chan%1024)
      c1 = alt->net_id; // HPGe crystal number
      angle = (int)scattering_angle_QEDGe(pos,c2,c1);
      ptr_ecal = (int)ptr->ecal;
      alt_ecal = (int)alt->ecal;

      ge_comp->Fill(ge_comp, ptr_ecal, alt_ecal, 1);
      if(ptr->esum>ptr_ecal){
        geadd_comp->Fill(geadd_comp, (int)ptr->esum, (int)(alt->esum), 1);
      }

      if( c1 >= 0 && c1 < 64 && c2 >= 0 && c2 < 1024 && ptr_ecal>5 && alt_ecal>5 && (abs_dt >= time_diff_gate_min[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL]) && (abs_dt <= time_diff_gate_max[SUBSYS_HPGE_A][SUBSYS_QED_PIXEL])){
        if(alt->esum>QED_GAMMA_ENERGY-QED_GAMMA_ENERGY_WINDOW && alt->esum<QED_GAMMA_ENERGY+QED_GAMMA_ENERGY_WINDOW){

          if(DEBUG_OUTPUT){ fprintf(stdout,"\nCOINC COMPTON: %d %d %ld | %.1f %.1f %.1f | %d vs %d %d %ld | %.1f %.1f %.1f | %d, dt=%ld, sumE=%.1f\n",ptr->subsys,ptr->chan,ptr->ts,ptr->ecal,ptr->alt_ecal,ptr->esum,ptr->net_id,alt->subsys,alt->chan,alt->ts,alt->ecal,alt->alt_ecal,alt->esum,alt->net_id,(ptr->ts-alt->ts),(ptr->ecal+alt->ecal)); }

          qedE_ge_theta_sum_t->Fill(qedE_ge_theta_sum_t, alt_ecal, angle, 1);
          qed_geE_theta_sum_t->Fill(qed_geE_theta_sum_t, ptr_ecal, angle, 1);
        }
      }

      break;
      case SUBSYS_DCOMPTONA: // ge-DCOMPTONA where is a coincidence between a DSSD pixel and a HPGE addback with sum energy of 511keV
      // QED DCOMPTONA EVENTS
      // Ge with addback is a Double Compton scatter (DSSD-Ge-Ge)
      // Identified as subsys==SUBSYS_DCOMPTONA
      // In DCOMPTONA event, ecal will be QED_PIXEL and alt_ecal will be Ge addback sum energy
      // In DCOMPTONA event, crystal_table[chan]=pos will be QED DSSD number [1-6], alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
      // In DCOMPTONA event, net_id will be first HPGe crystal number [1-64], alt2_chan will be second HPGe crystal number [1-64]
      ge_dcs->Fill(ge_dcs, (int)ptr->ecal, (int)(alt->esum), 1);
      if(ptr->esum>ptr->ecal){
        geadd_dcs->Fill(geadd_dcs, (int)ptr->esum, (int)(alt->esum), 1);
      }

      pos1 = crystal_table[alt->chan];
      qed1 = (alt->alt_chan%1024);
      c1 = alt->net_id;
      c2 = alt->alt2_chan;

      c3 = crystal_table[ptr->chan];
      c4 = crystal_table[ptr->alt_chan];
      if( c1 >= 0 && c1 < 64 &&  c2 >= 0 && c2 < 64 &&  c3 >= 0 && c3 < 64 &&  c4 >= 0 && c4 < 64 ){
        if(c1 != c3 && c1 != c4 && c2 != c3 && c2 != c4){
          initial_theta = scattering_angle_QEDGe(pos1, qed1, c1);
          omega = angular_diff_QEDGe(pos1,qed1, c3, 110);
          azimuthal = azimuthal_TCS_GeGe_SiGeGe(c3, c4, c1, c2);
          //fprintf(stdout,"TCS %d %d %d %d | %d %d | %0.1f %0.1f\n",pos1, qed1, ge1, c2, c3, c4,initial_theta,azimuthal);


          // Weighting factor for delta-phi plot
          // Here use a coincidence between DSSD-GE (511keV) with the Ge-Ge (1274keV)
          // This will be an isotropic relationship to essentially count the active DSSD-Ge combinations
          if(ptr->esum>1274-6 && ptr->esum>1274+6){
            qed_wf_omega->Fill(qed_wf_omega, (int)omega, 1);
            if(omega>90){
              qed_wf_dcs_azi->Fill(qed_wf_dcs_azi, (int)azimuthal, 1);
            }
          }

          dcsa_cs_omega_ge->Fill(dcsa_cs_omega_ge, (int)(omega), 1);
          if(omega>159){
            dcsa_theta_azi_ge->Fill(dcsa_theta_azi_ge, (int)(initial_theta), (int)(azimuthal), 1);
          }
        }
      }
      break;
      case SUBSYS_DCOMPTONB: // ge-DCOMPTONB where is a coincidence between a DSSD pixel and a HPGE addback with sum energy of 511keV
      // QED DCOMPTONB EVENTS
      //  (DSSD-DSSD-Ge)
      // Identified as subsys==SUBSYS_DCOMPTONB
      // In DCOMPTONB event, ecal will be the first QED_PIXEL, alt_ecal will be the second QED_PIXEL, alt2_ecal will be Ge energy
      // In DCOMPTONB event, crystal_table[chan]=pos will be  first QED DSSD number [1-6],  alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
      // In DCOMPTONB event, crystal_table[ tof]=pos will be second QED DSSD number [1-6], alt2_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
      // In DCOMPTONB event, net_id will be the HPGe crystal number [1-64]

      pos1 = crystal_table[alt->chan];
      qed1 = (alt->alt_chan&1023);
      pos2 = crystal_table[alt->tof];
      qed2 = (alt->alt2_chan&1023);
      c1 = alt->net_id;

      c2 = crystal_table[ptr->chan];
      c3 = crystal_table[ptr->alt_chan];
      if( c1 >= 0 && c1 < 64 &&  c2 >= 0 && c2 < 64 &&  c3 >= 0 && c3 < 64 ){
        if(c1 != c2 && c1 != c3 && c2 != c3){
          initial_theta = (int)scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
          omega = (int)angular_diff_QEDGe(pos1,qed1, c2, 110);
          azimuthal = (int)azimuthal_TCS_GeGe_SiSiGe(pos2, qed2, c1, c2, c3);

          dcsb_cs_omega_ge->Fill(dcsb_cs_omega_ge, omega, 1);
          if(omega>159){
            dcsb_theta_azi_ge->Fill(dcsb_theta_azi_ge, initial_theta, azimuthal, 1);
          }
        }
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
        dt_tacs_hist[c1]->Fill(dt_tacs_hist[c1], (int)(abs_dt+(DT_SPEC_LENGTH>>1)), 1);
        tac_lbl_ts_diff[c1]->Fill(tac_lbl_ts_diff[c1], (int)((alt->ts-ptr->ts)+(DT_SPEC_LENGTH>>1)), 1);
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
        tac_lbl_ts_diff[c1]->Fill(tac_lbl_ts_diff[c1], (int)((alt->ts-ptr->ts)+(DT_SPEC_LENGTH>>1)), 1);
      } break;
      case SUBSYS_TAC_ART:
      c1=crystal_table[alt->chan]-1;  // assign c1 as TAC number
      if(c1 >= 9 && c1 < N_TACS ){ // 8 LBL + 1 ZDS + 4 ARIES
        tac_lbl_ts_diff[c1]->Fill(tac_lbl_ts_diff[c1], (int)((alt->ts-ptr->ts)+(DT_SPEC_LENGTH>>1)), 1);
      } break;
    }
    return(0);
  }

  int frag_hist[PTR_BUFSIZE];
  int fill_coinc_histos(int win_idx, int frag_idx)
  {
    int global_window_size = (int)(sort_window_width>>1); // size in grif-replay should be double this
    Grif_event *alt, *ptr, *original_ptr = &grif_event[win_idx], *tmp;
    int dt, abs_dt,  pos, c1, c2, index, ptr_swap, delta_cfd;
    int ptr_subsys, alt_subsys, ptr_ecal, alt_ecal, sum_ecal;
    int pos1, qed1, ge1, pos2, qed2, pos3, qed3, ge2, ge3, angle, angle1, angle2;  // QED variables
    double omega, theta1, theta2, delta_theta, azimuthal, azimuthal2, initial_theta; // QED variables
    double energy_derived_theta1, energy_derived_theta2;
    TH2I *hist_ee; TH1I *hist_dt; TH1I *hist_dcfd;

    // histogram of coincwin-size
    dt = (frag_idx - win_idx + 2*PTR_BUFSIZE) % PTR_BUFSIZE; ++frag_hist[dt];

    while( win_idx != frag_idx ){ // check all conicidences in window
      if( ++win_idx == PTR_BUFSIZE ){ win_idx = 0; } // wrap
      ptr = original_ptr;
      alt = &grif_event[win_idx];
      if( ptr->subsys > alt->subsys ){ tmp = ptr; ptr = alt; alt = tmp; ptr_swap = 1; }

      abs_dt = dt = ptr->ts - alt->ts; if( dt < 0 ){ abs_dt = -1*dt; }
      if( abs_dt > global_window_size ){ break; }

      ptr_subsys = ptr->subsys; alt_subsys = alt->subsys;
      ptr_ecal = ptr->ecal; alt_ecal = alt->ecal;
      // the usual subsys-vs-subsys 1d-time-diff and 2d-EvsE
      if( (hist_dt = subsys_dt[ptr_subsys][alt_subsys]) != NULL ){
        hist_dt->Fill(hist_dt, (int)(abs_dt+(DT_SPEC_LENGTH>>1)), 1);
      }
      if( (hist_dcfd = subsys_dcfd[ptr_subsys][alt_subsys]) != NULL ){
        hist_dcfd->Fill(hist_dcfd, (int)((ptr->cfd>>4)-(alt->cfd>>4)+(DT_SPEC_LENGTH>>1)), 1);
      }
      if( (hist_ee = subsys_e_vs_e[ptr_subsys][alt_subsys]) != NULL ){
        if((abs_dt >= time_diff_gate_min[ptr_subsys][alt_subsys]) && (abs_dt <= time_diff_gate_max[ptr_subsys][alt_subsys]) ){
          /*
          if(ptr->subsys == SUBSYS_QED_PIXEL){

          //  (int)(c2/N_QED_STRIPS); // QED p strip number [0-31]
          //  (c2%N_QED_STRIPS); // QED n strip number [0-31]
          fprintf(stdout,"Breakpoint QED-QED [%d,%d,%d,%d] [%d,%d,%d,%d] %d, %d, %d %d\n",crystal_table[ptr->chan],ptr->alt_chan,(int)(ptr->alt_chan/N_QED_STRIPS),
          crystal_table[alt->chan],alt->alt_chan,alt->alt_chan%N_QED_STRIPS,
          ptr->trig_acc,alt->trig_acc,(int)ptr->ecal, (int)alt->ecal);
        }
        */
        hist_ee->Fill(hist_ee, ptr_ecal, alt_ecal, 1);
      }
    }
    switch(ptr_subsys){ // No Nested switch - use separate functions if needed
      case SUBSYS_HPGE_A: fill_ge_coinc_histos(ptr,   alt, abs_dt); break;
      case SUBSYS_LABR_L: fill_labr_coinc_histos(ptr, alt, abs_dt); break;
      case SUBSYS_BGO:
      if(alt_subsys == SUBSYS_BGO){
        if((abs_dt >= time_diff_gate_min[SUBSYS_BGO][SUBSYS_BGO]) && (abs_dt <= time_diff_gate_max[SUBSYS_BGO][SUBSYS_BGO]) ){
          c1 = crystal_table[ptr->chan];
          c2 = crystal_table[alt->chan];
          bgobgo_hit->Fill(bgobgo_hit, c1, c2, 1);
        }}
        break;
        case SUBSYS_RCMP:
        if(alt_subsys == SUBSYS_RCMP){
          if( (abs_dt >= time_diff_gate_min[SUBSYS_RCMP][SUBSYS_RCMP]) && (abs_dt <= time_diff_gate_max[SUBSYS_RCMP][SUBSYS_RCMP]) ){
            if((pos = crystal_table[ptr->chan]) == crystal_table[alt->chan] &&
            polarity_table[ptr->chan] != polarity_table[alt->chan] ){ // front and back of same DSSD
              c1 = element_table[ptr->chan];
              c2 = element_table[alt->chan];
              rcmp_fb[(pos-1)]->Fill(rcmp_fb[(pos-1)], ptr_ecal, alt_ecal, 1);
              if(polarity_table[ptr->chan]==0){ rcmp_hit[(pos-1)]->Fill(rcmp_hit[(pos-1)], c1, c2, 1);
              }else{
                rcmp_hit[(pos-1)]->Fill(rcmp_hit[(pos-1)], c2, c1, 1);
              }
            }}} break;
            case SUBSYS_ARIES_A:
            if( alt_subsys == SUBSYS_ARIES_A && ptr_ecal>0 && alt_ecal>0 ){
              if((abs_dt >= time_diff_gate_min[SUBSYS_ARIES_A][SUBSYS_ARIES_A]) && (abs_dt <= time_diff_gate_max[SUBSYS_ARIES_A][SUBSYS_ARIES_A])){
                c1 = crystal_table[ptr->chan];
                c2 = crystal_table[alt->chan];
                if( c1 >= 1 && c1 <=76 && c2 >= 1 && c2 <=76 ){
                  aa_hit->Fill(aa_hit, c1, c2, 1);
                }
              }}
              if( alt_subsys == SUBSYS_TAC_LABR && crystal_table[alt->chan] == 8 ){ // ARIES TAC
                // sum tac spectrum including all art
                tac_aries_art_sum->Fill(tac_aries_art_sum, alt_ecal, 1);
                c2 = crystal_table[ptr->chan]-1;
                if( c2>=0 && c2<N_ARIES ){ // tac spectrum per ART tiles
                  tac_aries_art[c2]->Fill(tac_aries_art[c2], alt_ecal, 1);
                }
              } break;
              case SUBSYS_ARIES_B:// ARIES Fast Output in CAEN
              if(alt_subsys == SUBSYS_DESWALL){ // aries-DSW
                art_dsw->Fill(art_dsw, ptr_ecal, (int)alt->alt_ecal, 1);
                desw_sum_e_b->Fill(desw_sum_e_b, alt_ecal, 1);
                desw_sum_tof_b->Fill(desw_sum_tof_b, (int)alt->alt_ecal, 1); // alt_ecal = corrected time-of-flight
              } break;
              case SUBSYS_ZDS_A: // grif16 zds
              if(alt_subsys == SUBSYS_ZDS_B ){ // ZDS GRIF-CAEN coincidence
                gc_hist->Fill(gc_hist, 3, 1);
                gc_hist->Fill(gc_hist, 5, 1);
                dt_hist[23]->Fill(dt_hist[23], (int)(abs(ptr->cfd - alt->cfd)+(DT_SPEC_LENGTH>>1)), 1);
              } break;
              case SUBSYS_DESWALL:
              if( alt_subsys == SUBSYS_DESWALL ){
                dt_hist[24]->Fill(dt_hist[24], (int)(abs(ptr->cfd - alt->cfd)+(DT_SPEC_LENGTH>>1)), 1);
                c1 = crystal_table[ptr->chan]-1; c2 = crystal_table[alt->chan]-1;
                if( c1 >= 0 && c1 < 60 && c2 >= 0 && c2 < 60 ){
                  dsw_hit->Fill(dsw_hit, c1, c2, 1);
                  dsw_dsw->Fill(dsw_dsw, (int)ptr->alt_ecal, (int)alt->alt_ecal, 1);
                  // Fold 2 sum spectra
                  desw_sum_e_nn->Fill(desw_sum_e_nn, ptr_ecal, 1);
                  desw_sum_tof_nn->Fill(desw_sum_tof_nn, (int)ptr->alt_ecal, 1);// alt_ecal = corrected time-of-flight
                  // DSW-DSW angular correlations
                  // Fill the appropriate angular bin spectrum with the corrected time-of-flight value
                  index = DSW_DSW_angles[c1][c2];
                  dsw_angcor[index]->Fill(dsw_angcor
                    [index],(int)ptr->alt_ecal,(int)alt->alt_ecal, 1);
                    // Fold 2, angle greater than 60 degrees, sum spectra
                    // index 13 = 58.555, index 14 = 61.535
                    if( index > 13 ){                                         // alt_ecal = corrected time-of-flight
                      desw_sum_e_nn_a->Fill(desw_sum_e_nn_a, ptr_ecal, 1);
                      desw_sum_tof_nn_a->Fill(desw_sum_tof_nn_a, (int)ptr->alt_ecal, 1);
                    }
                  }
                }
                if( alt_subsys == SUBSYS_ZDS_B ){ // ZDS-DSW
                  dt_hist[25]->Fill(dt_hist[25], (int)(abs(ptr->cfd - alt->cfd)+(DT_SPEC_LENGTH>>1)), 1);
                  desw_sum_e_b->Fill(desw_sum_e_b, ptr_ecal, 1);
                  desw_sum_tof_b->Fill(desw_sum_tof_b, (int)ptr->alt_ecal, 1); // alt_ecal = corrected time-of-flight

                  desw_q_tof->Fill(desw_q_tof, (int)ptr->q1, (int)ptr->alt_ecal, 1);
                  desw_cc_tof->Fill(desw_cc_tof, (int)ptr->cc_short, (int)ptr->alt_ecal, 1);
                } break;

                case SUBSYS_QED_PIXEL:
                if( alt_subsys == SUBSYS_QED_PIXEL ){  // QED-QED
                  pos1  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
                  qed1 = ptr->alt_chan; // QED pixel number [0-1023]
                  pos2  = crystal_table[alt->chan]; // QED DSSD number [1-6]
                  qed2 = alt->alt_chan; // QED pixel number [0-1023]
                  //  (int)(c2/N_QED_STRIPS); // QED p strip number [0-31]
                  //  (c2%N_QED_STRIPS); // QED n strip number [0-31]
                  if( qed1 >= 0 && qed1 < 1024 && qed2 >= 0 && qed2 < 1024 && ptr->ecal>5 && alt->ecal>5){ //&& (abs_dt >= time_diff_gate_min[SUBSYS_QED_PIXEL][SUBSYS_QED_PIXEL]) && (abs_dt <= time_diff_gate_max[SUBSYS_QED_PIXEL][SUBSYS_QED_PIXEL])){
                    if(pos1 == 2 && pos2 == 3){ // Back-to-back DSSD, QED2 and QED3
                      qed_qed_23->Fill(qed_qed_23, (int)ptr->ecal, (int)alt->ecal, 1);
                      delta_cfd = (ptr->cfd>>4) - (alt->cfd>>4);
                      qed_qed_23dt->Fill(qed_qed_23dt, delta_cfd+512, 1);
                      angle = (int)scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
                      sum_ecal = (int)(ptr->ecal+alt->cfd);
                      qed_qed_23_theta2->Fill(qed_qed_23_theta2, (int)ptr->ecal, angle, 1);
                      qed_qed_23_theta3->Fill(qed_qed_23_theta3, (int)alt->ecal, angle, 1);
                      qed_qed_23_totv2->Fill(qed_qed_23_totv2, (int)ptr->ecal, sum_ecal, 1);
                      qed_qed_23_totv3->Fill(qed_qed_23_totv3, (int)alt->ecal, sum_ecal, 1);
                    }else if(pos2 == 2 && pos1 == 3){
                      qed_qed_23->Fill(qed_qed_23, (int)alt->ecal, (int)ptr->ecal, 1);
                      delta_cfd = (ptr->cfd>>4) - (alt->cfd>>4);
                      qed_qed_23dt->Fill(qed_qed_23dt, delta_cfd+512, 1);
                      angle = (int)scattering_angle_QEDQED(pos2, qed2, pos1, qed1);
                      sum_ecal = (int)(ptr->ecal+alt->ecal);
                      qed_qed_23_theta2->Fill(qed_qed_23_theta2, (int)alt->ecal, angle, 1);
                      qed_qed_23_theta3->Fill(qed_qed_23_theta3, (int)ptr->ecal, angle, 1);
                      qed_qed_23_totv2->Fill(qed_qed_23_totv2, (int)alt->ecal, sum_ecal, 1);
                      qed_qed_23_totv3->Fill(qed_qed_23_totv3, (int)ptr->ecal, sum_ecal, 1);
                    }
                    if(pos1 == 1 && pos2 == 2){ // 90 degree DSSD, QED1 and QED2
                      qed_qed_12->Fill(qed_qed_12, (int)ptr->ecal, (int)alt->ecal, 1);
                      delta_cfd = (ptr->cfd>>4) - (alt->cfd>>4);
                      qed_qed_12dt->Fill(qed_qed_12dt, delta_cfd+512, 1);
                      angle = (int)scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
                      sum_ecal = (int)(ptr->ecal+alt->ecal);
                      qed_qed_12_theta1->Fill(qed_qed_12_theta1, (int)ptr->ecal, angle, 1);
                      qed_qed_12_theta2->Fill(qed_qed_12_theta2, (int)alt->ecal, angle, 1);
                      qed_qed_12_totv1->Fill(qed_qed_12_totv1, (int)ptr->ecal, sum_ecal, 1);
                      qed_qed_12_totv2->Fill(qed_qed_12_totv2, (int)alt->ecal, sum_ecal, 1);
                    }else if(pos2 == 1 && pos1 == 2){
                      qed_qed_12->Fill(qed_qed_12, (int)alt->ecal, (int)ptr->ecal, 1);
                      delta_cfd = (ptr->cfd>>4) - (alt->cfd>>4);
                      qed_qed_12dt->Fill(qed_qed_12dt, delta_cfd+512, 1);
                      angle = (int)scattering_angle_QEDQED(pos2, qed2, pos1, qed1);
                      sum_ecal = (int)(ptr->ecal+alt->ecal);
                      qed_qed_12_theta1->Fill(qed_qed_12_theta1, (int)alt->ecal, angle, 1);
                      qed_qed_12_theta2->Fill(qed_qed_12_theta2, (int)ptr->ecal, angle, 1);
                      qed_qed_12_totv1->Fill(qed_qed_12_totv1, (int)alt->ecal, sum_ecal, 1);
                      qed_qed_12_totv2->Fill(qed_qed_12_totv2, (int)ptr->ecal, sum_ecal, 1);
                    }
                    if(pos1 == 1 && pos2 == 4){ // Opposite DSSD, QED1 and QED4
                      qed_qed_14->Fill(qed_qed_14, (int)ptr->ecal, (int)alt->ecal, 1);
                      delta_cfd = (ptr->cfd>>4) - (alt->cfd>>4);
                      qed_qed_14dt->Fill(qed_qed_14dt, delta_cfd+512, 1);
                      angle = (int)scattering_angle_QEDQED(pos1, qed1, pos2, qed2);
                      sum_ecal = (int)(ptr->ecal+alt->ecal);
                      qed_qed_14_theta1->Fill(qed_qed_14_theta1, (int)ptr->ecal, angle, 1);
                      qed_qed_14_theta4->Fill(qed_qed_14_theta4, (int)alt->ecal, angle, 1);
                      qed_qed_14_totv1->Fill(qed_qed_14_totv1, (int)ptr->ecal, sum_ecal, 1);
                      qed_qed_14_totv4->Fill(qed_qed_14_totv4, (int)alt->ecal, sum_ecal, 1);
                    }else if(pos2 == 1 && pos1 == 4){
                      qed_qed_14->Fill(qed_qed_14, (int)alt->ecal, (int)ptr->ecal, 1);
                      delta_cfd = (ptr->cfd>>4) - (alt->cfd>>4);
                      qed_qed_14dt->Fill(qed_qed_14dt, delta_cfd+512, 1);
                      angle = (int)scattering_angle_QEDQED(pos2, qed2, pos1, qed1);
                      sum_ecal = (int)(ptr->ecal+alt->ecal);
                      qed_qed_14_theta1->Fill(qed_qed_14_theta1, (int)alt->ecal, angle, 1);
                      qed_qed_14_theta4->Fill(qed_qed_14_theta4, (int)ptr->ecal, angle, 1);
                      qed_qed_14_totv1->Fill(qed_qed_14_totv1, (int)alt->ecal, sum_ecal, 1);
                      qed_qed_14_totv4->Fill(qed_qed_14_totv4, (int)ptr->ecal, sum_ecal, 1);
                    }
                  }
                }
                break;
                case SUBSYS_COMPTON:
                if( alt_subsys == SUBSYS_COMPTON ){ // COMPTON-COMPTON where is a coincidence between a DSSD pixel and a HPGE with sum energy of 511keV
                  // QED COMPTON EVENTS
                  // Identified as subsys==SUBSYS_COMPTON
                  // pos is HPGE crystal number
                  // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL
                  // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
                  //  pos  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
                  //  c2 = (ptr->alt_chan%1024);        // Pixel number [0-1023]
                  //  c1 = ptr->net_id; // HPGe crystal number
                  pos1 = crystal_table[ptr->chan];
                  qed1 = (ptr->alt_chan&1023);
                  ge1 = ptr->net_id;
                  pos2 = crystal_table[alt->chan];
                  qed2 = (alt->alt_chan&1023);
                  ge2 = alt->net_id;
                  if(DEBUG_OUTPUT){ fprintf(stdout,"COMPTON-COMPTON: pos1,qed1,ge1 %d,%d,%d | pos2,qed2,ge2 %d,%d,%d | %.1f %.1f\n",pos1,qed1,ge1,pos2,qed2,ge2,ptr->esum,alt->esum); }
                  omega = angular_diff_QEDQED(pos1, qed1, pos2, qed2);
                  qed_dcs_omega->Fill(qed_dcs_omega, (int)omega, 1);
                  qed_dcs_omega_dt->Fill(qed_dcs_omega_dt, dt+512, (int)omega, 1);
                  qedx_dcs_omega_dt[pos1-1]->Fill(qedx_dcs_omega_dt[pos1-1], dt+512, (int)omega, 1);
                  if(pos1 != pos2 && ge1 != ge2){
                    qed_dcs_omega_dtx->Fill(qed_dcs_omega_dtx, dt+512, (int)omega, 1);
                    qed_dcs_omega_t->Fill(qed_dcs_omega_t, (int)omega, 1);
                    theta1 = scattering_angle_QEDGe(pos1, qed1, ge1);
                    theta2 = scattering_angle_QEDGe(pos2, qed2, ge2);
                    azimuthal = azimuthal_DCS(pos1, qed1, ge1, pos2, qed2, ge2);
                    delta_theta = (theta1>theta2) ? theta1-theta2 : theta2-theta1;
                    //    if(){  }
                    //  fprintf(stdout,"omega = %f, theta1 = %f, theta2 = %f, azimuth = %f, for [%d %d %d] - [%d %d %d]\n",omega,theta1,theta2,azimuthal,pos1, qed1, ge1, pos2, qed2, ge2);


                    if(omega>170){ // Require back-to-back coincidence
                      // Calculate the other angles
                      azimuthal2 = energy_corrected_azimuthal_DCS(pos1, qed1, ge1, ptr->alt_ecal,pos2, qed2, ge2, alt->alt_ecal);
                      energy_derived_theta1 = compton_angle(ptr->alt_ecal, 511.0);
                      energy_derived_theta2 = compton_angle(alt->alt_ecal, 511.0);

                      if(dt>-11 && dt<1){ // Prompt time coincidence: -10ns to -110ns

                        if(theta1<theta2){
                          qed_theta1_vs_theta2->Fill(qed_theta1_vs_theta2, (int)theta1, (int)theta2, 1);
                        }else{
                          qed_theta1_vs_theta2->Fill(qed_theta1_vs_theta2, (int)theta2, (int)theta1, 1);
                        }

                        // Histogram to calculate weighting factors post sorting
                        // Here we remember the associated Ge for each qed pixel that was in a true coincidence
                        qed_ge_weight->Fill(qed_ge_weight,((pos1-1)*1024)+qed1, ge1, 1);
                        qed_ge_weight->Fill(qed_ge_weight,((pos2-1)*1024)+qed2, ge2, 1);

                        qed_delta_theta1_theta2->Fill(qed_delta_theta1_theta2, (int)delta_theta, 1);
                        qed_sum_theta1_theta2->Fill(qed_sum_theta1_theta2, (int)(theta1+theta2), 1);
                        qed_theta1_azi->Fill(qed_theta1_azi, (int)theta1, (int)azimuthal, 1);
                        qed_theta2_azi->Fill(qed_theta2_azi, (int)theta2, (int)azimuthal, 1);

                        if(energy_derived_theta1<energy_derived_theta2){
                          qed2_theta1_vs_theta2->Fill(qed2_theta1_vs_theta2, (int)energy_derived_theta1, (int)energy_derived_theta2, 1);
                        }else{
                          qed2_theta1_vs_theta2->Fill(qed2_theta1_vs_theta2, (int)energy_derived_theta2, (int)energy_derived_theta1, 1);
                        }
                        qed2_theta1_azi->Fill(qed2_theta1_azi, (int)energy_derived_theta1, (int)azimuthal2, 1);
                        qed2_theta2_azi->Fill(qed2_theta2_azi, (int)energy_derived_theta2, (int)azimuthal2, 1);

                        qed_dcs_azi_t->Fill(qed_dcs_azi_t, (int)azimuthal, 1);
                        // Scattering angle 70 to 110
                        if(theta1>69 && theta1<111 && theta2>69 && theta2<111){
                          comp_comp->Fill(comp_comp, (int)ptr->esum, (int)(alt->esum), 1);
                          qed_dcs_azi->Fill(qed_dcs_azi, (int)azimuthal, 1);
                          // Scattering angle 93 to 103
                          if(theta1>92 && theta1<104 && theta2>92 && theta2<104){
                            qed_dcs_azi_tg->Fill(qed_dcs_azi_tg, (int)azimuthal, 1);
                          }
                        }

                        // I know this code is ugly. Sorry.
                        if(energy_derived_theta1>=0 && energy_derived_theta1<=180 && energy_derived_theta2>=0 && energy_derived_theta2<=180){
                          qed_dcs_azi_bins1->Fill(qed_dcs_azi_bins1, (int)(azimuthal2), 1);

                          if(energy_derived_theta1>=10 && energy_derived_theta1<=170 && energy_derived_theta2>=10 && energy_derived_theta2<=170){
                            qed_dcs_azi_bins2->Fill(qed_dcs_azi_bins2, (int)(azimuthal2), 1);

                            if(energy_derived_theta1>=20 && energy_derived_theta1<=160 && energy_derived_theta2>=20 && energy_derived_theta2<=160){
                              qed_dcs_azi_bins3->Fill(qed_dcs_azi_bins3, (int)(azimuthal2), 1);

                              if(energy_derived_theta1>=30 && energy_derived_theta1<=150 && energy_derived_theta2>=30 && energy_derived_theta2<=150){
                                qed_dcs_azi_bins4->Fill(qed_dcs_azi_bins4, (int)(azimuthal2), 1);

                                if(energy_derived_theta1>=40 && energy_derived_theta1<=140 && energy_derived_theta2>=40 && energy_derived_theta2<=140){
                                  qed_dcs_azi_bins5->Fill(qed_dcs_azi_bins5, (int)(azimuthal2), 1);

                                  if(energy_derived_theta1>=50 && energy_derived_theta1<=130 && energy_derived_theta2>=50 && energy_derived_theta2<=130){
                                    qed_dcs_azi_bins6->Fill(qed_dcs_azi_bins6, (int)(azimuthal2), 1);

                                    if(energy_derived_theta1>=60 && energy_derived_theta1<=120 && energy_derived_theta2>=60 && energy_derived_theta2<=120){
                                      qed_dcs_azi_bins7->Fill(qed_dcs_azi_bins7, (int)(azimuthal2), 1);

                                      if(energy_derived_theta1>=70 && energy_derived_theta1<=110 && energy_derived_theta2>=70 && energy_derived_theta2<=110){
                                        qed_dcs_azi_bins8->Fill(qed_dcs_azi_bins8, (int)(azimuthal2), 1);

                                        if(energy_derived_theta1>=93 && energy_derived_theta1<=103 && energy_derived_theta2>=93 && energy_derived_theta2<=103){
                                          qed_dcs_azi_bins8a->Fill(qed_dcs_azi_bins8a, (int)(azimuthal2), 1);
                                        }

                                        if(energy_derived_theta1>=80 && energy_derived_theta1<=100 && energy_derived_theta2>=80 && energy_derived_theta2<=100){
                                          qed_dcs_azi_bins9->Fill(qed_dcs_azi_bins9, (int)(azimuthal2), 1);

                                          if(energy_derived_theta1>=85 && energy_derived_theta1<=95 && energy_derived_theta2>=85 && energy_derived_theta2<=95){
                                            qed_dcs_azi_bins10->Fill(qed_dcs_azi_bins10, (int)(azimuthal2), 1);
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

                      }else{ // Time-random for weighting factors
                        // Build weighting factors here

                        // This filling has been moved to the pre_sort_qed_weights function
                        /*
                        qed_dcs_azi_TRWF_t->Fill(qed_dcs_azi_TRWF_t, (int)azimuthal, 1);
                        // Scattering angle 70 to 110
                        if(theta1>69 && theta1<111 && theta2>69 && theta2<111){
                        qed_dcs_azi_TRWF->Fill(qed_dcs_azi_TRWF, (int)azimuthal, 1);
                        // Scattering angle 93 to 103
                        if(theta1>92 && theta1<104 && theta2>92 && theta2<104){
                        qed_dcs_azi_TRWF_tg->Fill(qed_dcs_azi_TRWF_tg, (int)azimuthal, 1);
                      }
                    }
                    // I know this code is ugly. Sorry.
                    if(energy_derived_theta1>=0 && energy_derived_theta1<=180 && energy_derived_theta2>=0 && energy_derived_theta2<=180){
                    qed_dcs_azi_TRWF_bins1->Fill(qed_dcs_azi_TRWF_bins1, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=10 && energy_derived_theta1<=170 && energy_derived_theta2>=10 && energy_derived_theta2<=170){
                    qed_dcs_azi_TRWF_bins2->Fill(qed_dcs_azi_TRWF_bins2, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=20 && energy_derived_theta1<=160 && energy_derived_theta2>=20 && energy_derived_theta2<=160){
                    qed_dcs_azi_TRWF_bins3->Fill(qed_dcs_azi_TRWF_bins3, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=30 && energy_derived_theta1<=150 && energy_derived_theta2>=30 && energy_derived_theta2<=150){
                    qed_dcs_azi_TRWF_bins4->Fill(qed_dcs_azi_TRWF_bins4, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=40 && energy_derived_theta1<=140 && energy_derived_theta2>=40 && energy_derived_theta2<=140){
                    qed_dcs_azi_TRWF_bins5->Fill(qed_dcs_azi_TRWF_bins5, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=50 && energy_derived_theta1<=130 && energy_derived_theta2>=50 && energy_derived_theta2<=130){
                    qed_dcs_azi_TRWF_bins6->Fill(qed_dcs_azi_TRWF_bins6, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=60 && energy_derived_theta1<=120 && energy_derived_theta2>=60 && energy_derived_theta2<=120){
                    qed_dcs_azi_TRWF_bins7->Fill(qed_dcs_azi_TRWF_bins7, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=70 && energy_derived_theta1<=110 && energy_derived_theta2>=70 && energy_derived_theta2<=110){
                    qed_dcs_azi_TRWF_bins8->Fill(qed_dcs_azi_TRWF_bins8, (int)(azimuthal2), 1);

                    if(energy_derived_theta1>=93 && energy_derived_theta1<=103 && energy_derived_theta2>=93 && energy_derived_theta2<=103){
                    qed_dcs_azi_TRWF_bins8a->Fill(qed_dcs_azi_TRWF_bins8a, (int)(azimuthal2), 1);
                  }

                  if(energy_derived_theta1>=80 && energy_derived_theta1<=100 && energy_derived_theta2>=80 && energy_derived_theta2<=100){
                  qed_dcs_azi_TRWF_bins9->Fill(qed_dcs_azi_TRWF_bins9, (int)(azimuthal2), 1);

                  if(energy_derived_theta1>=85 && energy_derived_theta1<=95 && energy_derived_theta2>=85 && energy_derived_theta2<=95){
                  qed_dcs_azi_TRWF_bins10->Fill(qed_dcs_azi_TRWF_bins10, (int)(azimuthal2), 1);
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
*/

}

}// end of omega>170

}
}
break;
case SUBSYS_DCOMPTONA:
if( alt_subsys == SUBSYS_COMPTON ){ // COMPTON-DCOMPTONA where is a coincidence between a DSSD pixel and a HPGE with sum energy of 511keV
  // QED COMPTON EVENTS
  // Identified as subsys==SUBSYS_COMPTON
  // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL
  // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
  //  pos  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
  //  c2 = (ptr->alt_chan%1024);       // Pixel number [0-1023]
  //  c1 = ptr->net_id;                // HPGe crystal number
  // QED DCOMPTONA EVENTS
  // Ge with addback is a Double Compton scatter (DSSD-Ge-Ge)
  // Identified as subsys==SUBSYS_DCOMPTONA
  // In DCOMPTONA event, ecal will be QED_PIXEL and alt_ecal will be Ge addback sum energy
  // In DCOMPTONA event, crystal_table[chan]=pos will be QED DSSD number [1-6], alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
  // In DCOMPTONA event, net_id will be first HPGe crystal number [1-64], alt2_chan will be second HPGe crystal number [1-64]

  pos1 = crystal_table[alt->chan];
  qed1 = (alt->alt_chan&1023);
  ge1 = alt->net_id;

  pos2 = crystal_table[ptr->chan];
  qed2 = (ptr->alt_chan&1023);
  ge2 = ptr->net_id;
  ge3 = ptr->alt2_chan;

  if(pos1 != pos2 && ge1 != ge2 && ge1 != ge3){
    omega = angular_diff_QEDQED(pos1, qed1, pos2, qed2);
    initial_theta = scattering_angle_QEDGe(pos1, qed1, ge1);
    // Here theta1 could be calculated from the energies in the two Ge crystals.
    theta2 = scattering_angle_QEDGe(pos2, qed2, ge2);
    azimuthal = azimuthal_TCS_SiGe_SiGeGe(pos1, qed1, ge1, ge2, ge3);
    //  fprintf(stdout,"TCS %d %d %d | %d %d %d %d | %0.1f %0.1f\n",pos1, qed1, ge1, pos2, qed2, ge2, ge3,initial_theta,azimuthal);
    comp_dcs->Fill(comp_dcs, (int)ptr->esum, (int)(alt->esum), 1);

    dcsa_cs_omega->Fill(dcsa_cs_omega, (int)(omega), 1);
    if(omega>159){
      dcsa_theta_azi->Fill(dcsa_theta_azi, (int)(initial_theta), (int)(azimuthal), 1); // This one looks good. Endorsement of azimuthal_TCS_SiGe_SiGeGe
    }
  }
}
break;
case SUBSYS_DCOMPTONB:
if( alt_subsys == SUBSYS_COMPTON ){ // COMPTON-DCOMPTONB where is a coincidence between two DSSD pixels and a HPGE with total energy of 511keV
  // QED COMPTON EVENTS
  // Identified as subsys==SUBSYS_COMPTON
  // In COMPTON event, ecal will be Ge and alt_ecal will be QED_PIXEL
  // In COMPTON event, crystal_table[chan]=pos will be Ge, alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
  //  pos  = crystal_table[ptr->chan]; // QED DSSD number [1-6]
  //  c2 = (ptr->alt_chan%1024);       // Pixel number [0-1023]
  //  c1 = ptr->net_id;                // HPGe crystal number
  // QED DCOMPTONB EVENTS
  //  (DSSD-DSSD-Ge)
  // Identified as subsys==SUBSYS_DCOMPTONB
  // In DCOMPTONB event, ecal will be the first QED_PIXEL, alt_ecal will be the second QED_PIXEL, alt2_ecal will be Ge energy
  // In DCOMPTONB event, crystal_table[chan]=pos will be  first QED DSSD number [1-6],  alt_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
  // In DCOMPTONB event, crystal_table[ tof]=pos will be second QED DSSD number [1-6], alt2_chan will be [DSSD*PIXELnumber],[0-5*0-1023]
  // In DCOMPTONB event, net_id will be the HPGe crystal number [1-64]

  pos1 = crystal_table[alt->chan];
  qed1 = (alt->alt_chan&1023);
  ge1 = alt->net_id;

  pos2 = crystal_table[ptr->chan];
  qed2 = (ptr->alt_chan&1023);
  pos3 = crystal_table[ptr->tof];
  qed3 = (ptr->alt2_chan&1023);
  ge2 = ptr->net_id;

  if(pos1 != pos2 && pos1 != pos3 && pos2 != pos3 && ge1 != ge2){
    initial_theta = (int)scattering_angle_QEDQED(pos2, qed2, pos3, qed3);
    omega = (int)angular_diff_QEDQED(pos1, qed1, pos2, qed2);
    azimuthal = (int)azimuthal_TCS_SiGe_SiSiGe(pos1, qed1, ge1, pos3, qed3, ge2);

    dcsb_cs_omega->Fill(dcsb_cs_omega, omega, 1);
    if(omega>159){
      dcsb_theta_azi->Fill(dcsb_theta_azi, initial_theta, azimuthal, 1);
    }
  }
}
break;
} // end switch
}
return(0);
}
