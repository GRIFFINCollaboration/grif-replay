// 0) make sure user sort runs
//
// 1) reset bytes sorted
//     ** Runs are entries in filelist[+may contain multiple subruns] **
//     so handle subruns differently
//     status - subrun is x/y or 1/1 ****
//
// DONE 2) strings replaced with %20 - globals etc
//
// 3) getHistofileList - Histogram files need run title info [+date/time]
//                     [Already have odbstuff in config file]
//    ALSO remember first and last midas-timestamps
//      add config section: SortedFileDetails
//
// DONE 4) Add radware matrix script [optionally all histos]
//
// 5) Add 180-deg-coinc matrix - for summing calcns
//       56Co - 4400keV   1keV/chan
//
//    Add 2-photon decay 1dhisto of sum of pairs of crystals [2photonGe]
//       Add 2-photon Variable as well
//
// *** 6) sort gzipped files
//
// DONE 7) check root script with config file entry
//
// DONE 8) Reorder histo list - Energy Time Waveform Pulsehight
//
// DONE 9) Fix Gains from config not midas
//
// odb-error
//
// Sr90[]
//
/////////////////////////////////////////////////////////////////////////////

/* new js fn similar to odbget - will have to decode args in server.c */
/* Test commands ...

wget 'http://grifstore1:9093/?cmd=addDatafile&filename=/tig/grifstore1b/grifalt/schedule145/Dec2023/run21834_000.mid'

wget 'http://grifstore1:9093/?cmd=getDatafileList&dir=/tig/grifstore1b/grifalt/schedule145/Dec2023'

wget 'http://localhost:9093/?cmd=getDatafileList&dir=.'
wget 'http://localhost:9093/?cmd=getHistofileList&dir=.'
wget 'http://localhost:9093/?cmd=getSortStatus'
wget 'http://localhost:9093/?cmd=addDatafile&filename=/tig/grifstore1b/grifalt/schedule145/Dec2023/run21834_000.mid'
wget 'http://localhost:9093/?cmd=openHistofile&filename=histos.tar'


wget 'http://localhost:9093/?cmd=addGate&filename=histos.cfg'
wget 'http://localhost:9093/?cmd=saveConfig&filename=histos.cfg'



wget 'http://panther:9093/?cmd=addDatafile&filename=/tig/grifstore0b/griffin/schedule140/Calibrations-July2021/run18909_020.mid'

wget 'http://panther:9093/?cmd=endCurrentFile'

wget 'http://panther:9093/?cmd=getSpectrumList&filename=run18909_020.tar'

wget 'http://panther:9093/?cmd=callspechandler&spectum0=Hitpattern_Energy&spectrum1=GRG01BN00A_Energy'

*/
//////////////////////////////////////////////////////////////////////////////
//"spectrum1d/index.html"
/*    check for javascript cmds or send file               */
/*    file should be under custom pages - check url for .. */
// Format of CALLSPECHANDLER call is:
// host:PORT/?cmd=callSpectrumHandler&spectum0=firstSpec&spectrum1=secondSpec
// CURRENTLY DEFINED COMMANDS ...
//     /?cmd=getDatafileList&dir=XXXX
//     /?cmd=getHistofileList
//     /?cmd=getSortStatus
//     /?cmd=getSpectrumList
//     /?cmd=callSpectrumHandler
//     /?cmd=addDatafile&filename=XXXX
//
// ALSO FIXED URLs ...
//     /filter-status.html
//     /report
//     /*.css
//     /*.js
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>
#include "config.h"
#include "grif-replay.h"
#include "histogram.h"
#include "web_server.h"

///////////////////////////////////////////////////////////////////////////
//////////////         Url  Command  interpreter          /////////////////
///////////////////////////////////////////////////////////////////////////

extern int coinc_events_cutoff;

int handle_command(int fd, int narg, char url_args[][STRING_LEN])
{
   char *ptr = url_args[0]; // url_args already processed - dont need strncmp
   int i, j, xbins, xmin, xmax, ybins, ymin, ymax, value, val2;
   char *histodir, *configdir, *calsrc, *host, *expt, tmp[128];
   extern volatile int shutdown_server;
   Sort_status *sort_stat;
   static FILE *web_fp;
   FILE *tmp_fp;
   Config *cfg;

   if( strcmp(ptr, "getSpectrumList") == 0 ){ /* --- list spectra -- */
      if( strncmp(url_args[2],"filename",8) == 0 ){
         send_spectrum_list(url_args[3], fd);
      } else {
         send_spectrum_list("", fd);
      }
   } else
   if( strcmp(ptr, "getDatafileList") == 0 ){/* -- list datafiles -- */
      if( strcmp(url_args[2], "dir") == 0 ){
         send_datafile_list(url_args[3], fd, 0);
      } else {
         sprintf(tmp,"bad argument:%s in getdatafileList\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in getdatafileList\n", url_args[2]);
      }
   } else
   if( strcmp(ptr, "getDatafileDetails") == 0 ){/* -- datafile titles -- */
      if( strcmp(url_args[2], "dir") == 0 ){
         send_datafile_list(url_args[3], fd, 1);
      } else {
         sprintf(tmp,"bad argument:%s in getdatafileList\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in getdatafileList\n", url_args[2]);
      }
   } else
   if( strcmp(ptr, "getHistofileList") == 0 ){ /* ----list histos---- */
      if( strcmp(url_args[2], "dir") == 0 ){
         send_histofile_list(url_args[3], fd);
      } else {
         sprintf(tmp,"bad argument:%s in gethistofileList\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in gethistofileList\n", url_args[2]);
      }
   } else
   if( strcmp(ptr, "getConfigfileList") == 0 ){ /* ----list configs---- */
      if( strcmp(url_args[2], "dir") == 0 ){
         send_configfile_list(url_args[3], fd);
      } else {
         sprintf(tmp,"bad argument:%s in getConfigfileList\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in getConfigfileList\n", url_args[2]);
      }
   } else
   if( strcmp(ptr, "getSortStatus") == 0 ){ /* -----------------------*/
      send_sort_status(fd);
   } else
   if( strcmp(ptr, "terminateServer") == 0 ){ /* -----------------------*/
      send_header(fd, APP_JSON);
      shutdown_server = 1; return(0);
   } else
   if( strcmp(ptr, "setSortStatus") == 0 ){ /* -----------------------*/
      if( strcmp(url_args[2], "status") != 0 ){
         sprintf(tmp,"bad argument:%s in setSortStatus\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in setSortStatus\n", url_args[2]);
      }
      if( sscanf(url_args[3],"%d", &value) < 1 ){
         sprintf(tmp,"setSortStatus:cant Read value: %s\n",url_args[3]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"setSortStatus:cantRead val:%s\n",url_args[3]);
         return(-1);
      }
      send_header(fd, APP_JSON);
      sort_stat = get_sort_status();
      sort_stat->reorder       =  value     & 3;
      sort_stat->single_thread = (value>>2) & 1;
      sort_stat->sort_thread   = (value>>3) & 1;
      return(0);
   } else
   if( strcmp(ptr, "setCoincLimit") == 0 ){ /* -----------------------*/
      if( strcmp(url_args[2], "limit") != 0 ){
         sprintf(tmp,"bad argument:%s in setCoincLimit\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in setCoincLimit\n", url_args[2]);
      }
      if( sscanf(url_args[3],"%d", &value) < 1 ){
         sprintf(tmp,"setCoincLimit:cant Read value: %s\n",url_args[3]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"setCoincLimit:cantRead val:%s\n",url_args[3]);
         return(-1);
      }
      send_header(fd, APP_JSON);
      coinc_events_cutoff = value;
      return(0);
   } else
   if( strcmp(ptr, "endCurrentFile") == 0 ){ /* -----------------------*/
      send_header(fd, APP_JSON);
      end_current_sortfile(fd);
   } else
   if( strcmp(ptr, "openHistofile") == 0 ){ /* ---------------------- */
      if( strncmp(url_args[2],"filename",8) == 0 ){
         send_header(fd, APP_JSON);
         read_histofile(url_args[3], 0);
      } else {
         sprintf(tmp,"bad argument:%s in openHistoFile\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in openHistoFile\n", url_args[2]);
      }
   } else
   if( strcmp(ptr, "addDatafile") == 0 ){ /* ---- sort file ---------*/
      histodir = configdir = calsrc = host = expt = NULL;
      if( strcmp(url_args[2], "filename") != 0 ){
         sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[2]);
      }
      if( strcmp(url_args[3],"ONLINE") == 0 ){
         if( narg > 4 ){
            if( strcmp(url_args[4], "host") == 0 ){
               host = url_args[5];
            } else if( strcmp(url_args[4], "expt") == 0 ){
               expt = url_args[5];
            } else if( strcmp(url_args[4], "histodir") == 0 ){
               histodir = url_args[5];
            } else {
               sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[4]);
               send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
               fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[4]);
            }
         }
         if( narg > 6 ){
            if( strcmp(url_args[6], "host") == 0 ){
               host = url_args[7];
            } else if( strcmp(url_args[6], "expt") == 0 ){
               expt = url_args[7];
            } else if( strcmp(url_args[6], "histodir") == 0 ){
               histodir = url_args[7];
            } else {
               sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[6]);
               send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
               fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[6]);
            }
         }
         if( narg > 8 ){
            if( strcmp(url_args[8], "host") == 0 ){
               host = url_args[9];
            } else if( strcmp(url_args[8], "expt") == 0 ){
               expt = url_args[9];
            } else if( strcmp(url_args[8], "histodir") == 0 ){
               histodir = url_args[9];
            } else {
               sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[8]);
               send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
               fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[8]);
            }
         }
         send_header(fd, APP_JSON); // This header indicates success to the client
         add_sortfile(url_args[3], histodir, host, expt);
      } else {
         if( narg > 4 ){
            if( strcmp(url_args[4], "histodir") == 0 ){
               histodir = url_args[5];
            } else if( strcmp(url_args[4], "configdir") == 0 ){
               configdir = url_args[5];
            } else if( strcmp(url_args[4], "calibrationSource") == 0 ){
               calsrc = url_args[5];
            } else {
               sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[4]);
               send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
               fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[4]);
            }
         }
         if( narg > 6 ){
            if( strcmp(url_args[6], "histodir") == 0 ){
               histodir = url_args[7];
            } else if( strcmp(url_args[6], "configdir") == 0 ){
               configdir = url_args[7];
            } else if( strcmp(url_args[6], "calibrationSource") == 0 ){
               calsrc = url_args[7];
            } else {
               sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[6]);
               send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
               fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[6]);
            }
         }
         if( narg > 8 ){
            if( strcmp(url_args[8], "histodir") == 0 ){
               histodir = url_args[9];
            } else if( strcmp(url_args[8], "configdir") == 0 ){
               configdir = url_args[9];
            } else if( strcmp(url_args[8], "calibrationSource") == 0 ){
               calsrc = url_args[9];
            } else {
               sprintf(tmp,"bad argument:%s in addDataFile\n",url_args[8]);
               send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
               fprintf(stderr,"bad arg:%s in addDataFile\n", url_args[8]);
            }
         }
         send_header(fd, APP_JSON); // This header indicates success to the client
         add_sortfile(url_args[3], histodir, configdir, calsrc);
      }
   } else
   if( strcmp(ptr, "getDatainfo") == 0 ){ /* ---- file info ---------*/
      if( strcmp(url_args[2], "filename") != 0 ){
         sprintf(tmp,"bad argument:%s in getDatainfo\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg:%s in getDatainfo\n", url_args[2]);
      }
      send_file_details(url_args[3], fd);
   } else
   if( strcmp(ptr, "addGlobal") == 0 ){ /* ---------------------- */
      if( narg != 8 ){
         sprintf(tmp,"wrong number of arguments[%d] in addGlobal\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in addGlobal\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"globalname") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in addGlobal\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in addGlobal\n",url_args[2]); return(-1);
      }
      if( strcmp(url_args[4],"globalmin") != 0 ){
         sprintf(tmp,"Bad argument3[%s] in addGlobal\n",url_args[4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg3[%s] in addGlobal\n",url_args[4]); return(-1);
      }
      if( sscanf(url_args[5],"%d", &value) < 1 ){
         sprintf(tmp,"addGlobal: cant Read value: %s\n",url_args[5]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"addGlobal:cantRead val:%s\n",url_args[5]);return(-1);
      }
      if( strcmp(url_args[6],"globalmax") != 0 ){
         sprintf(tmp,"Bad argument4[%s] in addGlobal\n",url_args[6]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg4[%s] in addGlobal\n",url_args[6]); return(-1);
      }
      if( sscanf(url_args[7],"%d", &val2) < 1 ){
         sprintf(tmp,"addGlobal: cant Read value: %s\n",url_args[7]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"addGlobal:cantRead val:%s\n",url_args[7]);return(-1);
      }
      send_header(fd, APP_JSON);
      add_global(configs[0], url_args[3], value, val2);
   } else
   if( strcmp(ptr, "removeGlobal") == 0 ){ /* -------------------- */
      if( narg != 4 ){
         sprintf(tmp,"wrong number of arguments[%d] in removeGlobal\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in removeGlobal\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"globalname") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in removeGlobal\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in removeGlobal\n",url_args[2]);return(-1);
      }
      send_header(fd, APP_JSON);
      remove_global(configs[0], url_args[3]);
   } else
   if( strcmp(ptr, "addCond") == 0 ){ /* ---------------------- */
      if( narg != 10 ){
         sprintf(tmp,"wrong number of arguments[%d] in addCond\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in addCond\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"condname") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in addCond\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in addCond\n", url_args[2]); return(-1);
      }
      if( strcmp(url_args[4],"varname0") != 0 ){
         sprintf(tmp,"Bad argument3[%s] in addCond\n",url_args[4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg3[%s] in addCond\n", url_args[4]); return(-1);
      }
      if( strcmp(url_args[6],"op0") != 0 ){
         sprintf(tmp,"Bad argument4[%s] in addCond\n",url_args[6]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg4[%s] in addCond\n", url_args[6]); return(-1);
      }
      if( strcmp(url_args[8],"value0") != 0 ){
         sprintf(tmp,"Bad argument5[%s] in addCond\n",url_args[8]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg5[%s] in addCond\n", url_args[8]); return(-1);
      }
      if( sscanf(url_args[9],"%d", &value) < 1 ){
         sprintf(tmp,"addCond: cant Read value: %s\n",url_args[9]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"addCond:cantRead value:%s\n",url_args[9]);return(-1);
      }
      send_header(fd, APP_JSON);
      add_cond(configs[0], url_args[3], url_args[5], url_args[7], value);
   } else
   if( strcmp(ptr, "removeCond") == 0 ){ /* -------------------- */
      if( narg != 4 ){
         sprintf(tmp,"wrong number of arguments[%d] in removeCond\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in removeCond\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"condname") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in removeCond\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in removeCond\n",url_args[2]);return(-1);
      }
      send_header(fd, APP_JSON);
      remove_cond(configs[0], url_args[3]);
   } else
   if( strcmp(ptr, "addGate") == 0 ){ /* ---------------------- */
      if( narg < 10 ){
         sprintf(tmp,"wrong number of arguments[%d] in addGate\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in addGate\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"gatename") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in addGate\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in addGate\n", url_args[2]); return(-1);
      }
      add_gate(configs[0], url_args[3]);
      i = 4; while( i < narg ){
         if( strncmp(url_args[i],"varname", 7) != 0 ){
            sprintf(tmp,"Bad argument%d[%s] in addGate\n",i,url_args[2]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"bad arg%d[%s] in addGate\n", i, url_args[i]);
            return(-1);
         }
         if( strncmp(url_args[i+2],"op", 2) != 0 ){
            sprintf(tmp,"Bad argument%d[%s] in addGate\n",i+2,url_args[i+2]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"bad arg%d[%s] in addGate\n", i+2, url_args[i+2]);
            return(-1);
         }
         if( strncmp(url_args[i+4],"value", 5) != 0 ){
            sprintf(tmp,"Bad argument%d[%s] in addGate\n",i+4,url_args[i+4]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"bad arg%d[%s] in addGate\n", i+4, url_args[i+4]);
            return(-1);
         }
         if( sscanf(url_args[i+5],"%d", &value) < 1 ){
            sprintf(tmp,"addGate: cant Read value: %s\n",url_args[i+5]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addGate:cantRead value:%s\n", url_args[i+5]);
            return(-1);
         }
         if( (j = next_condname(configs[0]) ) == -1 ){ return(-1);}
         sprintf(tmp,"Cond%d", j);
         add_cond(configs[0], tmp, url_args[i+1], url_args[i+3], value);
         add_cond_to_gate(configs[0], url_args[3], tmp);
         i += 6;
     }
     send_header(fd, APP_JSON);
   } else
   if( strcmp(ptr, "removeGate") == 0 ){ /* -------------------- */
      if( narg != 4 ){
         sprintf(tmp,"wrong number of arguments[%d] in removeGate\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in removeGate\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"gatename") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in removeGate\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in removeGate\n",url_args[2]);return(-1);
      }
      send_header(fd, APP_JSON);
      remove_gate(configs[0], url_args[3]);
   } else
   if( strcmp(ptr, "applyGate") == 0 ){ /* ---------------------- */
      if( narg != 6 ){
         sprintf(tmp,"wrong number of arguments[%d] in applyGate\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in applyGate\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"gatename") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in applyGate\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in applyGate\n",url_args[2]);return(-1);
      }
      if( strcmp(url_args[4],"histoname") != 0 ){
         sprintf(tmp,"Bad argument3[%s] in applyGate\n",url_args[4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg3[%s] in applyGate\n",url_args[4]);return(-1);
      }
      send_header(fd, APP_JSON);
      apply_gate(configs[0], url_args[3], url_args[5]);
   } else
   if( strcmp(ptr, "unapplyGate") == 0 ){ /* -------------------- */
      if( narg != 6 ){
         sprintf(tmp,"wrong number of arguments[%d] in unapplyGate\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in unapplyGate\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"gatename") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in unapplyGate\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"unapplyGate:bad arg2[%s]\n",url_args[2]);return(-1);
      }
      if( strcmp(url_args[4],"histoname") != 0 ){
         sprintf(tmp,"Bad argument3[%s] in applyGate\n",url_args[4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"unapplyGate:bad arg3[%s]\n",url_args[4]);return(-1);
      }
      send_header(fd, APP_JSON);
      unapply_gate(configs[0], url_args[3], url_args[5]);
   } else
   if( strcmp(ptr, "addHistogram") == 0 ){ /* ---------------------- */
      if( narg < 15 ){
         sprintf(tmp,"wrong number of arguments[%d] in addHistogram\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"too few #arg[%d] in addHistogram\n",narg);return(-1);
      }
      if( strcmp(url_args[2],"name") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in addHistogram\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in addHisto\n",url_args[2]); return(-1);
      }
      if( strcmp(url_args[4],"title") != 0 ){
         sprintf(tmp,"Bad argument3[%s] in addHistogram\n",url_args[4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg3[%s] in addHisto\n",url_args[4]); return(-1);
      }
      if( strcmp(url_args[6],"path") != 0 ){
         sprintf(tmp,"Bad argument4[%s] in addHistogram\n",url_args[6]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg4[%s] in addHisto\n",url_args[6]); return(-1);
      }
      if( strcmp(url_args[8],"xvarname") != 0 ){
         sprintf(tmp,"Bad argument6[%s] in addHistogram\n",url_args[8]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg6[%s] in addHisto\n",url_args[8]); return(-1);
      }
      if( strcmp(url_args[10],"xbins") != 0 ){
         sprintf(tmp,"Bad argument8[%s] in addHistogram\n",url_args[10]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg8[%s] in addHisto\n",url_args[10]); return(-1);
      }
      if( sscanf(url_args[11],"%d", &xbins) < 1 ){
         sprintf(tmp,"addHistogram: cant Read xbins[%s]\n",url_args[11]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"addHisto:cantRd xbins[%s]\n",url_args[11]);return(-1);
      }
      if( strcmp(url_args[12],"xmin") != 0 ){
         sprintf(tmp,"Bad argument10[%s] in addHistogram\n",url_args[12]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg10[%s] in addHisto\n",url_args[12]);return(-1);
      }
      if( sscanf(url_args[13],"%d", &xmin) < 1 ){
         sprintf(tmp,"addHistogram: cant Read xmin[%s]\n",url_args[13]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"addHisto:cantRd xmin[%s]\n",url_args[13]);return(-1);
      }
      if( strcmp(url_args[14],"xmax") != 0 ){
         sprintf(tmp,"Bad argument12[%s] in addHistogram\n",url_args[14]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg12[%s] in addHisto\n",url_args[14]);return(-1);
      }
      if( sscanf(url_args[15],"%d", &xmax) < 1 ){
         sprintf(tmp,"addHistogram: cant Read xmax[%s]\n",url_args[15]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"addHisto:cantRd xmax[%s]\n",url_args[15]);return(-1);
      }
      ybins = 0;
      if( strcmp(url_args[16],"yvarname") == 0 ){
         if( strcmp(url_args[18],"ybins") != 0 ){
            sprintf(tmp,"Bad argument16[%s] in addHistogram\n",url_args[18]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addHisto:bad arg16[%s]\n",url_args[18]);return(-1);
         }
         if( sscanf(url_args[19],"%d", &ybins) < 1 ){
            sprintf(tmp,"addHistogram: cant Read y bins[%s]\n",url_args[19]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addHisto: ?ybins?[%s]\n",url_args[19]);return(-1);
         }
         if( strcmp(url_args[20],"ymin") != 0 ){
            sprintf(tmp,"Bad argument18[%s] in addHistogram\n",url_args[20]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addHisto:bad arg18[%s]\n",url_args[20]);return(-1);
         }
         if( sscanf(url_args[21],"%d", &ymin) < 1 ){
            sprintf(tmp,"addHistogram: cant Read ymin[%s]\n",url_args[21]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addHisto: ?ymin?[%s]\n",url_args[21]);return(-1);
         }
         if( strcmp(url_args[22],"ymax") != 0 ){
            sprintf(tmp,"Bad argument20[%s] in addHistogram\n",url_args[22]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addHisto:bad arg20[%s]\n",url_args[22]);return(-1);
         }
         if( sscanf(url_args[23],"%d", &ymax) < 1 ){
            sprintf(tmp,"addHistogram: cant Read ymax[%s]\n",url_args[23]);
            send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
            fprintf(stderr,"addHisto: ?ymax?[%s]\n",url_args[23]);return(-1);
         }
         add_histo(configs[0], url_args[3], url_args[5], url_args[7], xbins, url_args[9], xmin, xmax, ybins, url_args[17], ymin, ymax);
         i=24;
      } else {
         add_histo(configs[0], url_args[3], url_args[5], url_args[7], xbins, url_args[9], xmin, xmax, 0, "", 0, 0);
         i=16;
      }
      for(;i<narg; i+=2){
         if( sscanf(url_args[i], "gate%d", &j) < 1 ){
             sprintf(tmp,"addHistogram: bad arg%d[%s]\n",i,url_args[i]);
             send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
             fprintf(stderr,"addHisto: bad arg%d[%s]\n", i, url_args[i]);
             return(-1);
         }
         apply_gate(configs[0], url_args[3], url_args[i+1]);
      }
      send_header(fd, APP_JSON);
   } else
   if( strcmp(ptr, "removeHistogram") == 0 ){ /* -------------------- */
      if( narg != 4 ){
         sprintf(tmp,"wrong number of arguments[%d] in removeHistogram\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in removeHisto\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"histoname") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in removeHistogram\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in removeHist\n",url_args[2]);return(-1);
      }
      send_header(fd, APP_JSON);
      remove_histo(configs[0], url_args[3]);
   } else
   if( strcmp(ptr, "sumHistos") == 0 ){ /* -------------------- */
      queue_sum_histos(configs[0], narg, url_args, fd);
   } else
   if( strcmp(ptr, "setCalibration") == 0 ){ /* -------------------- */
      set_calibration(configs[0], narg, url_args, fd);
   } else
   if( strcmp(ptr, "setPileupCorrection") == 0 ){ /* -------------------- */
      set_pileup_correction(configs[0], narg, url_args, fd);
   } else
   if( strcmp(ptr, "setCrosstalkCorrection") == 0 ){ /* -------------------- */
      set_crosstalk_correction(configs[0], narg, url_args, fd);
   } else
   if( strcmp(ptr, "setDataDirectory") == 0 ){ /* -------------------- */
      send_header(fd, APP_JSON);
      set_directory(configs[0], "Data", url_args[3]);
   } else
   if( strcmp(ptr, "setHistoDirectory") == 0 ){ /* -------------------- */
      send_header(fd, APP_JSON);
      set_directory(configs[0], "Histo", url_args[3]);
   } else
   if( strcmp(ptr, "setConfigDirectory") == 0 ){ /* -------------------- */
      send_header(fd, APP_JSON);
      set_directory(configs[0], "Config", url_args[3]);
   } else
   if( strcmp(ptr, "saveConfig") == 0 ){ /* -------------------- */
      if( narg != 4 ){
         sprintf(tmp,"wrong number of arguments[%d] in saveConfig\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in saveConfig\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"filename") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in saveConfig\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in saveConfig\n",url_args[2]);return(-1);
      }
      send_header(fd, APP_JSON);
      save_config(configs[0], url_args[3], 1); // 1 => overwrite
   } else
   if( strcmp(ptr, "loadConfig") == 0 ){ /* -------------------- */
      if( narg != 4 ){
         sprintf(tmp,"wrong number of arguments[%d] in loadConfig\n",narg);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"wrong #arg[%d] in loadConfig\n", narg); return(-1);
      }
      if( strcmp(url_args[2],"filename") != 0 ){
         sprintf(tmp,"Bad argument2[%s] in loadConfig\n",url_args[2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"bad arg2[%s] in loadConfig\n",url_args[2]);return(-1);
      }
      send_header(fd, APP_JSON);
      load_config(configs[0], url_args[3], NULL);
   } else
   if( strcmp(ptr, "viewConfig") == 0 ){ /* -------------------- */
      if( web_fp == NULL ){
         if( (web_fp=fdopen(fd,"r+")) == NULL ){
            send_http_error_response(fd, STATUS_CODE_422,(char*)"viewConfig can't fdopen web fd");
            fprintf(stderr,"viewConfig can't fdopen web fd\n"); return(-1);
         }
      }
      if( strcmp(url_args[2],"filename") == 0 ){
         i = strlen(url_args[3]);
         if( strncmp(url_args[3]+i-5, ".json", 5) == 0 ){
            if( (cfg=add_config(url_args[3])) == NULL ){
               fprintf(stderr,"viewConfig: cant create config:%s\n", url_args[3]);
            }
            if( load_config(cfg, url_args[3], NULL) != 0 ){
               fprintf(stderr,"viewConfig: cant load config:%s\n", url_args[3]);
            }
         } else
         if( (cfg=read_histofile(url_args[3], 1)) == NULL ){
            send_http_error_response(fd, STATUS_CODE_404,(char*)"viewConfig can't read requested filename");
            fprintf(stderr,"viewConfig can't read requested filename, %s\n",url_args[3]);
           return(-1);
         }
         // web_fp did not work under certain conditions (see below)
         sprintf(tmp,"/tmp/tmp.json");
         if( (tmp_fp = fopen(tmp,"w+")) == NULL ){ // create if needed, truncate to zero
            send_http_error_response(fd, STATUS_CODE_500,(char*)"viewConfig can't open tmp file to write");
            fprintf(stderr,"can't open tmp file:%s to write\n", tmp);
           //return(-1);
         }
         send_header(fd, APP_JSON);
         write_config(cfg, tmp_fp); fseek(tmp_fp, 0, SEEK_SET);
         while( fgets(tmp, 128, tmp_fp) != NULL ){
            put_line(fd, tmp, strlen(tmp) );
         }
         fclose(tmp_fp);
         remove_config(cfg);
      } else {
         //write_config(configs[0], web_fp); fflush(web_fp); // empty if online??
         printf("ViewCurrentConfig\n");
         // hack workaround for now ...
         sprintf(tmp,"/tmp/tmp.json");
         if( (tmp_fp = fopen(tmp,"w+")) == NULL ){ // create if needed, truncate to zero
            send_http_error_response(fd, STATUS_CODE_500,(char*)"viewConfig can't open tmp file to write");
            fprintf(stderr,"can't open tmp file:%s to write\n", tmp);
           //return(-1);
         }
         send_header(fd, APP_JSON);
         write_config(configs[0], tmp_fp); fseek(tmp_fp, 0, SEEK_SET);
         while( fgets(tmp, 128, tmp_fp) != NULL ){
            put_line(fd, tmp, strlen(tmp) );
         }
         fclose(tmp_fp);
      }
   } else
   if( strcmp(ptr, "callspechandler") == 0 ){
      if( strcmp(url_args[2],"filename") == 0 ){
         send_spectrum( (narg-2)/2, url_args, url_args[3], fd);
      } else {
         send_spectrum( (narg-2)/2, url_args, NULL, fd);
      }
   } else
   if( strcmp(ptr, "callbinaryspechandler") == 0 ){
      if( strcmp(url_args[2],"filename") == 0 ){
         send_binary_spectrum( ((int)(narg-2)/2), url_args, url_args[3], fd);
      } else {
         send_binary_spectrum( ((int)(narg-2)/2), url_args, NULL, fd);
      }
   } else {
         sprintf(tmp,"Unknown Command: %s\n",ptr);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
      fprintf(stderr,"Unknown Command:%s\n", ptr);
   }
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////            Config  I/O            ///////////////////////
///////////////////////////////////////////////////////////////////////////
Config *configs[MAX_CONFIGS];

int write_config(Config *cfg, FILE *fp)
{
   Histogram *histo;  Cal_coeff *calib;  Global *global;
   Sortvar *var;      Cond *cond;        Gate *gate;
   int i, j, first;
   char tmp[64];

   fprintf(fp,"{\n   \"Analyzer\" : [\n");
   fprintf(fp,"      {\"Variables\" : [\n");
   for(i=0; i<cfg->nsortvar; i++){ var = &cfg->varlist[i];
      sprintf(tmp, "\"%s\"", var->name );
      fprintf(fp,"%9s{ \"name\" : %-16s,", "", tmp );
      fprintf(fp,"   \"title\" : \"%s\"}", var->title  );
      fprintf(fp,"%s", ( i<cfg->nsortvar-1 ) ? ",\n" : "\n" );
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Gates\" : [\n");
   for(i=0; i<cfg->ngates; i++){ gate = cfg->gatelist[i];
      fprintf(fp,"%9s{ \"name\" : \"%s\",\n", "", gate->name );
      fprintf(fp,"%12s\"gateCondition\" : [\n", "" );
      for(j=0; j<gate->nconds; j++){ cond = gate->conds[j];
         fprintf(fp,"%15s{\"indexID\" : %d, \"Variable\" : \"", "", j );
         fprintf(fp,"%s\" , \"Logic\" : \"%s\" , \"Value\" : %d",
                 cond->var->name, GATEOP(cond->op), cond->value );
         fprintf(fp, "%s", ( j<gate->nconds-1 ) ? "},\n" : "}\n" );
      }
      fprintf(fp,"            ]\n");
      fprintf(fp, "         %s", ( i<cfg->ngates-1 ) ? "},\n" : "}\n" );
   }
   fprintf(fp,"      ]},\n"); first = 1;
   fprintf(fp,"      {\"Histograms\" : [\n");
   for(i=0; i<cfg->nhistos; i++){ histo = cfg->histo_list[i];
      if( ! histo->user ){ continue; } // only save user histos in config
      if( first ){ first = 0; } else { fprintf(fp,",\n"); }
      fprintf(fp,"%9s{\"name\" : \"%s\",\n", "", histo->handle );
      fprintf(fp,"%12s\"path\" : \"%s\",\n", "", histo->path );
      fprintf(fp,"%12s\"Xvariable\" : \"%s\",\n", "", histo->xvar->name );
      fprintf(fp,"%12s\"Xmin\" : %d,", "", 0 );                              // This needs its own variable
      fprintf(fp,"%12s\"Xmax\" : %d,", "", histo->xbins );                   // This needs its own variable
      fprintf(fp,"%12s\"Xbins\" : %d", "", histo->xbins );
       if( histo->ybins != 0 ){
         fprintf(fp,",\n%12s\"Yvariable\" : \"%s\",\n","", histo->yvar->name );
         fprintf(fp,"%12s\"Ymin\" : %d,", "", 0 );                              // This needs its own variable
         fprintf(fp,"%12s\"Ymax\" : %d,", "", histo->ybins );                   // This needs its own variable
         fprintf(fp,"%12s\"Ybins\" : %d", "", histo->ybins );
       }
      fprintf(fp,",\n%12s\"histogramCondition\" : [\n", "");
      for(j=0; j<histo->num_gates; j++){
         fprintf(fp,"%15s{\"indexID\" : %d , \"Gate\" : \"%s\" ", "",
                                            j, histo->gate_names[j] );
         fprintf(fp, "%s", (j<histo->num_gates-1 ) ? "},\n" : "}\n" );
      }
      fprintf(fp,"]\n%9s}", "");
   }
   fprintf(fp,"\n      ]},\n");
   fprintf(fp,"      {\"Globals\" : [\n");
   for(i=0; i<cfg->nglobal; i++){ global = cfg->globals[i];
      fprintf(fp,"%9s{\"name\" : \"%s\" , \"min\" : %d , \"max\" : %d ",
              "", global->name, global->min, global->max );
      fprintf(fp, "%s", ( i<cfg->nglobal-1 ) ? "},\n" : "}\n" );
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Calibrations\" : [\n");
   for(i=0; i<cfg->ncal; i++){ calib = cfg->calib[i];
     fprintf(fp,"%9s{\"name\" : \"%s\" , \"address\" : %d , \"datatype\" : %d , \"offset\" : %f , \"gain\" : %f , \"quad\" : %e ", "", calib->name, calib->address, calib->datatype, calib->offset, calib->gain, calib->quad );
     if(strncmp(calib->name,"GRG",3)==0){
       if(calib->pileupk1[0] != -1){
         fprintf(fp,", \"pileupk1\" : [ %f , %f , %e , %e , %e , %e , %e ]",calib->pileupk1[0],calib->pileupk1[1],calib->pileupk1[2],calib->pileupk1[3],calib->pileupk1[4],calib->pileupk1[5],calib->pileupk1[6]);
         fprintf(fp,", \"pileupk2\" : [ %f , %f , %e , %e , %e , %e , %e ]",calib->pileupk2[0],calib->pileupk2[1],calib->pileupk2[2],calib->pileupk2[3],calib->pileupk2[4],calib->pileupk2[5],calib->pileupk2[6]);
         fprintf(fp,", \"pileupE1\" : [ %f , %f , %e , %e , %e , %e , %e ]",calib->pileupE1[0],calib->pileupE1[1],calib->pileupE1[2],calib->pileupE1[3],calib->pileupE1[4],calib->pileupE1[5],calib->pileupE1[6]);
       }else{
         fprintf(fp,", \"pileupk1\" : [ %d , %d , %d , %d , %d , %d , %d ]",1,0,0,0,0,0,0);
         fprintf(fp,", \"pileupk2\" : [ %d , %d , %d , %d , %d , %d , %d ]",1,0,0,0,0,0,0);
         fprintf(fp,", \"pileupE1\" : [ %d , %d , %d , %d , %d , %d , %d ]",0,0,0,0,0,0,0);
       }
       if(calib->crosstalk0[0] != -1){
         fprintf(fp,", \"crosstalk0\" : [ %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f ]",calib->crosstalk0[0],calib->crosstalk0[1],calib->crosstalk0[2],calib->crosstalk0[3],calib->crosstalk0[4],calib->crosstalk0[5],calib->crosstalk0[6],calib->crosstalk0[7],calib->crosstalk0[8],calib->crosstalk0[9],calib->crosstalk0[10],calib->crosstalk0[11],calib->crosstalk0[12],calib->crosstalk0[13],calib->crosstalk0[14],calib->crosstalk0[15]);
         fprintf(fp,", \"crosstalk1\" : [ %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f ]",calib->crosstalk1[0],calib->crosstalk1[1],calib->crosstalk1[2],calib->crosstalk1[3],calib->crosstalk1[4],calib->crosstalk1[5],calib->crosstalk1[6],calib->crosstalk1[7],calib->crosstalk1[8],calib->crosstalk1[9],calib->crosstalk1[10],calib->crosstalk1[11],calib->crosstalk1[12],calib->crosstalk1[13],calib->crosstalk1[14],calib->crosstalk1[15]);
         fprintf(fp,", \"crosstalk2\" : [ %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f , %f ]",calib->crosstalk2[0],calib->crosstalk2[1],calib->crosstalk2[2],calib->crosstalk2[3],calib->crosstalk2[4],calib->crosstalk2[5],calib->crosstalk2[6],calib->crosstalk2[7],calib->crosstalk2[8],calib->crosstalk2[9],calib->crosstalk2[10],calib->crosstalk2[11],calib->crosstalk2[12],calib->crosstalk2[13],calib->crosstalk2[14],calib->crosstalk2[15]);
       }else{
         fprintf(fp,", \"crosstalk0\" : [ %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d ]",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
         fprintf(fp,", \"crosstalk1\" : [ %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d ]",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
         fprintf(fp,", \"crosstalk2\" : [ %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d , %d ]",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
       }
     }
      fprintf(fp, "%s", ( i<cfg->ncal-1 ) ? "},\n" : "}\n" );
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Directories\" : [\n");
   {
      fprintf(fp,"%9s{\"name\" : \"Data\", ", "");
      fprintf(fp,"\"Path\" : \"%s\"},\n", cfg->data_dir);
      fprintf(fp,"%9s{\"name\" : \"Histo\", ", "");
      fprintf(fp,"\"Path\" : \"%s\"},\n", cfg->histo_dir);
      fprintf(fp,"%9s{\"name\" : \"Config\", ", "");
      fprintf(fp,"\"Path\" : \"%s\"}\n", cfg->config_dir); // NO COMMA
   }
   fprintf(fp,"      ]},\n");
   fprintf(fp,"      {\"Midas\" : [\n");
   {
      fprintf(fp,"%9s{\"name\" : \"Title\", ", "");
      fprintf(fp,"\"Value\" : \"%s\"},\n", cfg->midas_title);
      fprintf(fp,"%9s{\"name\" : \"StartTime\", ", "");
      fprintf(fp,"\"Value\" : \"%d\"},\n", cfg->midas_start_time);
      fprintf(fp,"%9s{\"name\" : \"Duration\", ", "");
      fprintf(fp,"\"Value\" : \"%d\"}\n", cfg->midas_runtime);  // NO COMMA
   }
   fprintf(fp,"      ]}\n"); // NO COMMA AFTER FINAL TABLE
   fprintf(fp,"   ]\n}\n"); // Analyser,File
   return(0);
}

static char config_data[1024*1024];
int load_config(Config *cfg, char *filename, char *buffer)
{
   int i,j, len, value, val2, val3, val4, val5, val6, address, type, instring;
   char *ptr, *name, *valstr, *title, *path, *var, *var2, op[8], tmp[80];
   float gain, offset, quad;
   float puk1[7], puk2[7], puE1[7];
   float ct0[16], ct1[16], ct2[16];
   // Initialize values to defaults
   // Values of -1 are ignored by edit_calibration - use this for all channels that are not HPGe to avoid bloating the size of the config
   float puk_reset[7]={1,0,0,0,0,0,0}, puE1_reset[7]={0,0,0,0,0,0,0}, pu_ignore[7]={-1,-1,-1,-1,-1,-1,-1};
   float ct_reset[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, ct_ignore[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
   Histogram *histo;
   Config *tmp_cfg;
   Cond *cond;
   FILE *fp;

   if( filename != NULL ){ len = strlen(filename);
      if( strncmp(filename+len-4, ".tar", 4) == 0 ){ // HISTOGRAM FILE
         if( (tmp_cfg=read_histofile(filename, 1)) == NULL ){ return(-1); }
         copy_config(tmp_cfg, cfg);
         remove_config(tmp_cfg);
         return(0);
      }
      if( (fp=fopen(filename,"r")) == NULL ){
         fprintf(stderr,"load_config: cant open %s to read\n", filename);
         return(-1);
      }
      fprintf(stderr,"load_config: reading %s\n", filename);
      instring = len = 0;// read line by line, copy to conf_data, skip space
      while( fgets(tmp, 80, fp) != NULL ){ // tmp always null-terminated
         for(i=0; i<strlen(tmp); i++){ // DO NOT SKIP SPACE WITHIN STRINGS
            if( tmp[i] == '"' ){ instring = 1-instring; }
            if( !isspace(tmp[i]) || instring ){ config_data[len++] = tmp[i]; }
         }
      }
      fclose(fp);
   } else if( buffer != NULL ){
      instring = len = 0; // read line by line, copy to conf_data, skip space
      for(i=0; i<strlen(buffer); i++){    // DO NOT SKIP SPACE WITHIN STRINGS
         if( buffer[i] == '"' ){ instring = 1-instring; }
         if( !isspace(buffer[i])||instring ){ config_data[len++]=buffer[i]; }
      }
   } else {
      fprintf(stderr,"load_config: no file or buffer specified\n");
   }
   clear_config(cfg); init_user_config(cfg); // setup hardcoded variables
   ptr=config_data;
   if( strncmp(ptr,"{\"Analyzer\":[", 13) != 0 ){
      fprintf(stderr,"load_config: err1 byte %ld\n", ptr-config_data);
       return(-1);
   } ptr += 13;
   if( strncmp(ptr,"{\"Variables\":[", 14) != 0 ){
      fprintf(stderr,"load_config: err2 byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 14;
   while( 1 ){ // variables
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: err3 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"title\":\"", 10) != 0 ){
         fprintf(stderr,"load_config: err4 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 10;
      title = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      // VARIABLES CURRENTLY HARDCODED
      //  - this section of config is only used for passing to viewer
      //cfg->lock=1; add_variable(cfg, name, title);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Gates\":[", 10) != 0 ){
      fprintf(stderr,"load_config: err5 byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 10;
   while( 1 ){ // Gates
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: err6 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=1; add_gate(cfg, name); cfg->lock=0;
      if( strncmp(ptr,",\"gateCondition\":[", 18) != 0 ){
         fprintf(stderr,"load_config: err7 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 18;
      while( 1 ){ // Gate-Conditions
         if( strncmp(ptr,"{\"indexID\":", 11) != 0 ){
            fprintf(stderr,"load_config: err8 byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 11; while( *ptr != ',' ){ ++ptr; /* skip index id */ }
         if( strncmp(ptr,",\"Variable\":\"", 13) != 0 ){
            fprintf(stderr,"load_config: err8b byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 13; var = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
         if( strncmp(ptr,",\"Logic\":\"", 10) != 0 ){
            fprintf(stderr,"load_config: err8c byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 10;
         if(       strncmp(ptr,"GE",2) == 0 ){ sprintf(op,">=");
         } else if(strncmp(ptr,"GT",2) == 0 ){ sprintf(op,">");
         } else if(strncmp(ptr,"LE",2) == 0 ){ sprintf(op,"<=");
         } else if(strncmp(ptr,"LT",2) == 0 ){ sprintf(op,"<");
         } else if(strncmp(ptr,"EQ",2) == 0 ){ sprintf(op,"=");
         } else if(strncmp(ptr,"RA",2) == 0 ){ sprintf(op,"RA");
         } else {
            fprintf(stderr,"load_config:err9 byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 3;
         if( strncmp(ptr,",\"Value\":", 9) != 0 ){
            fprintf(stderr,"load_config: err9b byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 9; valstr=ptr; while( isdigit(*ptr) ){ ++ptr; } *ptr++ = 0;
         if( sscanf( valstr, "%d", &value) < 1 ){
            fprintf(stderr,"load_config:errA byte %ld\n", ptr-config_data);
            return(-1);
         }
         if( (i = next_condname(cfg) ) == -1 ){ return(-1);}
         cfg->lock=1; sprintf(tmp,"Cond%d", i);
         add_cond(cfg, tmp, var, op, value);
         add_cond_to_gate(cfg, name, tmp); cfg->lock=0;
         if( *ptr++ == ',' ){ continue; }
         ++ptr; break; // skip ']'
      }
      if( *ptr++ == ',' ){ continue; }
      ptr += 2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Histograms\":[", 15) != 0 ){
      fprintf(stderr,"load_config: errB byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 15;
   while( 1 ){ // Histograms
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"",9) != 0 ){
         fprintf(stderr,"load_config: errC byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      //if( strncmp(ptr,"{\"title\":\"",10) != 0 ){
      //   fprintf(stderr,"load_config: errC byte %ld\n", ptr-config_data);
      //   return(-1);
      //} ptr += 10;
      //title = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      title = name;
      if( strncmp(ptr,",\"path\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errD byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      path = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Xvariable\":\"", 14) != 0 ){
         fprintf(stderr,"load_config: errE byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 14;
      var = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Xmin\":", 8) != 0 ){
         fprintf(stderr,"load_config: errFa byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 8;
      valstr = ptr; while( isdigit(*ptr) ){ ++ptr; } tmp[0]=*ptr; *ptr++=0;
      if( sscanf( valstr, "%d", &val2) < 1 ){
         fprintf(stderr,"load_config:errFb byte %ld\n", ptr-config_data);
            return(-1);
      }
      if( strncmp(ptr,"\"Xmax\":", 7) != 0 ){
         fprintf(stderr,"load_config: errFc byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 7;
      valstr = ptr; while( isdigit(*ptr) ){ ++ptr; } tmp[0]=*ptr; *ptr++=0;
      if( sscanf( valstr, "%d", &val3) < 1 ){
         fprintf(stderr,"load_config:errFd byte %ld\n", ptr-config_data);
            return(-1);
      }
      if( strncmp(ptr,"\"Xbins\":", 8) != 0 ){
         fprintf(stderr,"load_config: errFe byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 8;
      valstr = ptr; while( isdigit(*ptr) ){ ++ptr; } tmp[0]=*ptr; *ptr++=0;
      if( sscanf( valstr, "%d", &value) < 1 ){
         fprintf(stderr,"load_config:errFf byte %ld\n", ptr-config_data);
            return(-1);
      }
      if( strncmp(ptr,"\"Yvariable\":\"", 13) != 0 ){
         cfg->lock=1;
         add_histo(cfg, name, title, path, value, var, val2, val3, 0, "",0,0);
         cfg->lock=0;
      } else {
         ptr += 13;
         var2 = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
         if( strncmp(ptr,",\"Ymin\":", 8) != 0 ){
            fprintf(stderr,"load_config: errG byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 8;
         valstr = ptr; while( isdigit(*ptr) ){++ptr;} tmp[0]=*ptr;*ptr++=0;
         if( sscanf( valstr, "%d", &val5) < 1 ){
            fprintf(stderr,"load_config:errH byte %ld\n", ptr-config_data);
               return(-1);
         }
         if( strncmp(ptr,"\"Ymax\":", 7) != 0 ){
            fprintf(stderr,"load_config: errHa byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 7;
         valstr = ptr; while( isdigit(*ptr) ){++ptr;} tmp[0]=*ptr;*ptr++=0;
         if( sscanf( valstr, "%d", &val6) < 1 ){
            fprintf(stderr,"load_config:errHb byte %ld\n", ptr-config_data);
               return(-1);
         }
         if( strncmp(ptr,"\"Ybins\":", 8) != 0 ){
            fprintf(stderr,"load_config: errHc byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 8;
         valstr = ptr; while( isdigit(*ptr) ){++ptr;} tmp[0]=*ptr;*ptr++=0;
         if( sscanf( valstr, "%d", &val4) < 1 ){
            fprintf(stderr,"load_config:errI byte %ld\n", ptr-config_data);
               return(-1);
         }
         cfg->lock=1;
         add_histo(cfg, name, title, path, value, var, val2, val3, val4, var2,val5,val6); cfg->lock=0;
      }
      if( strncmp(ptr,"\"histogramCondition\":[", 22) != 0 ){
         fprintf(stderr,"load_config: errJ byte %ld\n", ptr-config_data);
            return(-1);
      } ptr += 22;
      while(1){  // Histo gates
         if( strncmp(ptr,"]}", 2) == 0 ){ ptr+=2; break; } // empty list
         if( strncmp(ptr,"{\"indexID\":", 11) != 0 ){
            fprintf(stderr,"load_config: errK byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 11; while( isdigit(*ptr) || *ptr=='-' ){++ptr;}
         if( strncmp(ptr,",\"Gate\":\"", 9) != 0 ){
            fprintf(stderr,"load_config: errKa byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 9;
         valstr = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
         ++ptr; // skip '}' - end of single condition
         cfg->lock=1; apply_gate(cfg, name, valstr); cfg->lock=0;
         if( *ptr == ',' ){ ++ptr; continue; } else { ptr +=2; break;}// ']}'
      }
      if( strncmp(ptr,",{\"", 3) == 0 ){ ++ptr; continue; }
      ptr += 3; break; // skip closing ]},
   }
   if( strncmp(ptr,"{\"Globals\":[", 12) != 0 ){
      fprintf(stderr,"load_config: errL byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 12;
   while( 1 ){ // Globals
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errM byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"min\":",7) != 0 ){
         fprintf(stderr,"load_config: errN byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 7;
      valstr = ptr; while( isdigit(*ptr) || *ptr=='-' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &value) < 1 ){
         fprintf(stderr,"load_config:errO byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"max\":",6) != 0 ){
         fprintf(stderr,"load_config: errP byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 6;
      valstr = ptr; while( isdigit(*ptr) || *ptr=='-' ){++ptr;} *ptr++ = 0;
      // last char written over with 0 was '}'
      if( sscanf( valstr, "%d", &val2) < 1 ){
         fprintf(stderr,"load_config:errQ byte %ld\n", ptr-config_data);
         return(-1);
      }
      cfg->lock=1; add_global(cfg, name, value, val2); cfg->lock=0;
      if( *ptr++ == ',' ){ continue; } // have skipped ']' if not
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Calibrations\":[", 17) != 0 ){
      fprintf(stderr,"load_config: errR byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 17;
   while( 1 ){ // Calibrations
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; fprintf(stdout,"Calibrations section empty so breaking here.\n"); break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errS byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"address\":",11) != 0 ){
         fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &address) < 1 ){
         fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"datatype\":",11) != 0 ){
         fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &type) < 1 ){
        fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"offset\":",9) != 0 ){
        fprintf(stderr,"load_config: errT byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 9; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &offset) < 1 ){
        fprintf(stderr,"load_config:errU byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"gain\":",7) != 0 ){
        fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 7; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &gain) < 1 ){
        fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"quad\":",7) != 0 ){
        fprintf(stderr,"load_config: errY byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 7; valstr = ptr;
      while(isdigit(*ptr)||*ptr=='.'||*ptr=='-'||*ptr=='+'||*ptr=='e'){++ptr;}
      *ptr++ = 0; // last char written over with 0 was '}'
      if( sscanf( valstr, "%f", &quad) < 1 ){
        fprintf(stderr,"load_config:errZ byte %ld\n", ptr-config_data);
        return(-1);
      }
      // The pileup correction parameters were introduced in Feb 2025.
      // Config files before this date will not have pileup corrections, and after this they are optional
      if( strncmp(ptr,"\"pileupk1\":",11) == 0 ){
        ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e], ", &puk1[0],&puk1[1],&puk1[2],&puk1[3],&puk1[4],&puk1[5],&puk1[6]) != 7 ){
          fprintf(stderr,"load_config:errPUA byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"pileupk2\":",11) != 0 ){
          fprintf(stderr,"load_config:errPUB byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e], ", &puk2[0],&puk2[1],&puk2[2],&puk2[3],&puk2[4],&puk2[5],&puk2[6]) != 7 ){
          fprintf(stderr,"load_config:errPUC byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"pileupE1\":",11) != 0 ){
          fprintf(stderr,"load_config:errPUD byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e]", &puE1[0],&puE1[1],&puE1[2],&puE1[3],&puE1[4],&puE1[5],&puE1[6]) != 7 ){
          fprintf(stderr,"load_config:errPUE byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!=']'){++ptr;}
        ptr+=2;
      }else if(strncmp(name,"GRG",3)==0){ // Only process pileup and crosstalk for HPGe
           memcpy(puk1,puk_reset, 7 * sizeof(float));
           memcpy(puk2,puk_reset, 7 * sizeof(float));
           memcpy(puE1,puE1_reset, 7 * sizeof(float));
      }else{
            memcpy(puk1,pu_ignore, 7 * sizeof(float));
            memcpy(puk2,pu_ignore, 7 * sizeof(float));
            memcpy(puE1,pu_ignore, 7 * sizeof(float));
      }
      // The crosstalk correction parameters were introduced in June 2025.
      // Config files before this date will not have crosstalk corrections, and after this they are optional
      if( strncmp(ptr,"\"crosstalk0\":",13) == 0 ){
        ptr += 13; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f], ", &ct0[0],&ct0[1],&ct0[2],&ct0[3],&ct0[4],&ct0[5],&ct0[6],&ct0[7],&ct0[8],&ct0[9],&ct0[10],&ct0[11],&ct0[12],&ct0[13],&ct0[14],&ct0[15]) != 16 ){
          fprintf(stderr,"load_config:errCTA byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"crosstalk1\":",13) != 0 ){
          fprintf(stderr,"load_config:errCTB byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 13; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f], ", &ct1[0],&ct1[1],&ct1[2],&ct1[3],&ct1[4],&ct1[5],&ct1[6],&ct1[7],&ct1[8],&ct1[9],&ct1[10],&ct1[11],&ct1[12],&ct1[13],&ct1[14],&ct1[15]) != 16 ){
          fprintf(stderr,"load_config:errCTC byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"crosstalk2\":",13) != 0 ){
          fprintf(stderr,"load_config:errCTD byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 13; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f]", &ct2[0],&ct2[1],&ct2[2],&ct2[3],&ct2[4],&ct2[5],&ct2[6],&ct2[7],&ct2[8],&ct2[9],&ct2[10],&ct2[11],&ct2[12],&ct2[13],&ct2[14],&ct2[15]) != 16 ){
          fprintf(stderr,"load_config:errCTE byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!=']'){++ptr;}
        ptr+=2;
      }else if(strncmp(name,"GRG",3)==0){ // Only process pileup and crosstalk for HPGe
            memcpy(ct0,ct_reset, 16 * sizeof(float));
            memcpy(ct1,ct_reset, 16 * sizeof(float));
            memcpy(ct2,ct_reset, 16 * sizeof(float));
      }else{
            memcpy(ct0,ct_ignore, 16 * sizeof(float));
            memcpy(ct1,ct_ignore, 16 * sizeof(float));
            memcpy(ct2,ct_ignore, 16 * sizeof(float));
      }
      cfg->lock=1; edit_calibration(cfg,name,offset,gain,quad,puk1,puk2,puE1,ct0,ct1,ct2,address,type,1); cfg->lock=0;
      if( *ptr++ == ',' ){ continue; } // have skipped ']' if not
      ptr+=2; break; // skip '},'
    }
   if( strncmp(ptr,"{\"Directories\":[", 16) != 0 ){
      fprintf(stderr,"load_config: errZA byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 16;
   while(1){
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZB byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Path\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZC byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      path = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=1; set_directory(cfg, name, path);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Midas\":[", 10) != 0 ){
      fprintf(stderr,"load_config: errZD byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 10;
   while(1){
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZE byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Value\":\"", 10) != 0 ){
         fprintf(stderr,"load_config: errZF byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 10;
      valstr = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=1; set_midas_param(cfg, name, valstr);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ++ptr; break; // skip '}'
   }
   if( strncmp(ptr,"]}", 2) != 0 || ptr+2-config_data != len ){
      fprintf(stderr,"load_config: errR near %ld\n", ptr-config_data);
      return(-1);
   }
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////          Config Commands          ///////////////////////
///////////////////////////////////////////////////////////////////////////
// config changes during sorting ... (sort probably running most of time)
// - sort will need to grab config and copy so any changes do not affect it
//


// init_config() call on server startup ...
//   get most recent config and initialise all variables, gates, histos
// *In case no recent config exists, a default config is setup first*
//  (this is usually immediately overwritten by the recent config)
//  (could change this to only setup the default, if recent fails)
int init_config()
{
  char hostname[32];
  int tmp;
   Config *cfg = configs[0];
   if( cfg == NULL ){ // not yet alloc'd live set
      if( (cfg=configs[0]=add_config("live")) == NULL ){ return(-1); }
      if( (configs[1]=add_config("sort")) == NULL ){ return(-1); }
      configs[0]->type = configs[1]->type = MEM_CONFIG;
   }
   init_default_config(cfg);  // populate default "test" config during testing
   load_config(cfg, DEFAULT_CONFIG, NULL); // attempt to load, ignore any error
   clear_calibrations(cfg); // Clear the calibrations to default values following server restart
   tmp = gethostname(hostname, 32);
   strtok(hostname, ".");
   fprintf(stdout,"Initial setup complete :-)\n\n");
   fprintf(stdout,"Now connect to grif-replay using a web browser at the following URL:\nhttps://griffincollaboration.github.io/SpectrumViewer/analyzerInterface.html?backend=%s&port=9093\n",hostname);
   return(0);
}

int clear_config(Config *cfg)
{
   time_t current_time = time(NULL);
   Histogram *histo;
   int i;
   //if( strcmp(cfg->name, "live") == 0 || strcmp(cfg->name, "sort") == 0 ){
   //   fprintf(stderr,"refusing to delete config:%s\n", cfg->name);
   //   return(-1);
   //}
   for(i=0; i<MAX_HISTOGRAMS; i++){ histo = &cfg->histo_array[i];
      if( histo->data != NULL ){ free(histo->data); histo->data = NULL; }
   }
   memset(cfg, 0, sizeof(Config));      // delete any current vars, gates etc.
   for(i=0; i<MAX_CALIB;     i++){cfg->calib[i]    = &cfg->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS;   i++){cfg->globals[i]  = &cfg->global_array[i]; }
   for(i=0; i<MAX_CONDS;     i++){cfg->condlist[i] = &cfg->cond_array[i]; }
   for(i=0; i<MAX_GATES;     i++){cfg->gatelist[i] = &cfg->gate_array[i]; }
   // used_sortvars, and user histos start empty and are filled randomly
   cfg->mtime = current_time;
   return(0);
}

int clear_calibrations(Config *cfg)
{
  float offset=0, gain=1, quad=0;
  int i, address=-1, datatype=-1;
  // Initialize values to defaults
  // Values of -1 are ignored by edit_calibration - use this for all channels that are not HPGe to avoid bloating the size of the config
  float puk_reset[7]={1,0,0,0,0,0,0}, puE1_reset[7]={0,0,0,0,0,0,0}, pu_ignore[7]={-1,-1,-1,-1,-1,-1,-1};
  float ct_reset[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, ct_ignore[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  // delete any calibration values
  for(i=0; i<cfg->ncal;     i++){
    if(strncmp(cfg->calib[i]->name,"GRG",3)==0){
      edit_calibration(cfg, cfg->calib[i]->name, offset, gain, quad, puk_reset, puk_reset, puE1_reset, ct_reset, ct_reset, ct_reset, address, datatype, 1);
    }else{
      edit_calibration(cfg, cfg->calib[i]->name, offset, gain, quad, pu_ignore, pu_ignore, pu_ignore, ct_ignore, ct_ignore, ct_ignore, address, datatype, 1);
    }
  }
  return(0);
}

int copy_config(Config *src, Config *dst)
{
   Histogram *histo, *srchist;  Gate *gate;  Cond *cond;  Global *global;
   Sortvar *srcvar, *dstvar;
   char *tmp, *tmp2, *ptr;
   long *tmp3;
   int i, j;

   src->lock = 1;
   delete_histograms(dst);              // delete current histograms
   memset(dst, 0, sizeof(Config));      // delete any current vars, gates etc.
   memcpy(dst, src, sizeof(Config));    // add all of above from live config
   dst->nhistos = 0; memset(dst->histo_array, 0, MAX_HISTOGRAMS*sizeof(Histogram) );
   // below is wrong - src array lists can contain holes if things were deleted
   //    use same offsets as in src arrays:
   //          offset = src->calib[i] - &src->calib_array[0];
   //   dst->calib[i] = offset + &dst->calib_array[0]
   for(i=0; i<MAX_CALIB;     i++){dst->calib[i]    = &dst->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS;   i++){dst-> globals[i] = &dst->global_array[i]; }
   for(i=0; i<MAX_CONDS;     i++){dst->condlist[i] = &dst->cond_array[i]; }
   for(i=0; i<MAX_GATES;     i++){dst->gatelist[i] = &dst->gate_array[i]; }

   // some of the arrays contain pointers: cond_array has var pointers
   //                                      gate_array has cond pointers
   // these have to be copied the long way
   dst->nconds=0; memset(dst->cond_array, 0, MAX_GATES*sizeof(Cond) );
   dst->ngates=0; memset(dst->gate_array, 0, MAX_GATES*sizeof(Gate) );
   for(i=0; i<src->nconds; i++){ cond = src->condlist[i];
      add_cond(dst, cond->name, cond->var->name, GATEOP(cond->op), cond->value);
   }
   for(i=0; i<src->ngates; i++){ gate = src->gatelist[i];
      add_gate(dst, gate->name);
      for(j=0; j<gate->nconds; j++){ cond = gate->conds[j];
         add_cond_to_gate(dst,gate->name , cond->name);
      }
   }

   // apply_gates below takes care of var->histo_list_x
   for(i=0; i<MAX_SORT_VARS; i++){dst->varlist[i].use_count_x = 0; }
   dst->nusedvar = dst->nuser = 0; // add_histos takes care of these
   // copy config histograms  ODB histos will follow later
   for(i=0; i<src->nhistos; i++){
      histo = src->histo_list[i];
      tmp = ( histo->ybins ) ? histo->yvar->name : NULL;
      if( add_histo(dst, histo->handle, histo->title, histo->path, histo->xbins, histo->xvar->name, 0, histo->xbins, histo->ybins, tmp, 0, histo->ybins) ){ return(-1); }
      // apply gates ...
      for(j=0; j<histo->num_gates; j++){
         apply_gate(dst, histo->handle, histo->gate_names[j]);
      }
   }
   /*
   // update user_histo_list to point to new histos
   for(i=0; i<src->nuser; i++){
      srchist = src->user_histos[i];
      if( (histo = find_histo(dst, srchist->handle) ) == NULL ){
            printf("copy_config: impossible error#1\n"); continue;
      }
      dst->user_histos[i] = histo;
   }
   // update used_vars list to point to dst->sortvars (from src->sortvars)
   for(i=0; i<src->nusedvar; i++){
      srcvar = src->usedvars[i];
      if( (dstvar = find_sortvar(dst, srcvar->name)) == NULL ){
         printf("copy_config: impossible error#2\n"); continue;
      }
      dst->usedvars[i] = dstvar;
   }
   // update sortvar histo list to point to dst histos
   for(i=0; i<src->nsortvar; i++){
      srcvar = &src->varlist[i];
      dstvar = &dst->varlist[i];
      for(j=0; j<srcvar->use_count_x; j++){
         srchist = srcvar->histo_list_x[i];
         if( (histo = find_histo(dst, srchist->handle) ) == NULL ){
            printf("copy_config: impossible error#3\n"); continue;
         }
         dstvar->histo_list_x[i] = histo;
      }
   }
   */
   /* tmp  = (char *)dst;
   tmp2 = (char *)src;
   for(i=0; i<sizeof(Config); i+=8){
      tmp3 = *(long *)(tmp+i);
      ptr = (char *)tmp3;
      if( ptr >= tmp2 && ptr <= tmp2+sizeof(Config) ){
         printf("offset:%d\n", i);
      }
   }
   */
   //dst->odb_daqsize = src->odb_daqsize;
   src->lock = 0; return(0);
}

// needed to update calibrations during sort
// dst is "sort" config, src is new stuff to merge
int merge_configs(Config *src, Config *dst)
{
   Cal_coeff *cal;
   Global *global;
   int i;

   dst->lock = 1;  src->lock = 1;
   for(i=0; i<src->ncal;      i++){ cal = src->calib[i];
      edit_calibration(dst, cal->name, cal->offset, cal->gain, cal->quad, cal->pileupk1, cal->pileupk2, cal->pileupE1, cal->crosstalk0, cal->crosstalk1, cal->crosstalk2, cal->address, cal->datatype, 1);
   }
   for(i=0; i<src->nglobal;   i++){ global = src->globals[i];
      add_global(dst, global->name, global->min, global->max);
   }
   src->lock = 0; dst->lock = 0; return(0);
}

Config *add_config(char *name)
{
   Config *cfg;
   int i, len;

   for(i=0; i<MAX_CONFIGS; i++){
      if( configs[i] == NULL ){ break; }
   }
   if( i == MAX_CONFIGS ){
      fprintf(stderr,"Exceeded Max histogram sets\n");
   }
   if( (cfg = configs[i] = calloc(1, sizeof(Config))) == NULL ){
      fprintf(stderr,"can't alloc new config\n");
   }
   if( (len=strlen(name)+1) >  SYS_PATH_LENGTH ){
      fprintf(stderr,"truncating configname: %s\n", name);
      len = SYS_PATH_LENGTH;
   }
   memcpy(cfg->name, name, len);
   for(i=0; i<MAX_CALIB;   i++){cfg->calib[i]    = &cfg->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS; i++){cfg->globals[i] = &cfg->global_array[i]; }
   for(i=0; i<MAX_CONDS;   i++){cfg->condlist[i] = &cfg->cond_array[i]; }
   for(i=0; i<MAX_GATES;  i++){cfg->gatelist[i] = &cfg->gate_array[i]; }
   init_user_config(cfg); // setup hardcoded variables
   return( cfg );
}

int remove_config(Config *cfg)
{
   int i;

   if( cfg == configs[0] || cfg == configs[0] ){ return(0); }// do not remove set#0 or #1
   for(i=1; i<MAX_CONFIGS; i++){
      if( cfg == configs[i] ){ break; }
   }
   if( i == MAX_CONFIGS ){
      fprintf(stderr,"can't find config to remove\n");
   }
   clear_config(cfg);
   free(cfg);
   configs[i] = NULL;
   return(0);
}

int save_config(Config *cfg, char *filename, int overwrite)
{
   FILE *fp;
   if( cfg->lock ){ return(-1); } // config file currently being read
   if( !overwrite ){
      if( (fp=fopen(filename,"r")) != NULL ){
         fprintf(stderr,"save_config: file %s already exists\n", filename);
         fclose(fp); return(-1);
      }
   }
   if( (fp=fopen(filename,"w")) == NULL ){
      fprintf(stderr,"save_config: cant open %s to write\n", filename);
      return(-1);
   }
   write_config(cfg, fp);
   fclose(fp);
   return(0);
}

/////////////////////////////   Variable   /////////////////////////////////

// THERE IS CURRENTLY NO WAY TO CALCULATE THE VALUE OF A VARIABLE
//   GIVEN ONLY ITS TEXT DESCRIPTION
//      SO THIS FUNCTION WOULD BE POINTLESS AT THE MOMENT

//int add_variable(Config *cfg, char *name, char *title)
//{
//   time_t current_time = time(NULL);
//   int i = cfg->nsortvar;
//   if( i == MAX_SORT_VARS ){
//      fprintf(stderr,"too many variables when adding %s\n", name); return(-1);
//   }
//   memcpy(cfg->varlist[i].name, name, strlen(name)+1);
//   memcpy(cfg->varlist[i].title, title, strlen(title)+1);
//   ++cfg->nsortvar;
//   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
//  return(0);
//}

Sortvar *find_sortvar(Config *cfg, char *name)
{
   int i;
   for(i=0; i<cfg->nsortvar; i++){
      if( strcmp(cfg->varlist[i].name, name) == 0 &&
          strlen(cfg->varlist[i].name) == strlen(name) ){
         return( &cfg->varlist[i] );
      }
   }
   return(NULL);
}

/////////////////////////  Global Condition   //////////////////////////////

// every array member has a pointer to itself
//    [the pointer list entries are swapped etc, but not cleared]
// same with gates, conditions, variables
// but *not* used-variables, user-histos
int add_global(Config *cfg, char *name, int value, int val2)
{
   time_t current_time = time(NULL);
   Global *global;
   int i, len;
   for(i=0; i<cfg->nglobal; i++){
      if( strcmp(cfg->globals[i]->name, name) == 0 &&
          strlen(cfg->globals[i]->name) == strlen(name) ){ break; }
   }
   global = cfg->globals[i];
   if( i == cfg->nglobal ){ // new global - i is first unused ptr
      if( cfg->nglobal >= MAX_GLOBALS ){
         fprintf(stderr,"too many globals for %s\n", name); return(-1);
      }
      if( (len=strlen(name)+1) > STRING_LEN ){
         fprintf(stderr,"truncating globalname: %s\n", name);
         len = STRING_LEN;
      }
      memcpy(global->name, name, len);
       ++cfg->nglobal;
   }
   global->min = value; global->max = val2;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int remove_global(Config *cfg, char *name)
{
   Global *global, *lastglobal = cfg->globals[cfg->nglobal];
   time_t current_time = time(NULL);
   int i;
   for(i=0; i<cfg->nglobal; i++){
      if( strcmp(cfg->globals[i]->name, name) == 0 &&
          strlen(cfg->globals[i]->name) == strlen(name) ){
         global = cfg->globals[i]; break;
      }
   }
   if( i == cfg->nglobal ){
      fprintf(stderr,"can't find global: %s to remove\n", name); return(-1);
   }
//   if( global->use_count != 0 ){
//      fprintf(stderr,"global[%s] still in use[%d]\n", name, global->use_count);
//      return(-1);
//   }
   // if removing final entry - no rearrangement needed
   if( i != cfg->nglobal-1 ){  // otherwise - swap pointers with last
      cfg->globals[i             ] = lastglobal;
      cfg->globals[cfg->nglobal-1] = global;
   }
   --cfg->nglobal;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

/////////////////////////////   GATE   /////////////////////////////////

// after editing gates, may end up with some unreferenced conditions ...
//    until config saved and reloaded
int add_gate(Config *cfg, char *name)
{
   Gate *gate = cfg->gatelist[cfg->ngates];
   time_t current_time = time(NULL);
   int i, len;
   for(i=0; i<cfg->ngates; i++){
      if( strncmp(name, cfg->gatelist[i]->name, strlen(name)) == 0 &&
                 strlen(cfg->gatelist[i]->name) == strlen(name)    ){
//         fprintf(stderr,"gateGate[%s] already exists\n", name); return(-1);
         gate = cfg->gatelist[i];
         gate->nconds = 0;
         cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
         return(0); // this is now acceptable for editing exisiting gates
      }
   }
   if( cfg->ngates >= MAX_GATES ){
      fprintf(stderr,"too many gategates: %d\n", MAX_GATES); return(-1);
   }
   if( (len=strlen(name)+1) > STRING_LEN ){
      fprintf(stderr,"truncating gatename: %s\n", name);
      len = STRING_LEN;
   }
   memcpy(gate->name, name, len);
   gate->nconds = 0; ++cfg->ngates;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int remove_gate(Config *cfg, char *name)
{
   Gate *gate, *lastgate = cfg->gatelist[cfg->ngates-1];
   time_t current_time = time(NULL);
   int i;
   for(i=0; i<cfg->ngates; i++){
      if( strcmp(cfg->gatelist[i]->name, name) == 0 &&
          strlen(cfg->gatelist[i]->name) == strlen(name) ){
         gate = cfg->gatelist[i]; break;
      }
   }
   if( i == cfg->ngates ){
      fprintf(stderr,"can't find gate: %s to remove\n", name); return(-1);
   }
   if( gate->use_count != 0 ){
      fprintf(stderr,"gate[%s] still in use[%d]\n", name, gate->use_count);
      return(-1);
   }
    // removing final entry - no rearrangement needed
   if( i != cfg->ngates-1 ){ // removing middle entry - swap pointers with last
      cfg->gatelist[i        ] = lastgate;
      cfg->gatelist[cfg->ngates-1] = gate;
   }
   --cfg->ngates;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

//  add additional gate of conditions to a histogram
int apply_gate(Config *cfg, char *histoname, char *gatename)
{
   Histogram *histo = cfg->histo_list[0];
   time_t current_time = time(NULL);
   Gate *gate;
   int i;

   for(i=0; i<cfg->ngates; i++){
      if( strcmp(cfg->gatelist[i]->name, gatename) == 0 &&
          strlen(cfg->gatelist[i]->name) == strlen(gatename) ){
         gate = cfg->gatelist[i]; break;
      }
   }
   if( i == cfg->ngates ){
      fprintf(stderr,"Apply gate: can't find gate[%s]\n", gatename);
      return(-1);
   }
   for(i=0; i<cfg->nhistos; i++){ histo = cfg->histo_list[i];
      if( strcmp(histo->handle, histoname) == 0 &&
          strlen(histo->handle) == strlen(histoname) ){ break; }
   }
   if( i == cfg->nhistos ){
      fprintf(stderr,"Apply gate can't find histogram: %s\n", histoname);
      return(-1);
   }
   if( (i = histo->num_gates) >= MAX_GATES ){
      fprintf(stderr,"Apply gate too many gates on histogram\n");
      return(-1);
   }
   histo-> gate_names[i] =  gate->name;
   histo->gate_passed[i] = &gate->passed;
   ++histo->num_gates;
   ++gate->use_count;
   for(i=0; i<gate->nconds; i++){ ++gate->conds[i]->use_count; }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int unapply_gate(Config *cfg, char *histoname, char *gatename)
{
   time_t current_time = time(NULL);
   Histogram *histo, *tmp;
   Gate *gate;
   int i;

   for(i=0; i<cfg->nhistos; i++){ tmp = cfg->histo_list[i];
      if( strcmp(tmp->handle, histoname) == 0 &&
          strlen(tmp->handle) == strlen(histoname) ){ histo = tmp; break; }
   }
   if( i == cfg->nhistos ){
      fprintf(stderr,"UnApply gate: can't find histo[%s]\n", histoname);
      return(-1);
   }
   for(i=0; i<cfg->ngates; i++){
      if( strcmp(cfg->gatelist[i]->name, gatename) == 0 &&
          strlen(cfg->gatelist[i]->name) == strlen(gatename) ){
         gate = cfg->gatelist[i]; break;
      }
   }
   if( i == cfg->ngates ){
      fprintf(stderr,"UnApply gate: can't find gate[%s]\n", gatename);
      return(-1);
   }
   for(i=0; i<histo->num_gates; i++){
      if( strcmp(histo-> gate_names[i], gatename) == 0 &&
          strlen(histo-> gate_names[i]) == strlen(gatename) ){ break; }
   }
   if( i == histo->num_gates ){
      fprintf(stderr,"Unapply gate: gate[%s] not applied to histo[%s]\n",
              gatename, histoname); return(-1);
   }
   if( i != (histo->num_gates)-1 ){ // not last one - rearrange
      histo->gate_names[i] = histo->gate_names[histo->num_gates-1];
      histo->gate_names[i] = histo->gate_names[histo->num_gates-1];
   }
   --histo->num_gates;
   --gate->use_count;
   for(i=0; i<gate->nconds; i++){ --gate->conds[i]->use_count; }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

///////////////////////////   CONDITION  /////////////////////////////////

// viewer submits unnamed gates - choose generic name gateXX
int next_condname(Config *cfg)
{
   char tmp[16];
   int i, j;

   for(j=0; j<MAX_CONDS; j++){ sprintf(tmp,"Cond%d", j);
      for(i=0; i<cfg->nconds; i++){
         if( strcmp(cfg->condlist[i]->name, tmp) == 0 &&
             strlen(cfg->condlist[i]->name) == strlen(tmp) ){ break; }
      }
      if( i == cfg->nconds ){ return(i); }
   }
   fprintf(stderr,"next_condname:all conds un use\n");
   return(-1);
}

int add_cond(Config *cfg, char *name, char *varname, char *op, int value)
//int add_cond(char *name, char *varname, int   op, int value)
{
   Cond *cond = cfg->condlist[cfg->nconds];
   time_t current_time = time(NULL);
   Sortvar *var;
   int i, len;

   if( cfg->nconds >= MAX_CONDS ){
      fprintf(stderr,"too many conds: %d\n", MAX_CONDS); return(-1);
   }
   for(i=0; i<cfg->nconds; i++){
      if( strncmp(name, cfg->condlist[i]->name, strlen(name)) == 0 &&
              strlen(name) == strlen(cfg->condlist[i]->name)       ){
         fprintf(stderr,"cond[%s] already exists\n", name); return(-1);
      }
   }
   for(i=0; i<cfg->nsortvar; i++){
      if( strcmp(cfg->varlist[i].name, varname) == 0 &&
          strlen(cfg->varlist[i].name) == strlen(varname) ){
         var = &cfg->varlist[i]; break;
      }
   }
   if( i == cfg->nsortvar ){
      fprintf(stderr,"add_cond: can't find variable[%s]\n", varname);
      return(-1);
   }
   if( (len=strlen(name)+1) > STRING_LEN ){
      fprintf(stderr,"truncating condname: %s\n", name);
      len = STRING_LEN;
   }
   memcpy(cond->name, name, len);
   if(        strcmp(op,"LT")  == 0 ){ cond->op = GATEOP_LT;
   } else if( strcmp(op,"LTE") == 0 ){ cond->op = GATEOP_LE;
   } else if( strcmp(op,"LE")  == 0 ){ cond->op = GATEOP_LE;
   } else if( strcmp(op,"GT")  == 0 ){ cond->op = GATEOP_GT;
   } else if( strcmp(op,"GE")  == 0 ){ cond->op = GATEOP_GE;
   } else if( strcmp(op,"GTE") == 0 ){ cond->op = GATEOP_GE;
   } else if( strcmp(op,"EQ")  == 0 ){ cond->op = GATEOP_EQ;
   } else if( strcmp(op,"<=")  == 0 ){ cond->op = GATEOP_LE;
   } else if( strcmp(op,"<" )  == 0 ){ cond->op = GATEOP_LT;
   } else if( strcmp(op,">=")  == 0 ){ cond->op = GATEOP_GE;
   } else if( strcmp(op,">" )  == 0 ){ cond->op = GATEOP_GT;
   } else if( strcmp(op,"=" )  == 0 ){ cond->op = GATEOP_EQ;
   } else if( strcmp(op,"RA")  == 0 ){ cond->op = GATEOP_RA;
   } else {
      fprintf(stderr,"unknown op:%s for cond %s\n", op, name); return(-1);
   }
   // check for veto-type condition (if zero passes condition)
   if(        cond->op == GATEOP_LT ){ cond->veto = (value >  0);
   } else if( cond->op == GATEOP_LE ){ cond->veto = (value >= 0);
   } else if( cond->op == GATEOP_GT ){ cond->veto = (value <  0);
   } else if( cond->op == GATEOP_GE ){ cond->veto = (value <= 0);
   } else if( cond->op == GATEOP_EQ ){ cond->veto = (value == 0);
   }
   cond->value = value; cond->var = var; cond->use_count = 0; ++cfg->nconds;    cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

// conds have a use-count - do not allow removal if non-zero
int remove_cond(Config *cfg, char *name)
{
   Cond *cond, *lastcond = cfg->condlist[cfg->nconds-1];
   time_t current_time = time(NULL);
   int i;
   for(i=0; i<cfg->nconds; i++){
      if( strcmp(cfg->condlist[i]->name, name) == 0 &&
          strlen(cfg->condlist[i]->name) == strlen(name) ){
         cond = cfg->condlist[i]; break;
      }
   }
   if( i == cfg->nconds ){
      fprintf(stderr,"can't find cond: %s to remove\n", name); return(-1);
   }
   if( cond->use_count != 0 ){
      fprintf(stderr,"cond[%s] still in use[%d]\n", name, cond->use_count);
      return(-1);
   }
   // if removing final entry - no rearrangement needed
   if( i == cfg->nconds-1 ){ // removing middle entry - swap pointers with last
      cfg->condlist[i        ] = lastcond;
      cfg->condlist[cfg->nconds-1] = cond;
   }
   --cfg->nconds;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int add_cond_to_gate(Config *cfg, char *gatename, char *condname)
{
   time_t current_time = time(NULL);
   Gate *gate;
   Cond *cond;
   int i;

   for(i=0; i<cfg->ngates; i++){
      if( strncmp(gatename,cfg->gatelist[i]->name,strlen(gatename)) == 0 &&
              strlen(cfg->gatelist[i]->name) == strlen(gatename) ){
         gate = cfg->gatelist[i]; break;
      }
   }
   if( i == cfg->ngates ){
      fprintf(stderr,"addCondToGate:can't find gate %s\n",gatename);
      return(-1);
   }
   if( gate->nconds >= MAX_GATE_CONDS ){
      fprintf(stderr,"too many conds in gate\n"); return(-1);
   }
   for(i=0; i<cfg->nconds; i++){
      if( strcmp(cfg->condlist[i]->name, condname) == 0 &&
          strlen(cfg->condlist[i]->name) == strlen(condname) ){
         cond = cfg->condlist[i]; break;
      }
   }
   if( i == cfg->nconds ){
      fprintf(stderr,"addCondToGate:can't find cond %s\n", condname);
      return(-1);
   }
   gate->conds[gate->nconds++] = cond;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

/////////////////////////////   HISTOGRAM   ///////////////////////////////

//int add_gate(Config *cfg, char *name)
////int add_gate(char *name, char *varname, int   op, int value)
//{
//}

// "user" histograms added by viewer are not currently filled during sort
// either loop over all user histos and fill
//     -> for each histo check var non-zero
// or     loop over non-zero variables, and fill dependant histos
//     -> for each non-zero var[at time it is calculated] fill its histos
// second does not waste any time, but requires cfgting up histo lists
//
// don't bother storing pointers to whole sortvar within histograms...
// when removing histo can loop over all vars to remove histo dependance
//
// add any gates one by one after this, using above functions
int add_histo(Config *cfg, char *name, char *title, char *path, int xbins, char *xvarname, int xmin, int xmax, int ybins, char *yvarname, int ymin, int ymax){
   time_t current_time = time(NULL);
   Sortvar *xvar, *yvar;
   Histogram *tmp;
   int i;

   for(i=0; i<cfg->nhistos; i++){ tmp = cfg->histo_list[i];
      if( strcmp(tmp->handle, name) == 0 &&
          strlen(tmp->handle) == strlen(name) ){
         //fprintf(stderr,"add_histo[%s] already exists\n", name); return(-1);
         remove_histo(cfg, name); break;
      }
   }
   for(i=0; i<cfg->nsortvar; i++){
      if( strcmp(cfg->varlist[i].name, xvarname) == 0 &&
          strlen(cfg->varlist[i].name) == strlen(xvarname) ){
         xvar = &cfg->varlist[i]; break;
      }
   }
   if( i == cfg->nsortvar ){
      fprintf(stderr,"add_histo: can't find x variable[%s]\n", xvarname);
      return(-1);
   }
   if( ybins == 0 ){ // 1D histograms
      if( (tmp=H1_BOOK(cfg,name,title,xbins,xmin,xmax)) == NULL){return(-1);}
   } else {          // 2D histograms
      for(i=0; i<cfg->nsortvar; i++){
         if( strcmp(cfg->varlist[i].name, yvarname) == 0 &&
             strlen(cfg->varlist[i].name) == strlen(yvarname) ){
            yvar = &cfg->varlist[i]; break;
         }
      }
      if( i == cfg->nsortvar ){
         fprintf(stderr,"add_histo: can't find y variable[%s]\n", yvarname);
         return(-1);
      }
      if( (tmp=(Histogram *)H2_BOOK(cfg,name,title,xbins,xmin,xmax,ybins,ymin,ymax)) == NULL ){ return(-1); }
       tmp->yvar = yvar;
   }
   memcpy(tmp->path, path, strlen(path)+1 );
   tmp->xvar = xvar; xvar->histo_list_x[xvar->use_count_x++] = tmp;
   tmp->user = 1;  cfg->user_histos[cfg->nuser++] = tmp;
   // add xvar to used-vars if not already there
   for(i=0; i<cfg->nusedvar; i++){
      if( strcmp(cfg->usedvars[i]->name, xvarname) == 0 &&
          strlen(cfg->usedvars[i]->name) == strlen(xvarname) ){ break; }
   }
   if( i == cfg->nusedvar ){ cfg->usedvars[cfg->nusedvar++] = xvar; }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

// histo->next - un-needed now have flat list?
//
//
int remove_histo(Config *cfg, char *name)
{
   time_t current_time = time(NULL);
   Histogram *histo, *tmp;
   Sortvar *var;
   int i, j, k;

   for(i=0; i<cfg->nhistos; i++){ tmp = cfg->histo_list[i];
      if( strcmp(tmp->handle, name) == 0 &&
          strlen(tmp->handle) == strlen(name) ){ histo = tmp; break; }
   }
   if( i == cfg->nhistos ){
      fprintf(stderr,"RemoveHisto: can't find histo[%s]\n", name);
      return(-1);
   }
   if( !histo->user ){
      fprintf(stderr,"histo[%s] not user-added, can't remove\n", name);
      return(-1);
   }
   while( histo->num_gates > 0 ){ // unapply all gates
      unapply_gate(cfg, name, histo->gate_names[0] );
   }
   // remove from xvar->histo_list_x
   var = histo->xvar;
   for(i=0; i<var->use_count_x; i++){
      if( var->histo_list_x[i] == histo ){
         if( i != var->use_count_x - 1 ){ // not last one - rearrange
            var->histo_list_x[i] = var->histo_list_x[var->use_count_x-1];
         }
         --var->use_count_x;
         break;
      }
   }
   // remove from cfg->user_histo list
   for(i=0; i<cfg->nuser; i++){
      if( cfg->user_histos[i] == histo ){
         if( i != cfg->nuser - 1 ){ // not last one - rearrange
            cfg->user_histos[i] = cfg->user_histos[cfg->nuser-1];
         }
         --cfg->nuser;
         break;
      }
   }
   // remove xvar from used_varlist (if not still in use by another histo)
   for(i=0; i<cfg->nuser; i++){
      if( cfg->user_histos[i]->xvar == var ){ break; }
   }
   if( i == cfg->nuser ){ // not still in use - remove
      for(i=0; i<cfg->nusedvar; i++){
         if( strcmp(cfg->usedvars[i]->name, var->name) == 0 &&
             strlen(cfg->usedvars[i]->name) == strlen(var->name) ){
            if( i != cfg->nusedvar - 1 ){ // not last one - rearrange
               cfg->usedvars[i] = cfg->usedvars[cfg->nusedvar - 1];
            }
            --cfg->nusedvar;
            break;
         }
      }
   }
   // finally remove histo from main histogram list
   cfg->histo_list[i] = cfg->histo_list[cfg->nhistos-1];//nop if i last
   cfg->histo_list[cfg->nhistos-1] = NULL;
   --cfg->nhistos;
   free(histo->data);
   cfg->folders_valid = 0;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);

   //for(i=0; i<sortvarnum; i++){ // check all vars
   //   for(j=0; j<varlist[i].use_count; j++){ // check all dep. histos
   //      if( varlist[i].histo_list[j] == histo ){ // remove it
   //         for(k=j+1; k<varlist[i].use_count; k++){ // memmove
   //            varlist[i].histo_list[k-1] = varlist[i].histo_list[k];
   //         }
   //         --varlist[i].use_count;
   //      }
   //   }
   //}
   return(0);
}

Histogram *find_histo(Config *cfg, char *name)
{
   int i;
   for(i=0; i<cfg->nhistos; i++){
      if( strcmp(cfg->histo_list[i]->handle, name) == 0 &&
          strlen(cfg->histo_list[i]->handle) == strlen(name) ){
         return( cfg->histo_list[i] );
      }
   }
   return(NULL);
}

/////////////////////////////   CALIBRATION   /////////////////////////////

int set_calibration(Config *cfg, int num, char url_args[][STRING_LEN], int fd)
{
   float offset, gain, quad;
   float puk1[7], puk2[7], puE1[7];
   float ct0[16], ct1[16], ct2[16];
   int i, address=-1, datatype=-1;
   char tmp[128];

   // Initialize values to -1
   for(i=0; i<7; i++){
     puk1[i] = puk2[i] = puE1[i] = -1;
   }
   for(i=0; i<16; i++){
     ct0[i]  = ct1[i]  = ct2[i]  = -1;
   }

   for(i=2; i<num; i+=8){
      if( strncmp(url_args[i], "channelName", 10) != 0 ){
         sprintf(tmp,"set_calibration: expected \"channelName\" at %s\n",url_args[i]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"channelName\" at %s\n", url_args[i]);
         return(-1);
      }
      if( strncmp(url_args[i+2], "quad", 4) != 0 ){
         sprintf(tmp,"set_calibration: expected \"quad\" at %s\n",url_args[i+2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"quad\" at %s\n", url_args[i+2]);
         return(-1);
      }
      if( sscanf(url_args[i+3], "%f", &quad) < 1 ){
         sprintf(tmp,"set_calibration: can't read quad value, %s\n",url_args[i+3]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read quad: %s\n", url_args[i+3]);
         return(-1);
      }
      if( strncmp(url_args[i+4], "gain", 4) != 0 ){
         sprintf(tmp,"set_calibration: expected \"gain\" at %s\n",url_args[i+4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"gain\" at %s\n", url_args[i+4]);
         return(-1);
      }
      if( sscanf(url_args[i+5], "%f", &gain) < 1 ){
         sprintf(tmp,"set_calibration: can't read gain value, %s\n",url_args[i+5]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read gain: %s\n", url_args[i+5]);
         return(-1);
      }
      if( strncmp(url_args[i+6], "offset", 6) != 0 ){
         sprintf(tmp,"set_calibration: expected \"offset\" at %s\n",url_args[i+6]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"offset\" at %s\n", url_args[i+6]);
         return(-1);
      }
      if( sscanf(url_args[i+7], "%f", &offset) < 1 ){
         sprintf(tmp,"set_calibration: can't read offset value, %s\n",url_args[i+7]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read offset: %s\n", url_args[i+7]);
         return(-1);
      }
//      if( strncmp(url_args[i+8], "address", 4) != 0 ){
//         fprintf(stderr,"expected \"address\" at %s\n", url_args[i+8]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+9], "%d", &address) < 1 ){
//         fprintf(stderr,"can't read address: %s\n", url_args[i+9]);
//         return(-1);
//      }
//      if( strncmp(url_args[i+10], "datatype", 6) != 0 ){
//         fprintf(stderr,"expected \"datatype\" at %s\n", url_args[i+10]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+11], "%d", &datatype) < 1 ){
//         fprintf(stderr,"can't read datatype: %s\n", url_args[i+11]);
//         return(-1);
//      }
      // Send the response header
      send_header(fd, APP_JSON);

      edit_calibration(cfg, url_args[i+1], offset, gain, quad, puk1, puk2, puE1, ct0, ct1, ct2,
                       address, datatype, 1);
   }
   return(0);
}

int set_pileup_correction(Config *cfg, int num, char url_args[][STRING_LEN], int fd)
{
   float offset=-1, gain=-1, quad=-1;
   float puk1[7], puk2[7], puE1[7];
   float ct0[16], ct1[16], ct2[16];
   int i, address=-1, datatype=-1;
   char tmp[128];

// Initialize values to defaults
for(i=0; i<7; i++){
  puk1[i] = puk2[i] = puE1[i] = 0;
}
puk1[0] = puk2[0] = 1; // set default factor as 1 not zero
for(i=0; i<16; i++){
  ct0[i]  = ct1[i]  = ct2[i]  = -1;
}

   for(i=2; i<num; i+=8){
      if( strncmp(url_args[i], "channelName", 10) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"channelName\" at %s\n",url_args[i]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"channelName\" at %s\n", url_args[i]);
         return(-1);
      }
      if( strncmp(url_args[i+2], "pileupk1", 8) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"pileupk1\" at %s\n",url_args[i+2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"pileupk1\" at %s\n", url_args[i+2]);
         return(-1);
      }
      if( sscanf(url_args[i+3], "%f,%f,%f,%f,%f,%f,%f", puk1,puk1+1,puk1+2,puk1+3,puk1+4,puk1+5,puk1+6) != 7 ){
         sprintf(tmp,"set_pileup_correction: can't read pileup k1 value (expected 7 values), %s\n",url_args[i+3]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read pileup k1, expected 7 values: %s\n", url_args[i+3]);
         return(-1);
      }
      if( strncmp(url_args[i+4], "pileupk2", 8) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"pileupk2\" at %s\n",url_args[i+4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"pileupk2\" at %s\n", url_args[i+4]);
         return(-1);
      }
      if( sscanf(url_args[i+5], "%f,%f,%f,%f,%f,%f,%f", puk2,puk2+1,puk2+2,puk2+3,puk2+4,puk2+5,puk2+6) != 7 ){
         sprintf(tmp,"set_pileup_correction: can't read pileup k2 value (expected 7 values), %s\n",url_args[i+5]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read pileup k2, expected 7 values: %s\n", url_args[i+5]);
         return(-1);
      }
      if( strncmp(url_args[i+6], "pileupE1", 8) != 0 ){
         sprintf(tmp,"set_pileup_correction: expected \"pileupE1\" at %s\n",url_args[i+6]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"pileupE1\" at %s\n", url_args[i+6]);
         return(-1);
      }
      if( sscanf(url_args[i+7], "%f,%f,%f,%f,%f,%f,%f", puE1,puE1+1,puE1+2,puE1+3,puE1+4,puE1+5,puE1+6) != 7 ){
         sprintf(tmp,"set_pileup_correction: can't read pileup E1 value (expected 7 values), %s\n",url_args[i+7]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read pileup E1, expected 7 values: %s\n", url_args[i+7]);
         return(-1);
      }
//      if( strncmp(url_args[i+14], "address", 4) != 0 ){
//         fprintf(stderr,"expected \"address\" at %s\n", url_args[i+14]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+15], "%d", &address) < 1 ){
//         fprintf(stderr,"can't read address: %s\n", url_args[i+15]);
//         return(-1);
//      }
//      if( strncmp(url_args[i+16], "datatype", 6) != 0 ){
//         fprintf(stderr,"expected \"datatype\" at %s\n", url_args[i+16]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+17], "%d", &datatype) < 1 ){
//         fprintf(stderr,"can't read datatype: %s\n", url_args[i+17]);
//         return(-1);
//      }
      // Send the response header
      send_header(fd, APP_JSON);

    edit_calibration(cfg, url_args[i+1], offset, gain, quad, puk1, puk2, puE1, ct0, ct1, ct2, address, datatype, 1);
   }
   return(0);
}

int set_crosstalk_correction(Config *cfg, int num, char url_args[][STRING_LEN], int fd)
{
   float offset=-1, gain=-1, quad=-1;
   float puk1[7], puk2[7], puE1[7];
   float ct0[16], ct1[16], ct2[16];
   int i, address=-1, datatype=-1;
   char tmp[128];

// Initialize values to defaults
for(i=0; i<7; i++){
  puk1[i] = puk2[i] = puE1[i] = -1;
}
for(i=0; i<16; i++){
  ct0[i]  = ct1[i]  = ct2[i]  = 0;
}

   for(i=2; i<num; i+=8){
      if( strncmp(url_args[i], "channelName", 10) != 0 ){
         sprintf(tmp,"set_crosstalk_correction: expected \"channelName\" at %s\n",url_args[i]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"channelName\" at %s\n", url_args[i]);
         return(-1);
      }
      if( strncmp(url_args[i+2], "crosstalk0", 8) != 0 ){
         sprintf(tmp,"set_crosstalk_correction: expected \"crosstalk0\" at %s\n",url_args[i+2]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"crosstalk0\" at %s\n", url_args[i+2]);
         return(-1);
      }
      if( sscanf(url_args[i+3], "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", ct0,ct0+1,ct0+2,ct0+3,ct0+4,ct0+5,ct0+6,ct0+7,ct0+8,ct0+9,ct0+10,ct0+11,ct0+12,ct0+13,ct0+14,ct0+15) != 16 ){
         sprintf(tmp,"set_crosstalk_correction: can't read crosstalk k1 value (expected 7 values), %s\n",url_args[i+3]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read crosstalk k1, expected 16 values: %s\n", url_args[i+3]);
         return(-1);
      }
      if( strncmp(url_args[i+4], "crosstalk1", 8) != 0 ){
         sprintf(tmp,"set_crosstalk_correction: expected \"crosstalk1\" at %s\n",url_args[i+4]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"crosstalk1\" at %s\n", url_args[i+4]);
         return(-1);
      }
      if( sscanf(url_args[i+5], "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", ct1,ct1+1,ct1+2,ct1+3,ct1+4,ct1+5,ct1+6,ct1+7,ct1+8,ct1+9,ct1+10,ct1+11,ct1+12,ct1+13,ct1+14,ct1+15) != 16 ){
         sprintf(tmp,"set_crosstalk_correction: can't read pileup k2 value (expected 7 values), %s\n",url_args[i+5]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read pileup k2, expected 16 values: %s\n", url_args[i+5]);
         return(-1);
      }
      if( strncmp(url_args[i+6], "crosstalk2", 8) != 0 ){
         sprintf(tmp,"set_crosstalk_correction: expected \"crosstalk2\" at %s\n",url_args[i+6]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"expected \"crosstalk2\" at %s\n", url_args[i+6]);
         return(-1);
      }
      if( sscanf(url_args[i+7], "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", ct2,ct2+1,ct2+2,ct2+3,ct2+4,ct2+5,ct2+6,ct2+7,ct2+8,ct2+9,ct2+10,ct2+11,ct2+12,ct2+13,ct2+14,ct2+15) != 16 ){
         sprintf(tmp,"set_crosstalk_correction: can't read pileup E1 value (expected 7 values), %s\n",url_args[i+7]);
         send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
         fprintf(stderr,"can't read pileup E1, expected 16 values: %s\n", url_args[i+7]);
         return(-1);
      }
//      if( strncmp(url_args[i+14], "address", 4) != 0 ){
//         fprintf(stderr,"expected \"address\" at %s\n", url_args[i+14]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+15], "%d", &address) < 1 ){
//         fprintf(stderr,"can't read address: %s\n", url_args[i+15]);
//         return(-1);
//      }
//      if( strncmp(url_args[i+16], "datatype", 6) != 0 ){
//         fprintf(stderr,"expected \"datatype\" at %s\n", url_args[i+16]);
//         return(-1);
//      }
//      if( sscanf(url_args[i+17], "%d", &datatype) < 1 ){
//         fprintf(stderr,"can't read datatype: %s\n", url_args[i+17]);
//         return(-1);
//      }
      // Send the response header
      send_header(fd, APP_JSON);

    edit_calibration(cfg, url_args[i+1], offset, gain, quad, puk1, puk2, puE1, ct0, ct1, ct2, address, datatype, 1);
   }
   return(0);
}

int edit_calibration(Config *cfg, char *name, float offset, float gain, float quad, float puk1[7], float puk2[7], float puE1[7], float ct0[16], float ct1[16], float ct2[16], int address, int type, int overwrite)
{
   time_t current_time = time(NULL);
   int i,j, len, arg;
   Cal_coeff *cal;
/*
if(strcmp(name,"GRG01WN00A")==0){
  fprintf(stdout,"edit_calibration for %s:\n",name);// debugger breakpoint
}
*/
   for(i=0; i<cfg->ncal; i++){ cal = cfg->calib[i];
      if( strncmp(name, cfg->calib[i]->name, strlen(name)) == 0 &&
          strlen(name) == strlen(cfg->calib[i]->name)       ){ break; }
   }
   if( i < cfg->ncal ){ // calib already exists
     if( overwrite ){
       if( offset != -1 ){ cal->offset = offset; }
       if( gain   != -1 ){ cal->gain = gain; }
       if( quad   != -1 ){ cal->quad = quad; }
       if(strncmp(name,"GRG",3)==0){ // Only process pileup and crosstalk for HPGe
         if( puk1[0] != -1 ){ for(j=0; j<7; j++){cal->pileupk1[j] = puk1[j];} }
         if( puk2[0] != -1 ){ for(j=0; j<7; j++){cal->pileupk2[j] = puk2[j];} }
         if( puE1[0] != -1 ){ for(j=0; j<7; j++){cal->pileupE1[j] = puE1[j];} }
         if( ct0[0] != -1 ){ for(j=0; j<16; j++){cal->crosstalk0[j] = ct0[j];} }
         if( ct1[0] != -1 ){ for(j=0; j<16; j++){cal->crosstalk1[j] = ct1[j];} }
         if( ct2[0] != -1 ){ for(j=0; j<16; j++){cal->crosstalk2[j] = ct2[j];} }
       }
     }
      if( address != -1 ){
         cal->address = address;  cal->datatype = type;
      }
   } else { // not already there ... add new one
      if( cfg->ncal >= MAX_CALIB ){
         fprintf(stderr,"too many calibs: %d\n", MAX_CALIB); return(-1);
      }
      cal = cfg->calib[i]; ++cfg->ncal;
      if( (len=strlen(name)+1) > 64 ){
         fprintf(stderr,"truncating calibname: %s\n", name);
         len = 64;
      }
      memcpy(cal->name, name, len);
      if( offset != -1 ){ cal->offset = offset; }else{ cal->offset = 0; }
      if( gain   != -1 ){ cal->gain = gain; }else{ cal->gain = 1; }
      if( quad   != -1 ){ cal->quad = quad; }else{ cal->quad = 0; }
      if(strncmp(name,"GRG",3)==0){ // Only process pileup and crosstalk for HPGe
        if( puk1[0] != -1 ){
          for(j=0; j<7; j++){ cal->pileupk1[j] = puk1[j]; }
        }else{
          cal->pileupk1[0]=1; cal->pileupk1[1]=0; cal->pileupk1[2]=0; cal->pileupk1[3]=0; cal->pileupk1[4]=0; cal->pileupk1[5]=0; cal->pileupk1[6]=0;
        }
        if( puk2[0] != -1 ){
          for(j=0; j<7; j++){cal->pileupk2[j] = puk2[j]; }
        }else{
          cal->pileupk2[0]=1; cal->pileupk2[1]=0; cal->pileupk2[2]=0; cal->pileupk2[3]=0; cal->pileupk2[4]=0; cal->pileupk2[5]=0; cal->pileupk2[6]=0;
        }
        if( puE1[0] != -1 ){
          for(j=0; j<7; j++){cal->pileupE1[j] = puE1[j]; }
        }else{
          cal->pileupE1[0]=0; cal->pileupE1[1]=0; cal->pileupE1[2]=0; cal->pileupE1[3]=0; cal->pileupE1[4]=0; cal->pileupE1[5]=0; cal->pileupE1[6]=0;
        }
        if( ct0[0] != -1 ){
          for(j=0; j<16; j++){ cal->crosstalk0[j] = ct0[j]; }
        }else{
          cal->crosstalk0[0]=0; cal->crosstalk0[1]=0; cal->crosstalk0[2]=0; cal->crosstalk0[3]=0; cal->crosstalk0[4]=0; cal->crosstalk0[5]=0; cal->crosstalk0[6]=0; cal->crosstalk0[7]=0;
          cal->crosstalk0[8]=0; cal->crosstalk0[9]=0; cal->crosstalk0[10]=0; cal->crosstalk0[11]=0; cal->crosstalk0[12]=0; cal->crosstalk0[13]=0; cal->crosstalk0[14]=0; cal->crosstalk0[15]=0;
        }
        if( ct1[0] != -1 ){
          for(j=0; j<16; j++){ cal->crosstalk1[j] = ct1[j]; }
        }else{
          cal->crosstalk1[0]=0; cal->crosstalk1[1]=0; cal->crosstalk1[2]=0; cal->crosstalk1[3]=0; cal->crosstalk1[4]=0; cal->crosstalk1[5]=0; cal->crosstalk1[6]=0; cal->crosstalk1[7]=0;
          cal->crosstalk1[8]=0; cal->crosstalk1[9]=0; cal->crosstalk1[10]=0; cal->crosstalk1[11]=0; cal->crosstalk1[12]=0; cal->crosstalk1[13]=0; cal->crosstalk1[14]=0; cal->crosstalk1[15]=0;
        }
        if( ct2[0] != -1 ){
          for(j=0; j<16; j++){ cal->crosstalk2[j] = ct2[j]; }
        }else{
          cal->crosstalk2[0]=0; cal->crosstalk2[1]=0; cal->crosstalk2[2]=0; cal->crosstalk2[3]=0; cal->crosstalk2[4]=0; cal->crosstalk2[5]=0; cal->crosstalk2[6]=0; cal->crosstalk2[7]=0;
          cal->crosstalk2[8]=0; cal->crosstalk2[9]=0; cal->crosstalk2[10]=0; cal->crosstalk2[11]=0; cal->crosstalk2[12]=0; cal->crosstalk2[13]=0; cal->crosstalk2[14]=0; cal->crosstalk2[15]=0;
        }
      }
      cal->address = address;  cal->datatype = type;
    }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);

/*
   if(strcmp(name,"GRG01WN00A")==0){
   fprintf(stdout,"ct0: [%f,%f,%f,%f,%f,%f,%f]\n",cal->crosstalk0[0],cal->crosstalk0[1],cal->crosstalk0[2],cal->crosstalk0[3],cal->crosstalk0[4],cal->crosstalk0[5],cal->crosstalk0[6]);
   fprintf(stdout,"ct1: [%f,%f,%f,%f,%f,%f,%f]\n",cal->crosstalk1[0],cal->crosstalk1[1],cal->crosstalk1[2],cal->crosstalk1[3],cal->crosstalk1[4],cal->crosstalk1[5],cal->crosstalk1[6]);
   fprintf(stdout,"ct2: [%f,%f,%f,%f,%f,%f,%f]\n",cal->crosstalk2[0],cal->crosstalk2[1],cal->crosstalk2[2],cal->crosstalk2[3],cal->crosstalk2[4],cal->crosstalk2[5],cal->crosstalk2[6]);
}
*/

   return(0);
}

int set_directory(Config *cfg, char *name, char *path)
{
   time_t current_time = time(NULL);
   int len;


   if( (len=strlen(path)) >= SYS_PATH_LENGTH ){
      fprintf(stderr,"set_directory: path too long[%s]\n", path);
      return(-1);
   }
   if(        strncmp(name, "Data",   4) == 0 ){
      memcpy(cfg->data_dir, path, len+1);
   } else if( strncmp(name, "Histo",  5) == 0 ){
      memcpy(cfg->histo_dir, path, len+1);
   } else if( strncmp(name, "Config", 6) == 0 ){
      memcpy(cfg->config_dir, path, len+1);
   } else {
      fprintf(stderr,"set_directory: Unknown directory:%s\n", name);
      return(-1);
   }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

int set_midas_param(Config *cfg, char *name, char *value)
{
   time_t current_time = time(NULL);
   int len;
   char clean_string[128], *tmp;

   if( (len=strlen(value)) >= SYS_PATH_LENGTH ){
      fprintf(stderr,"set_midas_param: value too long[%s]\n", value);
      return(-1);
   }
   if(        strncmp(name, "Title",   5) == 0 ){
      sprintf(clean_string, "%s", value);
      if( (tmp = strstr(clean_string, "\t"))>0 ){ // This illegal character cannot be handled in the browser
        while( (tmp = strstr(clean_string, "\t"))>0 ){
          strncpy(tmp, " ", 1); // keep the length the same, and make use of the terminating character already in clean_string
        }
        memcpy(cfg->midas_title, clean_string, len+1);
      }else{
        memcpy(cfg->midas_title, value, len+1);
      }
   } else if( strncmp(name, "StartTime",  9) == 0 ){
      if( sscanf(value, "%d", &cfg->midas_start_time) < 1 ){
         fprintf(stderr,"set_midas_param: can't read starttime: %s\n", value);
      }
   } else if( strncmp(name, "Duration", 8) == 0 ){
      if( sscanf(value, "%d", &cfg->midas_runtime) < 1 ){
         fprintf(stderr,"set_midas_param: can't read runtime: %s\n", value);
      }
   } else {
      fprintf(stderr,"set_midas_param: Unknown param:%s\n", name);
      return(-1);
   }
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

///////////////////////////// MISC SCRIPTS ETC ////////////////////////////
static Sortfile filelist[FILE_QLEN];
extern int sum_th1I(Config *dst_cfg, Config *src_cfg, TH1I *src);
int queue_sum_histos(Config *cfg, int num, char url_args[][STRING_LEN], int fd)
{
   int i, j, len, next, narg;
   Sort_status *arg;
   Sortfile *sort;
   FILE *fp;
   char tmp[128];

   if( strncmp(url_args[2], "outputfilename", 14) != 0 ){
      sprintf(tmp,"Sum histos: expected \"outputfilename\" at %s\n",url_args[2]);
      send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
      fprintf(stderr,"expected \"outputfilename\" at %s\n", url_args[2]);
      return(-1);
   }
   if( (fp=fopen(url_args[3],"r")) != NULL ){
      sprintf(tmp,"Sum histos: %s already exists NOT OVERWRITING. Histos will not be summed.\n",url_args[3]);
      send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
      fprintf(stderr,"%s already exists NOT OVERWRITING\n", url_args[3]);
      return(-1);
   }
   if( (fp=fopen(url_args[3],"w")) == NULL ){
      sprintf(tmp,"Sum histos: can't open %s to write\n",url_args[3]);
      send_http_error_response(fd, STATUS_CODE_400,(char*)tmp);
      fprintf(stderr,"can't open %s to write\n", url_args[3]);
      return(-1);
   }
   fclose(fp);

      // Send the response header
      send_header(fd, APP_JSON);

   arg = get_sort_status();
   sort = &filelist[arg->final_filenum]; next = arg->final_filenum + 1;
   if( next >= FILE_QLEN ){ next = 0; }
   if( next == arg->current_filenum ){
      fprintf(stderr,"FILE QUEUE FULL"); return(-1);
   }
   sort->narg = num/2 - 1;
   if( (sort->arg = calloc(sort->narg, sizeof(char *))) == NULL ){
      fprintf(stderr,"sum_histos: can't alloc memory for histo-list\n");
      return(-1);
   }
   if( (sort->arg[0] = malloc( strlen(url_args[3])+1 )) == NULL){
      fprintf(stderr,"sum_histos: can't alloc memory for histo-list\n");
      free_sortfile(sort); return(-1);
   }
   memcpy(sort->arg[0], url_args[3], strlen(url_args[3])+1);
   for(i=4; i<num; i+=2){
      if( strncmp(url_args[i], "filename", 8) != 0 ){
         fprintf(stderr,"expected \"filename\" at %s\n", url_args[i]);
         free_sortfile(sort); return(-1);
      }
      if( (fp=fopen(url_args[i+1],"r")) == NULL ){
         fprintf(stderr,"cant open histofile:%s to read\n", url_args[i+1]);
      }
      j = (i-2)/2; len = strlen(url_args[i+1])+1;
      if( (sort->arg[j] = malloc(len)) == NULL){
         fprintf(stderr,"sum_histos: can't alloc memory for histo-list\n");
         free_sortfile(sort); return(-1);
      }
      memcpy(sort->arg[j], url_args[i+1], len);
   }
   if( ++arg->final_filenum == FILE_QLEN ){ arg->final_filenum = 0; } //wrap
   return(0);
}

int sum_histos(Config *cfg, Sortfile *sort)  // cfg -> configs[0]
{
   Config *sum, *tmp, *tmp_conf;
   int i, j, arg;
   FILE *fp;
   if( sort->narg < 2 || sort->arg == NULL ){
      fprintf(stderr,"sum_histos: bad args\n"); return(-1);
   }
   if( (fp=fopen(sort->arg[0],"w")) == NULL ){
      fprintf(stderr,"can't open %s to write\n", sort->arg[0]);
      return(-1);
   }
   if( (sum=add_config(sort->arg[0])) == NULL ){ return(-1); }
   for(i=1; i<sort->narg; i++){ sort->carg = i;
      if( i == 1 ){
         tmp_conf=read_histofile(sort->arg[1], 1); copy_config(tmp_conf, sum);
         remove_config(tmp_conf);
      }
      if( (tmp = read_histofile(sort->arg[i],0)) == NULL ){ continue; }
      fprintf(stderr,"Adding histograms from %s ...\n", sort->arg[i]);
      for(j=0; j<tmp->nhistos; j++){
         if( tmp->histo_list[j]->data == NULL ){
            read_histo_data(tmp->histo_list[j], tmp->histo_fp );
            if( tmp->histo_list[j]->data == NULL ){
               fprintf(stderr, "sum_histos:cant alloc data for %s:%s\n",
                       sort->arg[i], tmp->histo_list[j]->handle );
               return(-1);
            }
            sum_th1I(sum, tmp, tmp->histo_list[j] );
            free(tmp->histo_list[j]->data); tmp->histo_list[j]->data = NULL;
         } else {
            //sum_th1I(sum,(TH1I *)sum->histo_list[j],(TH1I *)tmp->histo_list[j]);
            sum_th1I(sum, tmp, tmp->histo_list[j] );
         }
      }
      if( i != 1 ){ // update start time and duration for subsequent files
         tmp_conf=read_histofile(sort->arg[i], 1);
         if( tmp_conf->midas_start_time < sum->midas_start_time ){
            sum->midas_start_time = tmp_conf->midas_start_time;
         }
         sum->midas_runtime += tmp_conf->midas_runtime;
         remove_config(tmp_conf);
      }
      remove_config(tmp);
   }
   write_histofile(sum, fp);
   fprintf(stderr,"sum_histos: done\n");
   fclose(fp);
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////          Sorting Control          ///////////////////////
///////////////////////////////////////////////////////////////////////////

// have to run midas module as separate thread
//   it has to run rpc server and uses callbacks when events arrive

/* module directory usually needs to be in LD_LIBRARY_PATH */
#include <dlfcn.h>
int (*midas_module_main)(Sort_status *);
void *midas_module_handle;
char midas_host[64], midas_expt[64];
int load_midas_module(char *host, char *expt)
{
   char *symbol = "midas_module_main";
   char *file = "./midas_module.so";
   void *symbol_handle;

   // FnPtr: FnType *(*name)(arglist); // eg Char *(*func)(int);
   // cast:  (FnType *(*)(arglist))    //    ptr = (char *(*)(int))sym;

   if( host != NULL ){ strncpy(midas_host, host, 64); } else { midas_host[0] = 0; }
   if( expt != NULL ){ strncpy(midas_expt, expt, 64); } else { midas_expt[0] = 0; }
   midas_host[63] = midas_expt[63] = 0;

   if( (midas_module_handle=dlopen(file, RTLD_LAZY)) == NULL){
      printf("%s\n", dlerror() ); return(-1);
   }
   if( (symbol_handle=dlsym (midas_module_handle, symbol)) == NULL ){
      printf("%s\n", dlerror() ); dlclose( midas_module_handle ); return(-1);
   }
   midas_module_main = (int (*)(Sort_status *))symbol_handle;
   return(0);
}
void unload_midas_module()
{
   // must have deleted any pointers to module functions
   //   (will crash if theyre called after closing)
   dlclose(midas_module_handle); midas_module_handle = NULL;
}

static struct stat statbuf;
int send_sort_status(int fd)
{
   static time_t conftime;
   int i, mtime, subrun;
   Sort_status *arg;
   char tmpstr[256];
   Sortfile *tmp;
   long done;

   // Send the response header
   send_header(fd, APP_JSON);

   arg = get_sort_status();
   if( conftime == 0 ){ conftime = time(NULL); }
   if( arg->online_mode ){
      sprintf(tmpstr,"ONLINE %s %d", arg->run_in_progress ? "Running" : "Stopped", configs[0]->mtime );
      put_line(fd, tmpstr, strlen(tmpstr) );
      return(0);
   }
   if( arg->final_filenum == arg->current_filenum ){
      sprintf(tmpstr,"IDLE %d", configs[0]->mtime);
      put_line(fd, tmpstr, strlen(tmpstr) );
      return(0);
   }
   if( arg->sum_mode == 1 ){
      i = arg->current_filenum; tmp = &filelist[i];
      sprintf(tmpstr,"SUMMING HISTOS %d %d %d", configs[0]->mtime,
                                        tmp->carg-1, tmp->narg-1);
      put_line(fd, tmpstr, strlen(tmpstr) );
      return(0);
   }
   i = arg->current_filenum;
   while(i !=  arg->final_filenum){
      tmp = &filelist[i];
      done  = (i == arg->current_filenum) ? arg->midas_bytes  : 0;
      mtime = (i == arg->current_filenum) ? configs[0]->mtime : 0;
      sprintf(tmpstr,"%s %s %d %d %ld %ld %d,\n", tmp->data_dir, tmp->data_name, tmp->run, tmp->subrun, tmp->data_size, done, mtime );
      put_line(fd, tmpstr, strlen(tmpstr) );
      //printf("STATUS:%s", tmpstr);
      if( ++i >= FILE_QLEN ){ i = 0; }
   }
   return(0);
}

#define READ_LIMIT (5*1024*1024) // don't read more than this looking for odb
int read_datafile_info(Sortfile *sort, char *path)
{
   int i, expt=0, ppg=0, flt=0, done=0;
   char tmp[256], *ptr;
   FILE *fp;
   if( (fp=fopen(path,"r")) == NULL ){
      fprintf(stderr,"read_datafile_info: can't read file: %s\n", path);
      return(-1);
   }
   while( (ptr=fgets(tmp, 256, fp)) != NULL ){
      if( ftell(fp) > READ_LIMIT ){ break; }
      while( isspace(*ptr) ){ ++ptr; }
      if( strncmp(ptr,"</dir>",6) == 0 ){ expt = 0; }
      if( strncmp(ptr,"<dir name=\"Run parameters\">",27) == 0 ){ expt = 2; }
      // The Title/comment are in Run parameters
      //    Experiment/Status-items and Edit-on-start contain LINKS
      //if( strncmp(ptr,"<dir name=\"Experiment\">",23)     == 0 ){ expt = 1; }
      //if( strncmp(ptr,"<dir name=\"Edit on start\">", 26) == 0 ){ expt = 3; }
      if( expt && strncmp(ptr,"<key name=\"Run Title\"",21) == 0 ){
         ptr+=21; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[0][i++] = *ptr++; }
         sort->file_info[0][i]=0; done |= 1;
      }
      if( expt && strncmp(ptr,"<key name=\"Comment\"",19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[1][i++] = *ptr++; }
         sort->file_info[1][i]=0; done |= 2;
      }
      if( strncmp(ptr,"<dir name=\"PPG\">",16)            == 0 ){ ppg = 1; }
      if( ppg && strncmp(ptr,"<key name=\"Current\"", 19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[2][i++] = *ptr++; }
         sort->file_info[2][i]=0; ppg = 0; done |= 4;
      }
      if( strncmp(ptr,"<dir name=\"Filter\">",19)         == 0 ){ flt = 1; }
      if( flt && strncmp(ptr,"<key name=\"Current\"", 19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[3][i++] = *ptr++; }
         sort->file_info[3][i]=0; flt = 0; done |= 8;
      }
      if( done == 15 ){ break; } // 1st 4 bits of done all set
   }
   fclose(fp);
   return(0);
}

int send_file_details(char *path, int fd)
{
   char tmp2[256];
   Sortfile *tmp;
   if( (tmp = calloc(sizeof(Sortfile), 1)) == NULL){
      send_header(fd, APP_JSON);
      fprintf(stderr,"send_file_details: failed alloc\n");
      put_line(fd, " \n \n \n \n", 8 ); return(-1);
   }
   send_header(fd, APP_JSON);
   read_datafile_info(tmp, path);
   sprintf(tmp2,"%s\n%s\n%s\n%s\n", tmp->file_info[0], tmp->file_info[1],
                                   tmp->file_info[2], tmp->file_info[3] );
   put_line(fd, tmp2, strlen(tmp2) );
   free(tmp);
   return(0);
}

// split path into name,dir  then get histo and config names
// also do stat to get size
int run_number(Sort_status *arg, Sortfile *sort, char *name);
char *subrun_filename(Sortfile *sort, int subrun);
int add_sortfile(char *path, char *histodir, char *confdir, char *calsrc)
{
   int i, next, plen, dlen, ext_len, hlen, clen;
   char ptr, tmp[256], *fname;
   Sort_status *arg, tmp_stat;
   Config *cfg = configs[0];
   Sortfile *sort;

   arg = get_sort_status();
   if( strcmp(path,"ONLINE" ) == 0 ){
      if( load_midas_module(confdir, calsrc) ){ return(-1); }
      arg->online_mode = 1;
      if( ++arg->final_filenum == FILE_QLEN ){ arg->final_filenum = 0; }
      return(0);
   }
   memset(&tmp_stat, 0, sizeof(Sort_status) );
   sort = &filelist[arg->final_filenum]; next = arg->final_filenum + 1;
   if( next >= FILE_QLEN ){ next = 0; }
   if( next == arg->current_filenum ){
      fprintf(stderr,"FILE QUEUE FULL"); return(-1);
   }
   plen=strlen(path);
   ext_len = ( strncmp(path+plen-4, ".mid", 4) == 0 ) ? 4 : 0;
   for(i=plen; i>=0; i--){ if( path[i] == '/' ){ ++i; break; } }
   if( (dlen = i) == -1 ){ dlen = 0; } // no directory separator in path
   if( (sort->data_dir = malloc(dlen + 2)) == NULL ){
      fprintf(stderr,"can't alloc string for data_dir");
      free_sortfile(&filelist[arg->current_filenum]); return(-1);
   }
   if( dlen == 0 ){ sprintf(sort->data_dir, ".");
   } else {
      memcpy((char *)sort->data_dir, path, dlen-1);
      *(sort->data_dir+dlen-1) = 0; // overwrite trailing '/'
      set_directory(cfg, "Data", sort->data_dir);
   }
   if( (sort->data_name = malloc(plen-dlen+1)) == NULL ){
      fprintf(stderr,"can't alloc string for :%s", path);
      free_sortfile(&filelist[arg->current_filenum]); return(-1);
   }
   memcpy((char *)sort->data_name, path+dlen, plen-dlen);
   *(sort->data_name+plen-dlen) = 0;
   if( run_number(&tmp_stat, sort, sort->data_name) ){ return(-1); }
   most_recent_calib_file(sort->data_dir, sort->run, sort->recent_cal);
    if( strlen(sort->recent_cal) == 0 && strcmp(calsrc,"file") == 0 ){
       fprintf(stdout,"No recent calib file found. Switching to ODB parameters from this .mid file instead.\n");
       sprintf(calsrc,"midas");
    }
   memset(&statbuf, 0, sizeof(struct stat)); sort->data_size = 0;
   if( sort->subrun != -1 ){ // just single subrun
      fname = subrun_filename(sort, sort->subrun);
      if( stat(fname, &statbuf) != 0 ){
         fprintf(stderr,"can't stat single subrun: %s\n", path);
      }
      sort->data_size = (long)statbuf.st_size;
   } else {                  // all subruns
      for(i=0; ; i++){
         fname = subrun_filename(sort, i);
          if( stat(fname, &statbuf) != 0 ){
             fprintf(stderr,"can't stat multi-subrun: %s\n", path); break;
          }
          sort->data_size += (long)statbuf.st_size;
      }
      sort->num_subruns = i; sort->subrun = 0;
   }
   if( (sort->histo_name = malloc(plen-dlen-ext_len+5)) == NULL ){
      fprintf(stderr,"can't alloc string for histoname");
      free_sortfile(&filelist[arg->current_filenum]); return(-1);
   }
   if( (sort->conf_name = malloc(plen-dlen-ext_len+5+1)) == NULL ){ // .json = 5bytes
      fprintf(stderr,"can't alloc string for configname");
      free_sortfile(&filelist[arg->current_filenum]); return(-1);
   }
   memcpy((char *)sort->histo_name, path+dlen, plen-dlen-ext_len);
   sprintf(sort->histo_name+plen-dlen-ext_len,".tar");
   memcpy((char *)sort->conf_name, path+dlen, plen-dlen-ext_len);
   sprintf(sort->conf_name+plen-dlen-ext_len,".json");

   if( histodir == NULL ){
      sort->histo_dir = NULL;
   } else {
      plen=strlen(histodir);
      if( (sort->histo_dir = malloc(plen+1)) == NULL ){
         fprintf(stderr,"can't alloc string for histodir");
         free_sortfile(&filelist[arg->current_filenum]); return(-1);
      }
      memcpy((char *)sort->histo_dir, histodir, plen);
      *(sort->histo_dir+plen) = 0;
       set_directory(cfg, "Histo", histodir);
   }
   if( confdir == NULL ){
      sort->conf_dir = NULL;
   } else {
      plen=strlen(confdir);
      if( (sort->conf_dir = malloc(plen+1)) == NULL ){
         fprintf(stderr,"can't alloc string for configdir");
         free_sortfile(&filelist[arg->current_filenum]); return(-1);
      }
      memcpy((char *)sort->conf_dir, confdir, plen);
      *(sort->conf_dir+plen) = 0;
      set_directory(cfg, "Config", confdir);
   }
   if( calsrc == NULL ){
      sort->cal_src = NULL;
   } else {
      plen=strlen(calsrc);
      if( (sort->cal_src = malloc(plen+1)) == NULL ){
         fprintf(stderr,"can't alloc string for cal_src");
         free_sortfile(&filelist[arg->current_filenum]); return(-1);
      }
      memcpy((char *)sort->cal_src, calsrc, plen);
      *(sort->cal_src+plen) = 0;
   }
   if( ++arg->final_filenum == FILE_QLEN ){ arg->final_filenum = 0; } //wrap
   return(0);
}

//////////////////////  sorting loop /////////////////////
// if( currentfile == finalfile ){ sleep(1); continue; }
// open_sortfilelist();
// sortnextfile()
//     create midas thread
//     read gains + init histos
//     sort data
//     close thread
//     write histos
// close_sortfilelist();
//////////////////////////////////////////////////////////
int open_next_sortfiles(Sort_status *arg)
{
   Sortfile *sort = &filelist[arg->current_filenum];
   Config *tmp_cfg;
   char ptr, tmp[256];
   if( arg->online_mode ){ return(0); } // no files in this mode

   if( sort->data_name == NULL ){ // this is a histogram summing command
      arg->sum_mode = 1;
      sum_histos(configs[0], sort); return(1); // return non-zero => dont try to sort
   } else { arg->sum_mode = 0; }

   if( sort->histo_dir == NULL ){
      sprintf(tmp, "%s/%s", sort->data_dir, sort->histo_name);
   } else {
      sprintf(tmp, "%s/%s", sort->histo_dir, sort->histo_name);
   }
   // first check if target has already been written
   if( (arg->histo_fp=fopen(tmp,"r")) != NULL ){
      fprintf(stderr,"histo file %s already exits overwriting\n", tmp);
      //sleep(1);
      fclose(arg->histo_fp);
      unlink(tmp);
   }
   // then check target is writable
   if( (arg->histo_fp=fopen(tmp,"w+")) == NULL ){ // can't write - switch to current directory
      sprintf(tmp, "./%s", sort->histo_name);
   } else {
      fclose(arg->histo_fp);
   }
   if( (arg->histo_fp=fopen(tmp,"w")) == NULL ){ // can't write ?
      fprintf(stderr,"can't open %s to write\n", tmp); return(-1);
   }
   if( sort->num_subruns == 0 ){
      sprintf(tmp, "%s/%s", sort->data_dir, sort->data_name);
   } else {
      sprintf(tmp, "%s", subrun_filename(sort, 0) );
   }
   if( (arg->data_fp=fopen(tmp,"r")) == NULL ){
      fprintf(stderr,"can't open %s to read\n", tmp);  return(-1);
   }
   fprintf(stdout,"sorting file %d %s\n", arg->current_filenum, tmp);
   if( strcmp(sort->cal_src, "file") == 0 ){
      // first subrun - open cal file or most recent
      memcpy(tmp+strlen(tmp)-8, ".json", 6);
      if( (arg->cal_fp=fopen(tmp,"r")) == NULL ){
         fprintf(stdout,"No BOR calib file - ");
         sprintf(tmp, "%s/%s", sort->data_dir, sort->recent_cal);
         if( strlen(sort->recent_cal) != 0 && (arg->cal_fp=fopen(tmp,"r")) != NULL ){
            fprintf(stdout,"using most recent: %s\n", sort->recent_cal);
         } else {
            fprintf(stdout,"No recent calib file found either. Switching to ODB parameters from this .mid file instead.\n");
         }
       }
       if( arg->cal_fp != NULL ){ fclose(arg->cal_fp);
         if( (tmp_cfg=add_config(tmp)) != NULL ){
           if( load_config(tmp_cfg, tmp, NULL) == 0 ){
             merge_configs(tmp_cfg, configs[1]);
           } else {
             fprintf(stderr,"open sortfiles: cant load config:%s\nSwitching to ODB parameters from this .mid file instead.\n", tmp);
             sprintf(sort->cal_src,"midas");
           }
           remove_config(tmp_cfg);
         } else {
           fprintf(stderr,"open sortfiles: cant create config:%s\nSwitching to ODB parameters from this .mid file instead.\n", tmp);
           sprintf(sort->cal_src,"midas");
         }
       }
     }
   arg->midas_bytes = 0;
   if( strcmp(sort->cal_src, "midas") == 0 ){
      arg->cal_overwrite = 1;
   } else {
      arg->cal_overwrite = 0;  // cal src == "config" or "file"
   }
   return(0);
}

// Only first subrun contains odb event
// *could* just sort single subrun if subrun of specified data file is nonzero
//         otherwise sort all subruns
int open_next_subrun(Sort_status *arg)
{
   Sortfile *sort = &filelist[arg->current_filenum];
   char *filename;
   fclose(arg->data_fp); arg->data_fp = NULL;

   if( sort->subrun >= sort->num_subruns-1 ){ // last or no subruns
      fprintf(stderr,"Final Subrun[%d] completed\n", sort->subrun); return(-1);
   }
   filename = subrun_filename(sort, ++sort->subrun);
   if( (arg->data_fp=fopen(filename,"r")) == NULL ){
      fprintf(stderr,"can't open %s to read\n", filename);  return(-1);
   }
   fprintf(stdout,"sorting subrun %d\n", sort->subrun);
   return(0);
}

char *subrun_filename(Sortfile *sort, int subrun)
{
   static char name[256];
   int len, digits;
   char tmp[64];

   sprintf(name, "%s/run", sort->data_dir);

   sprintf(tmp,"%d", sort->run); digits = strlen(tmp);
   len = strlen(name);
   while( digits++ < sort->run_digits ){ name[len] = '0'; name[++len+1] = 0; }
   sprintf(name+strlen(name),"%d_", sort->run);

   sprintf(tmp,"%d", subrun);  digits = strlen(tmp);
   len = strlen(name);
   while( digits++ < sort->subrun_digits ){
      name[len++]='0'; name[len] = 0;
   }
   sprintf(name+strlen(name),"%d.mid", subrun);

   return(name);
}

int close_sortfiles(Sort_status *arg)
{  // data_fp usually closed in nxtsubrun
   if( arg->online_mode ){ return(0); } // no files in this mode
   if( arg->data_fp != NULL ){ fclose(arg->data_fp); }
   fclose(arg->histo_fp);
   free_sortfile(&filelist[arg->current_filenum]);
   return(0);
}

int free_sortfile(Sortfile *sort)
{
   if( sort->histo_dir  != NULL ){ free(sort->histo_dir);  }
   if( sort->histo_name != NULL ){ free(sort->histo_name); }
   if( sort->data_dir   != NULL ){ free(sort->data_dir);   }
   if( sort->data_name  != NULL ){ free(sort->data_name);  }
   if( sort->conf_dir   != NULL ){ free(sort->conf_dir);   }
   if( sort->conf_name  != NULL ){ free(sort->conf_name);  }
   while( sort->narg > 0 && sort->arg != NULL ){
      if( sort->arg[sort->narg-1] != NULL ){ free(sort->arg[sort->narg-1]); }
      --sort->narg;
   }
   free(sort->arg);
   memset((char *)sort, 0, sizeof(Sortfile));
   return(0);
}

// read sun/subrun from filename, then count #digits in each
int run_number(Sort_status *arg, Sortfile *sort, char *name)
{
   char *ptr = name, fmt[16], tmp[256];
   FILE *fp;
   int i;

   if( strncmp(ptr, "run", 3) != 0 ){
      fprintf(stderr,"datafilename:%s does not being with \"run\"\n", name);
      return(-1);
   } ptr += 3;
   while( 1 ){
      if( !isdigit(*ptr) ){ sort->run_digits = ptr-name-3; break; }
      ++ptr;
   }
   if( *ptr != 0 ){               // name contains stuff after run number ...
      if( *ptr == '_' ){
         if( sscanf(name, "run%d_%d.mid", &sort->run, &sort->subrun) != 2 ){
            fprintf(stderr,"can't read run and subrun number in %s\n", name);
            return(-1);
         }
         sort->subrun_digits = -1; ++ptr; while( 1 ){
            if( !isdigit(*ptr) ){
               sort->subrun_digits = ptr-name-4-sort->run_digits;  break;
            }
            ++ptr;
         }
         if( sort->subrun_digits == -1 ){
            fprintf(stderr,"can't read subrun number in %s\n", name);
            return(-1);
         }
         if( strncmp(ptr, ".mid", 4) == 0 ){ return(0); }
         else {
            fprintf(stderr,"no .mid extension in datafile: %s\n", name);
            return(-1);
         }
      } else if( strncmp(ptr, ".mid", 4) == 0 ){
         fprintf(stderr,"subrun number missing in datafile: %s\n", name);
         return(-1);
      } else {
         fprintf(stderr,"bad data filename format in %s\n", name);
         return(-1);
      }
   } else {               // name only contains run number - sort all subruns
      sort->subrun = -1;
      if( sscanf(name, "run%d.mid", &sort->run) != 1 ){
         fprintf(stderr,"can't read run number in %s\n", name);
         return(-1);
      }
      for(i=1; i<5; i++){ // look for up to 5 subrun digits (usually 3)
         sprintf(fmt, "%%s/%%s_%%0%dd.mid", i);
         sprintf(tmp, fmt, sort->data_dir, name, 0 );
         if( (fp = fopen(tmp,"r")) == NULL ){ continue; }
         sort->subrun_digits = i; fclose(fp); return(0);
      }
      fprintf(stderr,"can't open subrun0 for datafile:%s\n", name);
      return(-1);
   }
   return(0);
}

int end_current_sortfile(int fd)
{
   Sort_status *arg;
   Sortfile *sort;

   arg = get_sort_status();
   arg->end_of_data = 1; //arg-> shutdown_midas = 1;
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////          Directory reading          /////////////////////
///////////////////////////////////////////////////////////////////////////
#include <dirent.h>

int most_recent_calib_file(char *data_dir, int data_run, char *result)
{
   int run, subrun, closest_run=-1, last_subrun=-1;
   struct dirent *d_ent;
   DIR *d;
   char extn[16];

   result[0]=0;
   if( (d=opendir(data_dir)) == NULL ){
     return(-1);
   }
   while( (d_ent = readdir(d)) != NULL ){
     if( strncmp(d_ent->d_name, ".", 1) == 0 ){ continue; } // Ignore
     if( sscanf(d_ent->d_name, "run%d_%d.%15s", &run, &subrun, extn)
     != 3 ){
       if( sscanf(d_ent->d_name, "run%d.%15s", &run, extn) != 2 ){
         continue; }
         subrun=-1;
       }
       if( strcmp(extn, "json") != 0 ){ continue; }
       if( run >= data_run ){ continue; } // after current file
       if( run < closest_run ){ continue; } // already seen more recent
       if( run > closest_run ){ last_subrun=-1; }
       closest_run = run;
       if( subrun >= last_subrun ){
         sprintf(result, "%s", d_ent->d_name );
         last_subrun = subrun;
       }
     }
     return( closest_run == -1 );
}

int send_datafile_list(char *path, int fd, int type)
{
   char tmp[256];  Sortfile *tmp_srt;
   int nlen, run, subrun, entry=0;
   struct dirent *d_ent;
   DIR *d;

   if( (d=opendir(path)) == NULL ){
      sprintf(tmp,"can't open directory %s\n",path);
      send_http_error_response(fd, STATUS_CODE_404,(char*)tmp);
      fprintf(stderr,"can't open directory %s\n", path);
      return(-1);
   }
   set_directory(configs[0], "Data", path);
   if( type == 1 ){
      if( (tmp_srt  = calloc(sizeof(Sortfile),    1)) == NULL ){
         fprintf(stderr,"send_datafile_list: failed alloc\n");
      }
   } else { tmp_srt = NULL; }
   send_header(fd, APP_JSON);
   put_line(fd, " [ \n", 4 );
   while( (d_ent = readdir(d)) != NULL ){
      //fprintf(stdout,"File[%s] ...\n", d_ent->d_name);
      if( strncmp(d_ent->d_name, ".", 1) == 0 ){
         continue; // Ignore
      }
      if( sscanf(d_ent->d_name, "run%d_%d.mid", &run, &subrun) != 2 ){
         if( sscanf(d_ent->d_name, "run%d.mid", &run) != 1 ){
            continue; // Not Midas Data File
         }
         subrun=0;
      }
      nlen = strlen(d_ent->d_name);
      if( strncmp(d_ent->d_name+nlen-4, ".mid", 4)     != 0 ){ // &&
          //strncmp(d_ent->d_name+nlen-7, ".mid.gz", 7)  != 0 &&
          //strncmp(d_ent->d_name+nlen-8, ".mid.bz2", 8) != 0 ){
         continue; // Not Midas DataFilename Extension
      }
      sprintf(tmp,"%s/%s", path, d_ent->d_name);
      if( stat(tmp, &statbuf) != 0 ){
         fprintf(stderr,"can't stat %s\n", tmp); statbuf.st_size = 1;
      }
      if( entry++ != 0 ){ put_line(fd, " , \n ", 5 ); }
      put_line(fd, d_ent->d_name, strlen(d_ent->d_name) );
      sprintf(tmp," , %ld ", (long)statbuf.st_size);
      put_line(fd, tmp, strlen(tmp) );
      if( (entry % 1000) == 0 ){ printf("Entry: %d\n", entry); }
      if( type == 0 ){ continue; }
      if( tmp_srt == NULL ){ put_line(fd," ,  ",4); continue; }

      if( subrun == 0 ){
         sprintf(tmp,"%s/%s", path, d_ent->d_name);
         read_datafile_info(tmp_srt, tmp);
      } else { tmp_srt->file_info[0][0] = tmp_srt->file_info[1][0] = 0; }
      if( strlen(tmp_srt->file_info[0]) > 0 ){
         sprintf(tmp," , %s ", tmp_srt->file_info[0] );
      } else {
         sprintf(tmp," , %s ", tmp_srt->file_info[1] );
      }
      put_line(fd, tmp, strlen(tmp) );
   }
   put_line(fd, " ]\n", 3 );
   if( tmp_srt != NULL ){ free(tmp_srt); }
   return(0);
}

int send_histofile_list(char *path, int fd)
{
   int nlen, run, subrun, first_entry=1;
   char tmp[128];
   struct dirent *d_ent;
   DIR *d;
   if( (d=opendir(path)) == NULL ){
      sprintf(tmp,"can't open directory, %s\n",path);
      send_http_error_response(fd, STATUS_CODE_404,(char*)tmp);
      fprintf(stderr,"can't open directory %s\n", path);
      return(-1);
   }

   send_header(fd, APP_JSON);
   put_line(fd, " [ \n", 4 );
   while( (d_ent = readdir(d)) != NULL ){
      //fprintf(stdout,"File[%s] ...\n", d_ent->d_name);
      if( strncmp(d_ent->d_name, ".", 1) == 0 ){
         continue; // Ignore
      }
      //if( sscanf(d_ent->d_name, "run%d.tar",    &run         ) != 1 &&
      //    sscanf(d_ent->d_name, "run%d_%d.tar", &run, &subrun) != 2 ){
      //   continue; // Not Midas Histogram File
      //}
      nlen = strlen(d_ent->d_name);
      if( strncmp(d_ent->d_name+nlen-4, ".tar", 4)     != 0 ){
         continue; // Not Midas HistoFilename Extension
      }
      if( first_entry == 1 ){
         first_entry = 0;
      } else {
         put_line(fd, " , \n ", 5 );
      }
      put_line(fd, d_ent->d_name, strlen(d_ent->d_name) );
   }
   put_line(fd, " ]\n", 3 );
   return(0);
}

int send_configfile_list(char *path, int fd)
{
   int nlen, run, subrun, first_entry=1;
   char tmp[128];
   struct dirent *d_ent;
   DIR *d;
   if( (d=opendir(path)) == NULL ){
      sprintf(tmp,"can't open directory, %s\n",path);
      send_http_error_response(fd, STATUS_CODE_404,(char*)tmp);
      fprintf(stderr,"can't open directory %s\n", path);
      return(-1);
   }

   send_header(fd, APP_JSON);
   put_line(fd, " [ \n", 4 );
   while( (d_ent = readdir(d)) != NULL ){
      //fprintf(stdout,"File[%s] ...\n", d_ent->d_name);
      if( strncmp(d_ent->d_name, ".", 1) == 0 ){ continue; } // Ignore
      nlen = strlen(d_ent->d_name);
      if( strncmp(d_ent->d_name+nlen-5, ".json", 5) != 0 ){
         continue; // Not Configfile Extension
      }
      if( first_entry == 1 ){
         first_entry = 0;
      } else {
         put_line(fd, " , \n ", 5 );
      }
      put_line(fd, d_ent->d_name, strlen(d_ent->d_name) );
   }
   put_line(fd, " ]\n", 3 );
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////      spectrum list and contents      ////////////////////
///////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------
// json format ...
//
// {} = object start/end, comma-sep list of zero or more "string":value
// [] = array start/end,  comma-sep list of zero or more          value
//
//     value: null, string, number[int/float], bool, array or object
//
// *NOTE* spec says NO-TRAILING-COMMAS (although they are usually accepted)
// ------------------------------------------------------------------------

int send_spectrum_list(char *cfgname, int fd)
{
   static int first_elem=1;
   Config *cfg = configs[1];
   char tmp[TITLE_LENGTH];
   int i, type, ascend;
   char *name;

     if( strlen(cfgname) != 0 ){
       for(i=0; i<MAX_CONFIGS; i++){ if( configs[i] == NULL ){ continue; }
       if(strcmp(configs[i]->name, cfgname) == 0){ cfg = configs[i]; break; }
     }
     if( i < MAX_CONFIGS ){ // Found this histogram file already open
       // If the histogram file was open, clear this config in case the file on disk has been modified
       remove_config(configs[i]);
     }
     // Now read the histogram file
     if( (cfg=read_histofile(cfgname,0)) == NULL ){
       fprintf(stderr,"send_spec_list: can't find or read:%s\n", cfgname);
       return(-1);
     }
   }
   if( (name = next_histotree_item(cfg, 1, &type, &ascend )) == NULL ){
      fprintf(stderr,"send_spectrum_list: EMPTY LIST\n");
      return(-1);
   }
   send_header(fd, APP_JSON);
   put_line(fd, "{\n", 2 );
   while( (name = next_histotree_item(cfg, 0, &type, &ascend )) != NULL ){
      if( ascend < 0 ){
         put_line(fd, " ], [ \n", 7 ); ascend = 0; first_elem=1;
      }
      while( ascend-- > 0 ){
	 first_elem=1; sprintf(tmp, " ]%s \n", ascend>0 ? " " : ", ");
         put_line(fd, tmp, strlen(tmp) );
      }
      if(        type == 0 ){ first_elem=1;       // new subfolder
	 sprintf(tmp, " \"%s\" : [ ", name);
         put_line(fd, tmp, strlen(tmp) );
      } else if (type == 1 ){                   // new histo
 	 sprintf(tmp, "%s\"%s\"\n", first_elem ? " " : ", ", name );
         put_line(fd, tmp, strlen(tmp) ); first_elem = 0;
      } else if (type == 2 ){ first_elem=1;     // exit subfolder
         put_line(fd, " ], \n", 5 );
      } else if (type == 3 ){ first_elem=1;     // new array
         sprintf(tmp, " [ \"%s\", \n", name);
         put_line(fd, tmp, strlen(tmp) );
      } else {
	fprintf(stderr,"send_spectrum_list: Unknown itemtype:%d\n", type);
      }
   }
   while( ascend-- > 0 ){
      sprintf(tmp, " ] \n");
      put_line(fd, " ] ", 3 );
   }
   put_line(fd, "}", 1 );
   return(0);
}

// next send as binary -  8bitlength 4bitType4bitWidth [only need 1+2]
// list:1bytechan,WidthbytesVal OR  WidthBytesChan[s]WidthbytesVal.....
// array is easy                    [2nd/3rd/4thpackedchan=0+>ignore]

#define HIST_HDR "{"  // 1
#define HIST_SEP ", " // 2
#define HIST_TRL "}"  // 1

static int submatrix_type[512*512];
int send_spectrum(int num, char url_args[][STRING_LEN], char *name, int fd)
{
  int i, j, k, l, m, pos, count, val, max, xbins, ybins, lastType=1, emptyCount=0, *hist_data;
  int list_max[4], list_valueSize[4], val1, val2, val3, val4, valCount, listIndex, matrixMaximum=16;
  Config *cfg = configs[1];
  char tmp[4096];  // 70 bit values is 12 characters. 12*256=3072
  char coordString[512]; // 1 byte coordinates. 1*64
  char valueString[4096]; // List: 70 bit values is 12 characters. 12*64=768. Array: 70 bit values is 12 characters. 12*256=3072
  TH1I *hist;

  int bitMask[6]={         0x3F, //  6 bits
                           0xFC0, // 12 bits
                         0x3F000, // 18 bits
                        0xFC0000, // 24 bits
                      0x3F000000, // 30 bits
                      0xC0000000, // 32 bits
                    // The variables for the matrix are currently 32 bits.
                    /*
                     0xFC0000000, // 36 bits
                   0x3F000000000, // 42 bits
                  0xFC0000000000, // 48 bits
                0x3F000000000000, // 54 bits
               0xFC0000000000000, // 60 bits
             0x3F000000000000000, // 66 bits
            0xFC0000000000000000  // 70 bits
            */
          };

  int bitShift[6]={ 0, 6, 12, 18, 24, 30};// 6 to 32 bits
  //int bitShift[12]={ 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66};// 6 to 70 bits


  if( name != NULL ){
    for(i=0; i<MAX_CONFIGS; i++){
      if( configs[i] == NULL ){ continue; }
      if( strcmp(name, configs[i]->name) == 0 ){ cfg = configs[i]; break; }
    }
    if( i == MAX_CONFIGS ){
      if( (cfg=read_histofile(name,0)) == NULL ){
        fprintf(stderr,"send_spec_list: can't find/read:%s\n", name);
        return(-1);
      }
    }
  }
  if( cfg == NULL || cfg->nhistos == 0 ){
    sprintf(tmp,"%s %s", HIST_HDR, HIST_TRL );
    send_header(fd, APP_JSON);
    put_line(fd, tmp, strlen(tmp) ); return(-1);
  }
  send_header(fd, APP_JSON);
  put_line(fd, HIST_HDR, strlen(HIST_HDR) );
  j = (name == NULL) ? 0 : 1;
  for(; j<num; j++){
    if( j > ((name == NULL) ? 0 : 1) ){
      put_line(fd, HIST_SEP, strlen(HIST_SEP) );
    }
    if( (hist = hist_querytitle(cfg, url_args[2*(j+1)+1])) == NULL ){ // don't have it
    sprintf(tmp,"\'%s\':NULL", name );
    put_line(fd, tmp, strlen(tmp) );
   } else {                            // do have this - send contents
    if( hist->data == NULL ){ read_histo_data(hist, cfg->histo_fp ); }
    if( hist->data == NULL ){ // we were not able to read data
       sprintf(tmp,"\'%s\':NULL", name );
       put_line(fd, tmp, strlen(tmp) );
    } else if( hist->type == INT_1D ){
      sprintf(tmp,"\'%s\':[", hist->title );
      put_line(fd, tmp, strlen(tmp) );
      for(i=0; i<hist->xbins; i++){
        if( i > 0){ put_line(fd, ",", 1 ); }
        sprintf(tmp,"%d", (int)hist->data[i] );
        put_line(fd, tmp, strlen(tmp) );
      }
      put_line(fd, "]", 1 );
    } else if( hist->type == INT_2D || hist->type == INT_2D_SYMM ){
      xbins = hist->xbins; if( xbins > 8192 ){ xbins = 8192; }
      ybins = hist->ybins; if( ybins > 8192 ){ ybins = 8192; }
      if( (hist_data = malloc( xbins*ybins*sizeof(int) )) == NULL){
         fprintf(stderr,"can't alloc memory for sending 2d histo\n");
         continue;
      }
         memset(hist_data, 0, xbins*ybins*sizeof(int) );
      if( hist->symm == 0 ){
         memcpy( hist_data, hist->data, xbins*ybins*sizeof(int) );
      } else { // symmetrize here before sending
         for(k=0; k<ybins; k++){
            for(i=0; i<=k; i++){
              hist_data[i+k*xbins] = hist->data[i+k*xbins] + hist->data[k+i*xbins];
            //  hist_data[k+i*xbins] = hist->data[i+k*xbins] + hist->data[k+i*xbins];
            }
         }
      }
      sprintf(tmp,"\"name\" : \"%s\",", hist->title );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"XaxisLength\" : %d,", xbins );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"XaxisMin\" : %d,", hist->xmin );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"XaxisMax\" : %d,", hist->xmax );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"YaxisLength\" : %d,", ybins  );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"YaxisMin\" : %d,", hist->ymin );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"YaxisMax\" : %d,", hist->ymax );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"symmetrized\" : %s,", hist->symm ? "true" : "false" );
      put_line(fd, tmp, strlen(tmp) );
      pos = 0;
      for(k=0; k<ybins; k+=16){
        for(i=0; i<xbins; i+=16){
          count = 0; max=0;
          for(m=0; m<16; m++){
            for(l=0; l<16; l++){
              if( hist_data[i+l + (k+m)*xbins]!=0){
                ++count;
                if( hist_data[i+l + (k+m)*xbins]>max){max = hist_data[i+l + (k+m)*xbins];}
              }
            }
          }
          // the Type will be communicated as an integer value between 0 and 8 (1 byte)
          // the Type is always element zero of each submatrix array
          // 0 = empty
          // 1 = list type using 6-bit values
          // 2 = list type using 12-bit values
          // 3 = list type using 18-bit values
          // 4 = list type using between 18 and 70-bit values
          // 5 = array type using 6-bit values
          // 6 = array type using 12-bit values
          // 7 = array type using 18-bit values
          // 8 = array type using between 18 and 70-bit values
          if( count == 0 ){ submatrix_type[pos] = 0; } // empty
          else if( count < 200 ){ // list
            submatrix_type[pos] = 1;
          }
          else { // array
            submatrix_type[pos] = 5;
          }

          // Set the submatrix type based on the maximum value
          // Set the list_valueSize[0] which will be used by Array to determine the value size.
          // The List mode calculates the list_valueSize for each string separately later.
          if(max>0){ list_valueSize[0] = 0; } // 6 bits
          if(max>0x3F){ list_valueSize[0]=1; ++submatrix_type[pos]; } // 12 bits
          if(max>0xFFF){ list_valueSize[0]=2; ++submatrix_type[pos]; } // 18 bits
          if(max>0x3FFFF){ list_valueSize[0]=3; ++submatrix_type[pos]; } // 24 bits
          if(max>0xFFFFFF){ list_valueSize[0]=4; ++submatrix_type[pos]; } // 30 bits
          if(max>0x3FFFFFFF){ list_valueSize[0]=5; ++submatrix_type[pos]; } // 32 bits
          // The variables for the matrix are currently 32 bits.
          if(max>0xFFFFFFFF){ fprintf(stdout,"Maximum value requires more than 32 bits to transmit!, %d\n",max); } // too big!
          /*
          if(max>0x3FFFFFFF){ list_valueSize[0]=5; } // 36 bits
          if(max>0xFFFFFFFFF){ list_valueSize[0]=6; } // 42 bits
          if(max>0x3FFFFFFFFFF){ list_valueSize[0]=7; } // 48 bits
          if(max>0xFFFFFFFFFFFF){ list_valueSize[0]=8; } // 54 bits
          if(max>0x3FFFFFFFFFFFFF){ list_valueSize[0]=9; } // 60 bits
          if(max>0xFFFFFFFFFFFFFFF){ list_valueSize[0]=10; } // 66 bits
          if(max>0x3FFFFFFFFFFFFFFFF){ list_valueSize[0]=11; } // 70 bits
          if(max>0xFFFFFFFFFFFFFFFFFF){ fprintf(stdout,"Maximum value requires more than 70 bits to transmit!, %d\n",max); } // Too big!
*/
          ++pos;
          if(max>matrixMaximum){ matrixMaximum = max; }
        }
      }

      sprintf(tmp,"\"ZaxisMin\" : %d,", 0 );
      put_line(fd, tmp, strlen(tmp) );
      sprintf(tmp,"\"ZaxisMax\" : %d,", matrixMaximum );
      put_line(fd, tmp, strlen(tmp) );
      put_line(fd,"\"data2\":[", 9);
      pos = 0;
      for(k=0; k<ybins; k+=16){
        for(i=0; i<xbins; i+=16){ count = 0;
          if( i+k != 0 && emptyCount==0){ put_line(fd, ",", 1); }
          if( submatrix_type[pos] == 0 ){
            ++pos; lastType=0; ++emptyCount; continue; // empty
          } else {
            if(lastType==0){
              // This is the first non-empty case so read out the empty list
              sprintf(tmp,"[0,%d],", emptyCount );
              put_line(fd, tmp, strlen(tmp) ); // empty
              emptyCount=0;
            }
            lastType=submatrix_type[pos]; // remember this submatrix_type
            if( submatrix_type[pos] < 5 ){ // list
              sprintf(tmp,"[%d,", submatrix_type[pos] );
              put_line(fd, tmp, strlen(tmp) ); // list

              // Determine the maximum values for each of the four sets of strings
              // Only make the value string as big as it needs to be for the maximum value in that subsubmatrix
              list_max[0] = list_max[1] = list_max[2] = list_max[3] = 0;
              list_valueSize[0] = list_valueSize[1] = list_valueSize[2] = list_valueSize[3] = 0;
              for(m=0; m<16; m++){
                for(l=0; l<16; l++){
                  listIndex = floor((16*m+l)/64);
                  if( hist_data[i+l + (k+m)*xbins]>list_max[listIndex]){
                    list_max[listIndex] = hist_data[i+l + (k+m)*xbins];

                    // Set the value size required for the largest number in each value string
                    if(list_max[listIndex]>0){ list_valueSize[listIndex] = 0; } // 6 bits
                    if(list_max[listIndex]>0x3F){ list_valueSize[listIndex]=1; } // 12 bits
                    if(list_max[listIndex]>0xFFF){ list_valueSize[listIndex]=2; } // 18 bits
                    if(list_max[listIndex]>0x3FFFF){ list_valueSize[listIndex]=3; } // 24 bits
                    if(list_max[listIndex]>0xFFFFFF){ list_valueSize[listIndex]=4; } // 30 bits
                    if(list_max[listIndex]>0x3FFFFFF){ list_valueSize[listIndex]=5; } // 32 bits
                    /*
                    if(list_max[listIndex]>0x3FFFFFFF){ list_valueSize[listIndex]=5; } // 36 bits
                    if(list_max[listIndex]>0xFFFFFFFFF){ list_valueSize[listIndex]=6; } // 42 bits
                    if(list_max[listIndex]>0x3FFFFFFFFFF){ list_valueSize[listIndex]=7; } // 48 bits
                    if(list_max[listIndex]>0xFFFFFFFFFFFF){ list_valueSize[listIndex]=8; } // 54 bits
                    if(list_max[listIndex]>0x3FFFFFFFFFFFFF){ list_valueSize[listIndex]=9; } // 60 bits
                    if(list_max[listIndex]>0xFFFFFFFFFFFFFFF){ list_valueSize[listIndex]=10; } // 66 bits
                    if(list_max[listIndex]>0x3FFFFFFFFFFFFFFFF){ list_valueSize[listIndex]=11; } // 70 bits
                    */
                  }
                }
              }

              // Loop through this submatrix and build the strings
              for(m=0; m<16; m++){
                for(l=0; l<16; l++){val=hist_data[i+l+(k+m)*xbins];
                  // List type (8 comma-separated strings which requires 16 quotes and 7 commas).
                  // All strings must be present (the quotes) but can have zero length.
                  // String of single-character coordinates for 0-63.
                  // String of values. All values within the string have the same size. (divide length by length of first string to give number of characters per value).
                  // These two strings are repeated four times for the four groups of 64 elements in the 256 element submatrix.
                  if((16*m+l)%64 == 0 && emptyCount==0){
                    // This is the start of a new string within the list type.
                    // Print the string of coordinates first. Then the values.
                    if((16*m+l)>0){
                      // Write out the strings from the previous quarter of this submatrix
                      strcat(coordString,"\",");
                      strcat(valueString,"\",");
                      put_line(fd, coordString, strlen(coordString) );
                      put_line(fd, valueString, strlen(valueString) );
                    }
                    // Start the new strings
                    sprintf(coordString,"\"");
                    sprintf(valueString,"\"");
                  }

                  // Skip over zero values, but check if this is the last coordinate of this submatrix
                  if( val == 0 ){
                    if((16*m+l)>254){
                      // End of this submatrix so write out the final strings
                      strcat(coordString,"\",");
                      strcat(valueString,"\"");
                      put_line(fd, coordString, strlen(coordString) );
                      put_line(fd, valueString, strlen(valueString) );
                    }
                    continue;
                  }

                  // Add this coordinate to the string
                  val1 = (16*m+l)%64 + 64;
                  if(val1 < 64 || val1 > 127){ fprintf(stdout,"List type %d: Illegal coordinate character ASCII code, %d\n",submatrix_type[pos], val1); }
                  if(val1==92){ val1=60; }
                  sprintf(tmp,"%c", val1 );
                  strcat(coordString, tmp );

                  // Add this value to the string
                  listIndex = floor((16*m+l)/64);
                  valCount = list_valueSize[listIndex];
                  while(valCount>=0){
                    val1 = ((val & bitMask[valCount])>>bitShift[valCount]) + 64;
                    if(val1 < 64 || val1 > 127){ fprintf(stdout,"Loop List type %d: Illegal character ASCII code, %d\n",submatrix_type[pos],val1); }
                    if(val1==92){ val1=60; }
                    sprintf(tmp,"%c", val1 );
                    strcat(valueString, tmp );
                    --valCount;
                  }

                  ++count;
                  if((16*m+l)>254){
                    // End of this submatrix so write out the final strings
                    strcat(coordString,"\",");
                    strcat(valueString,"\"");
                    put_line(fd, coordString, strlen(coordString) );
                    put_line(fd, valueString, strlen(valueString) );
                  }
                }// end of for l
              }// end of for m


            } else{ // array
              sprintf(valueString,"[%d,\"", submatrix_type[pos] );


              for(m=0; m<16; m++){
                for(l=0; l<16; l++){val=hist_data[i+l+(k+m)*xbins];
              // Array type
              // Encode 0-63 values as a character using the ASCII code.
              // ASCII code 0-31 are non-printable control characters (not permitted in JSON strings).
              // ASCII code 128-255 are variable definitions on different operating systems (will not use here).
              // ASCII codes 34, 47, 92 require a preceeding reverse solidus in JSON strings. (so we will avoid these)
              // Encoding shall be: val + 64 = ASCII code. So characters 64 to 127 are used with the exception of 92 which is replaced by 60.
              // Value of 92 will be made 92-32=60 instead of 92 to avoid the requirement for a preceeding character.

              // Add this value to the string
              valCount = submatrix_type[pos]-5;
              while(valCount>=0){
                val1 = ((val & bitMask[valCount])>>bitShift[valCount]) + 64;
                if(val1 < 64 || val1 > 127){ fprintf(stdout,"Loop Array type %d: Illegal character ASCII code, %d\n",submatrix_type[pos],val1); }
                if(val1==92){ val1=60; }
                sprintf(tmp,"%c", val1 );
                strcat(valueString, tmp );
                --valCount;
              }
              ++count;

            }// end of for l
          }// end of for m

          // Add the trailing quotation marks and send the string
          strcat(valueString,"\"");
          put_line(fd, valueString, strlen(valueString) );

        }


        put_line(fd, "]", 1);
        ++pos;
      }
    }
  }
  if(lastType==0){
    // This is the end so read out the empty list if its non-zero
    sprintf(tmp,"[0,%d]", emptyCount );
    put_line(fd, tmp, strlen(tmp) ); // empty
  }
  put_line(fd, "]", 1);
  free(hist_data);
    }
    }
}
put_line(fd, HIST_TRL, strlen(HIST_TRL) );
return(0);
}

#define MAX_SUBMATRICES 524287 // 19 bits
int send_binary_spectrum(int num, char url_args[][STRING_LEN], char *name, int fd)
{
  int i, j, k, l, m, pos, count, val, index, xbins, ybins, *hist_data;
  int list_max[4], list_valueSize[4], valCount;
  int num_submatrices, num_nonempty_submatrices, transfer_method, coord_size, last_type;
  int list_coordinates[4][256], list_values[4][256], list_value_counts[4], size[2];
  int list_of_indexes[4][255],list_of_values[4][255],list_of_counts[4][255];
  Config *cfg = configs[1];
  TH1I *hist;

  // The variables for the matrix are currently 32 bits.
  int bitMask_coord[4] ={ 0x0000007F, 0x00007F00, 0x007F0000, 0x7F000000 }; // 7, 15, 23, 31 bits
  int bitMask[4]       ={ 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000 }; // 8, 16, 24, 32 bits
  int bitShift[4]      ={          0,          8,         16,         24 }; // 8, 16, 24, 32 bits

  // Binary spectrum data transfer format
  // All transfered as 8-bit integers
  // Histogram Header format:
  //    "name"               // <80 characters ending with string termination
  //    "XaxisLength"        // 16 bits (<65,535)
  //    "YaxisLength"        // 16 bits (<65,535)
  //    "symmetrized"        // 1 bit, 0=non-symmeterized, 1=symmeterized
  //    "XaxisMin"           // 15 bits (<16,383)
  //    "XaxisMax"           // 16 bits (<65,535)
  //    "transfer method"    // 1 bit, 0=submatrix id and type, 1=submatrix type only header.
  //    "YaxisMin"           // 15 bits (<16,383)
  //    "YaxisMax"           // 16 bits (<65,535)
  //    "submatrix type map" // Either submatrix id and type, or submatrix type only.
  //    "submatrix data"     // Histogram contents follows
  //
  // Following the Histogram header is
  // transfer method is 1: a Submatrix type map, a 2-bit type for all submatrices of the histogram.
  // transfer method is 0: a list of the id numbers with type for non-empty Submatrices.
  // Then follows the data for each submatrix.
  //
  // Submatrix type map: is a 2-bit type per submatrix in an ordered list.
  // 0 = empty, 1 = list, 2 = array, 3 = sequence (all subsequent submatrix types the same as the previous one).
  // 4096x4096 example. 256x256 = 65,536 submatrices. 2 bits is 131,072 bits, 32,768 bytes, 16,384 characters.
  // A 'sequence' type would end the header and means all subsequent submatrices are the same as the previous type.
  //
  // Submatrix id and type map:
  // Only coordinates (id number) of non-empty submatrices appears in the header word, all others are empty type.
  // The coordinate size in bytes (1, 2, 3 or 4) is calculated from the total number of submatrices.
  // The upper 1 bit is used to encode the type. 0=list, 1=array
  // At the start of the list is the count for how many non-empty submatrices appear in the header list.
  // The count is the same number of bytes as required for each coorindate.
  //
  // Submatrix data formats by type:
  // Nothing appears for an empty type.
  //
  // List type header:
  // Four 8-bit characters represent the count of each data size values; 8, 16, 24, 32 bits
  // Coordinates given as 8-bit characters in the order of 8, 16, 24, 32 bit value sizes
  // Data values are given in order of size type (8, 16, 24, 32 bit), in the order of the coordinates given.
  // List type is smaller than array type if count is 160 which is when the submatrix is 62.5% full.
  //
  // Array type header:
  // 2 bits per entry for 256 entries. 2 bits represents type 0,1,2,3 for empty(missing),8-bit,16-bit,32-bit
  // Size of header is 512 bits = 128 bytes = 64 characters for coordinates and type.
  // This variable size of the values reduces the size of array type.
  // Data values are given in sequential order in the size given in the header. Empty values are absent.
  //

  if( name != NULL ){
    for(i=0; i<MAX_CONFIGS; i++){
      if( configs[i] == NULL ){ continue; }
      if( strcmp(name, configs[i]->name) == 0 ){ cfg = configs[i]; break; }
    }
    if( i == MAX_CONFIGS ){
      if( (cfg=read_histofile(name,0)) == NULL ){
        fprintf(stderr,"send_binary_spectrum: can't find/read:%s\n", name);
        return(-1);
      }
    }
  }
  if( cfg == NULL || cfg->nhistos == 0 ){
    send_http_error_response(fd, STATUS_CODE_404,(char*)"send_binary_spectrum can't read requested filename.");
    return(-1);
  }
  j = (name == NULL) ? 0 : 1;
  for(; j<num; j++){
    if( (hist = hist_querytitle(cfg, url_args[2*(j+1)+1])) == NULL ){ // don't have it
    send_http_error_response(fd, STATUS_CODE_404,(char*)"send_binary_spectrum can't read requested filename.");
  } else { // do have this - send contents
    if( hist->data == NULL ){ read_histo_data(hist, cfg->histo_fp ); }
    if( hist->data == NULL ){ // we were not able to read data
      send_http_error_response(fd, STATUS_CODE_404,(char*)"send_binary_spectrum can't read requested filename.");
    } else if( hist->type == INT_1D ){
      send_http_error_response(fd, STATUS_CODE_400,(char*)"send_binary_spectrum is only for 2d histograms.");
    } else if( hist->type == INT_2D || hist->type == INT_2D_SYMM ){
      xbins = hist->xbins; if( xbins > 8192 ){ xbins = 8192; }
      ybins = hist->ybins; if( ybins > 8192 ){ ybins = 8192; }
      if( (hist_data = malloc( xbins*ybins*sizeof(int) )) == NULL){
        fprintf(stderr,"can't alloc memory for sending 2d histo\n");
        send_http_error_response(fd, STATUS_CODE_500,(char*)"send_binary_spectrum can't alloc memory for sending 2d histogram.");
        continue;
      }
      memset(hist_data, 0, xbins*ybins*sizeof(int) );
      if( hist->symm == 0 ){
        memcpy( hist_data, hist->data, xbins*ybins*sizeof(int) );
      } else { // symmetrize here before sending
        for(k=0; k<ybins; k++){ for(i=0; i<=k; i++){
          hist_data[i+k*xbins] = hist->data[i+k*xbins] + hist->data[k+i*xbins];
        }}
      }
      // Send success header and start the histogram sending
      send_header(fd, APP_OCT);

      // Determine the number of submatrices that contain data.
      // This will determine which transfer method is used for the submatrix ids and types
      num_submatrices = ceil(xbins/16) * ceil(ybins/16);
        // Set maximum sizes for 524,287 submatrices (19 bits)

              // Determine the size required to store submatrix coordinates (transfer method 0)
              // Upper bit is used for the list/array type, so the coordinate is 7, 15, 23, or 31 bits
              if(      num_submatrices <= 0x80       ){  coord_size = 1; } //  7-bit values (<128 = axis lengths: 176x176)
              else if( num_submatrices <= 0x8000     ){  coord_size = 2; } // 15-bit values (256-32,767 = axis lengths: 2896x2896)
              else if( num_submatrices <= 0x800000   ){  coord_size = 3; } // 23-bit values (32,768-8,388,607 = axis lengths: 46,336x46,336)
              else if( num_submatrices <= 0x80000000 ){  coord_size = 4; } // 31-bit values (8,388,608-2,147,483,647 = axis lengths: 741,440x741,440)
              else                                    {  fprintf(stdout,"Problem with coordinate size for %d submatrices. Needs to be less than 31 bits\n",num_submatrices); }

      //  Allocate memory to the arrays, and set all entries to zero
      int *submatrix_ids = calloc(num_submatrices, sizeof(int));
      short *submatrix_type = calloc(num_submatrices, sizeof(short));
      short *submatrix_count = calloc(num_submatrices, sizeof(short));
      int8_t *binaryArray = calloc(((num_submatrices+1)*coord_size), sizeof(int8_t));

      // Determine the number of non-zero values in each submatrix.
      num_nonempty_submatrices = 0;
      pos = 0;
      for(k=0; k<ybins; k+=16){
        for(i=0; i<xbins; i+=16){ count=0;
          for(m=0; m<16; m++){ for(l=0; l<16; l++){
            if( hist_data[i+l + (k+m)*xbins]!=0){ count++; }
          }}
          if( count ==  0 ){ submatrix_type[pos] = 0; } // empty
          // List type size is Header word (8 bits), coordinates (count * 8 bits), values (count * 8, 16, 24 or 32 bits)
          // Array type size is Header word (512 bits), values (count * 8, 16, 24 or 32 bits)
          // List = 8 + (count * 8) + (list_value_counts[0]*8) + (list_value_counts[1]*16) + (list_value_counts[2]*24) + (list_value_counts[3]*32);
          // Array = 512 + ((256-count)*8) + (list_value_counts[0]*8) + (list_value_counts[1]*16) + (list_value_counts[2]*24) + (list_value_counts[3]*32);
          // So if (512-8=504). If (count*8)<(504+((256-count)*8)) then list type is smaller than array type.
          // This happens when the count is 160 which is when the submatrix is 62.5% full.
          else if( count<160 ){ submatrix_type[pos] = 1; submatrix_ids[num_nonempty_submatrices]=pos; num_nonempty_submatrices++; } // list
            else { submatrix_type[pos] = 2; submatrix_ids[num_nonempty_submatrices]=pos; num_nonempty_submatrices++; } // array
              submatrix_count[pos] = count;
              ++pos; // Advance to next submatrix
              if(pos>=num_submatrices){ break; }
            }
          }

          // Determine the transfer method for the submatrix types, and ids
          // transfer_method=0 requires the number of bytes per coordinate per non-empty submatrix
          // transfer_method=1 requires 2 bits per submatrix, so 1 byte per 4 submatrices for full histogram.
          if( num_nonempty_submatrices==0 ){ transfer_method = 0; }
          else if((num_nonempty_submatrices*coord_size) > (num_submatrices/4)){ transfer_method = 1; }
          else{ transfer_method = 0; }

          // Determine the first instance of the last submatrix type in a series and mark it
          // This will terminate the submatrix_type header at this point for transfer_method=1=submatrix type only.
          last_type = submatrix_type[pos-1]; pos-=2;
          if( submatrix_type[pos] == last_type && transfer_method==1){
            while( submatrix_type[pos] == last_type ){ pos--; if(pos==0){ break; } }
            if(pos+1<num_submatrices){ submatrix_type[pos+1] = 3; }
          }

          // Build the Histogram Header
          count = strlen(hist->title)+1; i=0;
          for(i=0; i<count; i++){
            binaryArray[i] = hist->title[i];                              // "name"            // <80 characters ending with string termination
          }
          binaryArray[i++] = (xbins & 0xFF00) >> 8;                       // "XaxisLength"     // 16 bits (<65,535), upper bits
          binaryArray[i++] = (xbins & 0xFF);                              // "XaxisLength"     // 16 bits (<65,535), lower bits
          count +=2;
          binaryArray[i++] = (ybins & 0xFF00) >> 8;                       // "YaxisLength"     // 16 bits (<65,535), upper bits
          binaryArray[i++] = (ybins & 0xFF);                              // "YaxisLength"     // 16 bits (<65,535), lower bits
          count +=2;
          binaryArray[i] = hist->symm ? 0x80 : 0;                         // "symmetrized"     // 1 bit, 0=submatrix id numbers, 1=submatrix type header.
          binaryArray[i] = binaryArray[i] | ((hist->xmin & 0x7F00) >> 8); // "XaxisMin"        // 15 bits (<16,383), upper bits
          i++;
          binaryArray[i++] = (hist->xmin & 0xFF);                         // "XaxisMin"        // 15 bits (<16,383), lower bits
          count +=2;
          binaryArray[i++] = (hist->xmax & 0xFF00) >> 8;                  // "XaxisMax"        // 16 bits (<65,535), upper bits
          binaryArray[i++] = (hist->xmax & 0xFF);                         // "XaxisMax"        // 16 bits (<65,535), lower bits
          count +=2;
          binaryArray[i] = transfer_method ? 0x80 : 0;                    // "transfer method" // 1 bit
          binaryArray[i] = binaryArray[i] | ((hist->ymin & 0x7F00) >> 8); // "YaxisMin"        // 15 bits (<16,383), upper bits
          i++;
          binaryArray[i++] = (hist->ymin & 0xFF);                         // "YaxisMin"        // 16 bits (<65,535), lower bits
          count +=2;
          binaryArray[i++] = (hist->ymax & 0xFF00) >> 8;                  // "YaxisMax"        // 16 bits (<65,535), upper bits
          binaryArray[i++] = (hist->ymax & 0xFF);                         // "YaxisMax"        // 16 bits (<65,535), lower bits
          count +=2;
          put_binary(fd, binaryArray, count ); // Send header
          // Header complete. Now process the submatrix type information.

          // DEBUG
          binaryArray[0] = 0xFF;
          binaryArray[1] = 0xFF;
          put_binary(fd, binaryArray, 2 ); // Send DEBUG
          // DEBUG

          // Submatrix type header word
          index=0;
          memset(binaryArray, 0, ((num_submatrices+1)*coord_size) );
          if( transfer_method ){
            // First handle completely empty histograms
            if(num_nonempty_submatrices==0){
              binaryArray[index] = 0x30;           // Completely empty histogram
              put_binary(fd, binaryArray, 1 ); // Send Submatrix type header word
              free(submatrix_ids); free(submatrix_type); free(submatrix_count); // Free the submatrix arrays when done
              free(binaryArray); // Free the binaryArray when done
              free(hist_data);                 // Free the memory
              return(0);                       // End now as no submatrices to transfer
            }

            // Submatrix type header is the type given as 2 bits for all submatrices
            for(i=0; i<num_submatrices; i+=4){
              binaryArray[index] = submatrix_type[i] << 6;
              binaryArray[index] = binaryArray[index] | (submatrix_type[i+1] << 4);
              binaryArray[index] = binaryArray[index] | (submatrix_type[i+2] << 2);
              binaryArray[index] = binaryArray[index] |  submatrix_type[i+3];
              index++;
              if(submatrix_type[i]  ==3){ submatrix_type[i]  =submatrix_type[i-1]; break; }
              if(submatrix_type[i+1]==3){ submatrix_type[i+1]=submatrix_type[i];   break; }
              if(submatrix_type[i+2]==3){ submatrix_type[i+2]=submatrix_type[i+1]; break; }
              if(submatrix_type[i+3]==3){ submatrix_type[i+3]=submatrix_type[i+2]; break; }
            }
            put_binary(fd, binaryArray, index ); // Send submatrix type header word
          }else{
            if(num_nonempty_submatrices==0){ continue; }
            // Submatrix type header is a list of submatrix id numbers (coorindates) with the top bit indicating the type
            // This header starts with the first entry being the number of coorindates which follow in the Submatrix type header.
            valCount=coord_size-1;
            while(valCount>=0){
              binaryArray[index] = ((num_nonempty_submatrices & bitMask[valCount])>>bitShift[valCount]);
              --valCount; index++;
            }
            // Now send each submatrix coordinate and type
            // Type is the upper 1 bit. 0=list, 1=array. The lower 7, 15, 23 or 31 bits are the submatrix id number (coordinate)
            for(i=0; i<num_nonempty_submatrices; i++){
              valCount=coord_size-1;
              binaryArray[index] = ((submatrix_type[submatrix_ids[i]]-1) << 7);
              while(valCount>=0){
                if(valCount==(coord_size-1)){ binaryArray[index] = binaryArray[index] | ((submatrix_ids[i] & bitMask_coord[valCount])>>bitShift[valCount]); }
                else{                         binaryArray[index] = binaryArray[index] | ((submatrix_ids[i] & bitMask[valCount])>>bitShift[valCount]); }
                --valCount; index++;
              }
            }
            put_binary(fd, binaryArray, index ); // Send submatrix type header word
          }
          // Header and submatrix type information complete. Now process the submatrix data.

                    // DEBUG
                    binaryArray[0] = 0xFF;
                    binaryArray[1] = 0xFF;
                    put_binary(fd, binaryArray, 2 ); // Send DEBUG
                    // DEBUG

          // Now begin transfer of the submatrices; Each submatrix begins with a header word followed by contents
          pos = 0;
          for(k=0; k<ybins; k+=16){
            for(i=0; i<xbins; i+=16){ count = 0;
              if( submatrix_type[pos] == 0 ){ ++pos; continue; // empty type
              } else {

                if( submatrix_type[pos] ==1 ){ // list
                  // Loop through this submatrix initially to count value sizes
                  memset(list_values, 0, 4*256*sizeof(int) );
                  memset(list_value_counts, 0, 4*sizeof(int) );
                  memset(list_coordinates, 0, 4*256*sizeof(int) );
                  count=0; // Use for index of coordinates
                  for(m=0; m<16; m++){
                    for(l=0; l<16; l++){val=hist_data[i+l+(k+m)*xbins];
                      if(val==0){ continue; }
                      if     ( val <= 0x100 && list_value_counts[0]<254 ){  list_values[0][list_value_counts[0]]=val; list_coordinates[0][list_value_counts[0]]=(m*16)+l; list_value_counts[0]++; count++; } //  8-bit values
                      else if( val <= 0x10000 && list_value_counts[1]<254  ){  list_values[1][list_value_counts[1]]=val; list_coordinates[1][list_value_counts[1]]=(m*16)+l; list_value_counts[1]++; count++;  } // 16-bit values
                      else if( val <= 0x1000000 && list_value_counts[2]<254 ){  list_values[2][list_value_counts[2]]=val; list_coordinates[2][list_value_counts[2]]=(m*16)+l; list_value_counts[2]++; count++;  } // 24-bit values
                      else                       {  list_values[3][list_value_counts[3]]=val; list_coordinates[3][list_value_counts[3]]=(m*16)+l; list_value_counts[3]++; count++;  } // 32-bit values
                    }
                  }
                  // Build the list submatrix header
                  index=0;
                  binaryArray[index] = list_value_counts[index];
                  if(list_value_counts[0]<254){ index++; binaryArray[index] = list_value_counts[index]; }
                  if((list_value_counts[0]+list_value_counts[1])<254){ index++; binaryArray[index] = list_value_counts[index]; }
                  if((list_value_counts[0]+list_value_counts[1]+list_value_counts[2])<254){ index++; binaryArray[index] = list_value_counts[index]; }
                  index++;
                  // Loop through this submatrix to build the coordinates. Skip empty values
                  for(m=0; m<4; m++){
                    for(l=0; l<list_value_counts[m]; l++){
                      // Add the coordinate
                      binaryArray[index] = list_coordinates[m][l];
                      index++;
                    }
                  }
                  // Loop through this submatrix to build the values, list type
                  for(m=0; m<4; m++){
                    for(l=0; l<list_value_counts[m]; l++){ val=list_values[m][l];
                      valCount = m;
                      while(valCount>=0){
                        binaryArray[index] = ((val & bitMask[valCount])>>bitShift[valCount]);
                        --valCount; index++;
                      }
                    }
                  }
                  // Send this submatrix
                  put_binary(fd, binaryArray, index ); // list type

                }else if( submatrix_type[pos] == 2 ){ // array type
                  // Build the list of data value sizes
                  // Array types, 0 (8-bit), 1 (16-bit), 2 (24-bit), 3 (32-bit)
                  // First find the maximum value per 4 subsubmatrices
                  index=0;
                  list_max[0]=list_max[1]=list_max[2]=list_max[3]=0;
                  for(m=0; m<16; m++){
                    for(l=0; l<16; l++){val=hist_data[i+l+(k+m)*xbins];
                      if(val>list_max[(int)(index/64)]){ list_max[(int)(index/64)]=val; }
                      index++;
                    }
                  }
                  // Assign value size for each subsubmatrix based on maximum value
                  for(m=0; m<4; m++){
                    if     ( list_max[m] <= 0x100     ){ list_valueSize[m]=0; } //  8-bit values
                    else if( list_max[m] <= 0x10000   ){ list_valueSize[m]=1; } // 16-bit values
                    else if( list_max[m] <= 0x1000000 ){ list_valueSize[m]=2; } // 24-bit values
                    else                               { list_valueSize[m]=3; } // 32-bit values
                  }
                  // Build array submatrix header word.
                  index=0;
                  binaryArray[index] = list_valueSize[0]<<6;
                  binaryArray[index] = binaryArray[index] | list_valueSize[1]<<4;
                  binaryArray[index] = binaryArray[index] | list_valueSize[2]<<2;
                  binaryArray[index] = binaryArray[index] | list_valueSize[3];
                  index++;

                  // Build the submatrix data values, array type
                  for(m=0; m<16; m++){
                    for(l=0; l<16; l++){val=hist_data[i+l+(k+m)*xbins];
                      valCount=list_valueSize[(int)((16*m+l)/64)];
                      while(valCount>=0){
                        binaryArray[index] = ((val & bitMask[valCount])>>bitShift[valCount]);
                        --valCount; index++;
                      }
                    }
                  }
                  // send the data for this submatrix
                  put_binary(fd, binaryArray, index ); // array type
                } ++pos; // Advance to next submatrix
              }
            }
          }
          free(submatrix_ids); free(submatrix_type); free(submatrix_count); // Free the submatrix arrays when done
          free(binaryArray); // Free the binaryArray when done
          free(hist_data); // Free the histogram memory when done
        } // End of else if( hist->type == INT_2D || hist->type == INT_2D_SYMM )
      }
    }
    return(0);
  }
