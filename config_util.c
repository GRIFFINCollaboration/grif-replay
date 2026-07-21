#include <stdio.h>
#include <stdlib.h> // alloc/free
#include <string.h>
#include <time.h>   // config last modification time
#include "config.h"
#include "grif-replay.h"
#include "histogram.h"

// Histograms/Gates etc can be added at any time via the web interface
// => config can change during sorting (sort probably running most of time)
// => Each sort will need to take current config and make its own copy,
//     so any subsequent changes do not affect the in-progress sort
///////////////////////////////////////////////////////////////////////////
/////////////////          Config Commands          ///////////////////////
///////////////////////////////////////////////////////////////////////////
// init_config(), clear_config(), copy_config(), merge_config()
// add_config(), remove_config(), save_config()

// init_config() call on server startup ...
//    get most recent config and initialise all variables, gates, histos
//    if no saved config exists, just start with an empty config
int init_config(int webport)
{
   char hostname[32];
   int tmp;
   Config *cfg = configs[0];
   if( cfg == NULL ){ // not yet alloc'd live set
      if( (cfg=configs[0]=add_config("live")) == NULL ){ return(-1); }
      if( (configs[1]=add_config("sort")) == NULL ){ return(-1); }
      configs[0]->type = configs[1]->type = MEM_CONFIG;
   }
   load_config(cfg, DEFAULT_CONFIG, NULL); // attempt to load, ignore any error
   clear_calibrations(cfg); // Clear the calibrations to default values following server restart
   tmp = gethostname(hostname, 32);
   strtok(hostname, ".");
   fprintf(stdout,"Initial setup complete :-)\n\n");
   fprintf(stdout,"Now connect to grif-replay using a web browser at the following URL:\nhttps://griffincollaboration.github.io/SpectrumViewer/analyzerInterface.html?backend=%s&port=%d\n",hostname,webport);
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

   dst->nuser = 0; // add_histos takes care of these
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
