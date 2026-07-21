#include <stdio.h>
#include <stdlib.h> // alloc/free
#include <string.h>
#include <time.h>   // config last modification time
#include "config.h"
#include "grif-replay.h"
#include "histogram.h"

///////////////////////////////////////////////////////////////////////////
///////////////         Config Usersort Support        ////////////////////
///////////////////////////////////////////////////////////////////////////
// add_variable(),  find_sortvar(),
// add_global(),    remove_global()
// add_gate(),      remove_gate(), apply_gate(), unapply_gate(),
// next_condname(), add_cond(),   remove_cond(), add_cond_to_gate(),
// add_histo(), remove_histo(), find_histo(), set_directory(),
// set_midas_param(), 


/////////////////////////  Sort variables   //////////////////////////////
// NOTE: There is currently no way to calculate the value of a variable
//   given only it's text description (so this function cannot be implemented)
//int add_variable(Config *cfg, char *name, char *title)
//{
//   time_t current_time = time(NULL);
//   int i = cfg->nsortvar;
//   if( i == MAX_SORT_VARS ){
//     fprintf(stderr,"too many variables when adding %s\n", name); return(-1);
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
   for(i=0; i<cfg->ngates; i++){ // check if gate already exists
      if( strncmp(name, cfg->gatelist[i]->name, strlen(name)) == 0 &&
                    strlen(cfg->gatelist[i]->name) == strlen(name) ){
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
   histo->gatelist[i] = gate;
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
      histo->gatelist[i] = histo->gatelist[histo->num_gates-1];
   }
   histo->gatelist[histo->num_gates-1] = NULL;
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
// add user histogram, add any gates one by one after this, using above fns
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
   tmp->xvar = xvar;
   tmp->user = 1;  cfg->user_histos[cfg->nuser++] = tmp;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);
   return(0);
}

// histo->next - un-needed now have flat list?
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
   // finally remove histo from main histogram list
   cfg->histo_list[i] = cfg->histo_list[cfg->nhistos-1];//nop if i last
   cfg->histo_list[cfg->nhistos-1] = NULL;
   --cfg->nhistos;
   free(histo->data);
   cfg->folders_valid = 0;
   cfg->mtime = current_time;  save_config(cfg, DEFAULT_CONFIG, 1);

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
