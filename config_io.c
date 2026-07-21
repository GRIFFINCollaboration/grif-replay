#include <stdio.h>
#include <string.h>
#include "config.h"
#include "grif-replay.h"
#include "histogram.h"

///////////////////////////////////////////////////////////////////////////
/////////////////            Config  I/O            ///////////////////////
///////////////////////////////////////////////////////////////////////////
// write_config(),  load_config() 
Config *configs[MAX_CONFIGS];

#define MAX_CONFIG_SIZE 1024*1024
static char config_data[MAX_CONFIG_SIZE+1];

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
            if( !isspace(tmp[i]) || instring ){
               if( len >= MAX_CONFIG_SIZE ){ --len; }
               config_data[len++] = tmp[i];
            }
         }
         if( len >= MAX_CONFIG_SIZE ){
            fprintf(stderr,"load_config: file too large, truncated - increase MAX_CONFIG_SIZE [currently %d] in config_io.c\n", MAX_CONFIG_SIZE);
            break;
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
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section or end of section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errS byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"address\":",11) != 0 ){
         fprintf(stderr,"load_config: errV1 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &address) < 1 ){
         fprintf(stderr,"load_config:errX1 byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"datatype\":",11) != 0 ){
         fprintf(stderr,"load_config: errV2 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &type) < 1 ){
         fprintf(stderr,"load_config:errX2 byte %ld\n", ptr-config_data);
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
         fprintf(stderr,"load_config: errV3 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 7; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &gain) < 1 ){
         fprintf(stderr,"load_config:errX3 byte %ld\n", ptr-config_data);
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
      } else if(strncmp(name,"GRG",3)==0){ // Only process pileup and crosstalk for HPGe
         memcpy(puk1,puk_reset, 7 * sizeof(float));
         memcpy(puk2,puk_reset, 7 * sizeof(float));
         memcpy(puE1,puE1_reset, 7 * sizeof(float));
      } else {
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
