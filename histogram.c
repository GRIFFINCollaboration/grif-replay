//    small+simple replacement for root histograms
//    Keep same code (add this pointers), style is useful when having mutiple
//    histogram types/dimensions (otherwise clearer to convert to procedural)
// file IO is separate from histogram access: analyser dumps files at endofrun
#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> // memset
#include <time.h>   // time()
#include <math.h>   // NAN
#include <zlib.h>   //

#include "config.h"
#include "histogram.h"

static int check_folder(char *folder); // return valid length (or -1) if invalid
int open_folder(Config *cfg, char* folder)
{
   char *path = cfg->current_path;
   int plen, flen;

   if( (flen = check_folder(folder)) <= 0 ){ return(-1); }
   if( flen + (plen = strlen(path)) >= HISTO_FOLDER_LENGTH ){
      fprintf(stderr,"folder path longer than %d\n", HISTO_FOLDER_LENGTH );
      return(-1);
   }
   if( path[0] == 0 ){ memcpy(path, folder, flen); path[flen]=0; }
   else { sprintf(&path[plen], "/%s", folder); }
   return(0);
}

int close_folder(Config *cfg) // work back from end to next slash
{
   char *path = cfg->current_path;
   int i, len = strlen(path);
   if( path[len-1] == '/' ){ --len; } // remove trailing /
   for(i=len-1; i>0; i--){
      if( path[i] == '/' ){ break; }
   }
   path[i] = '\0'; return(0);
}

int Zero_Histograms(Config *cfg)
{
   int i;  TH1I *ptr;
   for(i=0; i<cfg->nhistos; i++){
      ptr = (TH1I *)cfg->histo_list[i];
      ptr->Reset(ptr);
   }
   return(0);
}

/* if this starts being slow - add names to hash table */
TH1I *hist_querytitle(Config *cfg, char *name)
{
   TH1I *ptr;   int i;
   for(i=0; i<cfg->nhistos; i++){
     ptr = (TH1I *)cfg->histo_list[i];
      if( strcmp(name, ptr->title) == 0 ){ return( ptr ); }
   }
   return( NULL );
}
TH1I *hist_queryhandle(Config *cfg, char *name)
{
   TH1I *ptr;   int i;
   for(i=0; i<cfg->nhistos; i++){
     ptr = (TH1I *)cfg->histo_list[i];
      if( strcmp(name, ptr->handle) == 0 ){ return( ptr ); }
   }
   return( NULL );
}
//TH1I *hist_querynum(int num){}

char *next_histotitle(Config *cfg, int reset)
{
   static char filename[HISTO_FOLDER_LENGTH+TITLE_LENGTH+1];
   static int index;
   if( reset ){ index = 0; }
   if( index >= cfg->nhistos ){ return(NULL); }
   sprintf(filename, "%s/%s", ((TH1I *)cfg->histo_list[index])->path,
                              ((TH1I *)cfg->histo_list[index])->title );
   ++index;
   return(filename);
}

TH1I *H1_BOOK(Config *cfg, char *name, char *title, int nbins, int xmin, int xmax)
{
   int tlen, hlen; TH1I *result = &cfg->histo_array[cfg->nhistos];

   if( cfg == NULL ){ return(NULL); }
   if( cfg->nhistos >= MAX_HISTOGRAMS ){
      fprintf(stderr,"H1_BOOK: max number of histograms:%d exceeded when handling %s\n",
	 MAX_HISTOGRAMS, name );
      return(NULL);
   }
   // always allocate the data for sorting histograms
   // skip allocation for large histos read from disk (only read when needed)
   if( nbins <= SMALL_HISTO_BINS ||
       cfg == configs[0] || cfg == configs[1] ){
      if( (result->data = (int *)malloc(nbins*sizeof(int))) == NULL){
         fprintf(stderr,"H1_BOOK: data malloc failed\n");
         free(result); return(NULL);
      }
      memset(result->data, 0, nbins*sizeof(int) );
   } else {
      result->data = NULL;
   }
   if( (tlen=strlen(title)) >= TITLE_LENGTH  ){ tlen = TITLE_LENGTH-1; }
   if( (hlen=strlen(name))  >= HANDLE_LENGTH ){ hlen = HANDLE_LENGTH-1; }

   memcpy(result->path, cfg->current_path, strlen(cfg->current_path)+1 );
   memcpy(result->handle, name, hlen+1);
   memcpy(result->title, title, tlen+1);
   result->underflow     = 0;
   result->overflow      = 0;
   result->entries       = 0;
   result->xbins         = nbins;
   result->xmin          = xmin;
   result->xmax          = xmax;
   result->xrange        = xmax-xmin;
   result->valid_bins    = nbins;
   result->Reset         = &TH1I_Reset;
   result->Fill          = &TH1I_Fill;
   result->SetBinContent = &TH1I_SetBinContent;
   result->GetBinContent = &TH1I_GetBinContent;
   result->SetValidLen   = &TH1I_SetValidLen;
   result->type          = INT_1D;
   cfg->histo_list[cfg->nhistos++] = (void *)result;
   cfg->folders_valid = 0;
   return(result);
}

int TH1I_Reset(TH1I *this)
{

  // fprintf(stdout,"Reset TH1I histogram, %s\n",this->title);
   memset(this->data, 0, this->xbins*sizeof(int)); return(0);
   this->valid_bins    = this->xbins;
   this->underflow     = 0;
   this->overflow      = 0;
   this->entries       = 0;
}

int TH1I_Fill(TH1I *this, int xval, int count)
{
  float bin;
  if(this->xmin != 0 || this->xbins != this->xrange){
    bin = (1.0*xval-this->xmin) * this->xbins/(1.0*this->xrange);
  }else{
    bin = xval;
  }
//   (this->entries)++; // entries is not actually used for anything
   if( bin <            0 ){ (this->underflow)++; return(0); }
   if( bin >= this->xbins ){ (this-> overflow)++; return(0); }
   (this->data[(int)bin])+=count;
   return(0);
}

int TH1I_SetBinContent(TH1I *this, int bin, int value)
{
   if( bin < 0 || bin >= this->xbins ){ return(-1); }
   (this->data[bin])=value; return(0);
}

int TH1I_GetBinContent(TH1I *this, int bin)
{
   if( bin < 0 || bin >= this->xbins ){ return(0); }
   return( (this->data[bin]) );
}

int TH1I_SetValidLen(TH1I *this, int bins)
{
   if( bins < 0 || bins >= this->xbins ){ return(0); }
   (this->valid_bins)=bins; return(0);
}

TH2I *H2_BOOK(Config *cfg, char *name, char *title, int xbins, int xmin, int xmax, int ybins, int ymin, int ymax)
{
   int tlen, hlen; TH2I *result = (TH2I *)&cfg->histo_array[cfg->nhistos];

   if( cfg == NULL ){return(NULL); }
   if( cfg->nhistos >= MAX_HISTOGRAMS ){
      fprintf(stderr,"H2_BOOK: max number of histograms:%d exceeded when handling %s\n",
	 MAX_HISTOGRAMS, name );
      return(NULL);
   }
   if( ybins == SYMMETERIZE ){
     // This matrix will be symmetrized. Set the flag and the ybins to equal xbins.
     result->symm = 1; // Symmetric matrix
     result->type = INT_2D_SYMM;
     ybins = xbins; ymin = xmin; ymax = xmax;
   }else{
     result->symm = 0; // Non-symmetric matrix
     result->type = INT_2D;
   }
   // always allocate the data for sorting histograms
   // skip allocation for large histos read from disk (only read when needed)
   if( xbins*ybins <= SMALL_HISTO_BINS /* ||
      cfg == configs[0] || cfg == configs[1]*/ ){
      if( (result->data = (int *)malloc(xbins*ybins*sizeof(int))) == NULL){
         fprintf(stderr,"H2_BOOK: data malloc failed\n");
         return(NULL);
      }
      memset(result->data, 0, xbins*ybins*sizeof(int) );
   } else {
      result->data = NULL;
   }
   if( (tlen=strlen(title)) >= TITLE_LENGTH  ){ tlen = TITLE_LENGTH-1; }
   if( (hlen=strlen(name))  >= HANDLE_LENGTH ){ hlen = HANDLE_LENGTH-1; }

if(strncmp(title,"LBL-LBL_vs_TAC",14)==0){
fprintf(stdout,"create %s\n",title);
}

   memcpy(result->path, cfg->current_path, strlen(cfg->current_path)+1 );
   memcpy(result->handle, name, hlen+1);
   memcpy(result->title, title, tlen+1);
   result->underflow     = 0;
   result->overflow      = 0;
   result->entries       = 0;
   result->xbins         = xbins;
   result->xmin          = xmin;
   result->xmax          = xmax;
   result->ybins         = ybins;
   result->xrange        = xmax-xmin;
   result->ymin          = ymin;
   result->ymax          = ymax;
   result->yrange        = ymax-ymin;
   result->valid_bins    = xbins*ybins;
   result->Reset         = &TH2I_Reset;
   result->Fill          = &TH2I_Fill;
   result->SetBinContent = &TH2I_SetBinContent;
   result->GetBinContent = &TH2I_GetBinContent;
   result->SetValidLen   = &TH2I_SetValidLen;
   cfg->histo_list[cfg->nhistos++] = (void *)result;
   cfg->folders_valid = 0;
   return(result);
}

int TH2I_Reset(TH2I *this)
{
  // fprintf(stdout,"Reset TH1I histogram, %s\n",this->title);
   if( this->data != NULL ){
      memset(this->data, 0, this->xbins*this->ybins*sizeof(int));
   } return(0);
   this->valid_bins    = this->xbins*this->ybins;
   this->underflow     = 0;
   this->overflow      = 0;
   this->entries       = 0;
}

int TH2I_Fill(TH2I *this, int xval, int yval, int count)
{
  float xbin, ybin;
  if(this->xmin != 0 || this->xbins != this->xrange){
    xbin = (1.0*xval-this->xmin) * this->xbins/(1.0*this->xrange);
  }else{ xbin = xval; }
  if(this->ymin != 0 || this->ybins != this->yrange){
    ybin = (1.0*yval-this->ymin) * this->ybins/(1.0*this->yrange);
  }else{ ybin = yval; }
  if( xbin <            0 ){ (this->underflow)++; return(0); }
  if( xbin >= this->xbins ){ (this-> overflow)++; return(0); }
  if( ybin <            0 ){ (this->underflow)++; return(0); }
  if( ybin >= this->ybins ){ (this-> overflow)++; return(0); }
  if( this->data == NULL ){
     if( (this->data=(int *)malloc(this->xbins*this->ybins*sizeof(int)))==NULL){
        fprintf(stderr,"TH2I_Fill: data malloc failed for %s\n",this->handle);
        return(-1);
     }
     memset(this->data, 0, this->xbins*this->ybins*sizeof(int) );
  }
  (this->data[(int)xbin + (int)(ybin)*this->xbins])+=count;
  return(0);
}

int TH2I_SetBinContent(TH2I *this, int xbin, int ybin, int value)
{
   if( xbin < 0 || xbin >= this->xbins ){ return(-1); }
   if( ybin < 0 || ybin >= this->ybins ){ return(-1); }
   if( this->data == NULL ){
      if( (this->data=(int *)malloc(this->xbins*this->ybins*sizeof(int)))==NULL){
         fprintf(stderr,"TH2I_Fill: data malloc failed for %s\n",this->handle);
         return(-1);
      }
      memset(this->data, 0, this->xbins*this->ybins*sizeof(int) );
   }
   (this->data[xbin+ybin*this->xbins])=value; return(0);
}

int TH2I_GetBinContent(TH2I *this, int xbin, int ybin)
{
   if( xbin < 0 || xbin >= this->xbins ){ return(-1); }
   if( ybin < 0 || ybin >= this->ybins ){ return(-1); }
   if( this->data == NULL ){ return(0); }
   return( (this->data[xbin+ybin*this->xbins]) );
}

int TH2I_SetValidLen(TH2I *this, int bins)
{
   if( bins < 0 || bins >= this->xbins*this->ybins ){ return(0); }
   (this->valid_bins)=bins; return(0);
}

/////////////////////////////////////////////////////////////////////////////
//////////////////////     histogram file IO etc.  //////////////////////////
// ----------------------------------------------------------------------------
// if not using root, need another way of storing histograms on disk
//    do not want old method of 1 file per histogram
//    want multi-histogram archive format
//    tar,zip,7zip,rar,cpio[=rpm],...,... don't like any except tar
// histo file format - tar?
//   can chain-seek to read single files
//   could also add index file at end - not really needed though
//   **by chain-seeking to eof to get index, will already have obtained it!!
//   **if uncompressed** can live update disk file
//      even if compressed - can add new entry at end, and mark old invalid
//         will have to sort out when closing file - not that important
//         so could defer this (e.g. if analyzer crashes)
// compression (not only gzip etc, use cubesort methods)
//    [could use more working code for compressing spartan config file!]
//    empty spectra - zero length!
//    sparse - store chan/count pairs
//    rle - for non-sparse but mainly same value
//    *maybe* finally gzip if full (and large)
//       - not worth it if small - will gain nothing
//       - may not gain enough even if big
//    (large spectra need blocks with different methods for each block)
// histo format ...
//    series of bins, over+underflow bins are also useful
//    axis-scale different to bin-scale seemed useful but almost never used
//    int16/int32/float32, 1d or 2d, don't bother with 3+d
//griffin histogram xbins=x_xxx_xxx ybins=x_xxx_xxx comp=0 binformat=float32\n
//             V2                                        1           int16__\n
//                                                                   ascii3col
// first file "griffin histogram archive file"[zero length]
//   other fields for communication between programs having file open
//   memory map this entry - doesn't actually reduce disk-reads when checking?
//     could do similar with unused bits of indiv. spec headers (unnecessary?)
//
// use name[100] for [path]/handle [=>short name in tar content listings]
// use uid/gid for xbin/ybin - also shown in listing
// use real size and mtime, linkflag=0 to allow extract file, title->linkname
// encode bin-format in mode [shown in listing], but also store name in ownr
// store compression format in group, over/underflow in devmaj/min?
//
// have tree since store pathnames, also don't need to store directory entries
// as histo folders don't have any attributes that need saving
//
//    **also need handle, title, maybe calibration etc? (probably not)
//    **this plus mtime etc can be stored in tar header**
//    file itself has owner/group - could use these fields for bins
// extension .grif not .tar - to distinguish
// tar headers ...     **size OCTAL**    linkname[100] follows type -> 257bytes
//    name[100] mode[8] uid[8] gid[8] size[12] mtime[12] cksum[8] linktype[1]
//    at offset=257 are extension fields (ignored by earlier versions)
//  "ustar "[6] version[2] owner[32] group[32] devmaj[6] devmin[6] prefix[155]
//  [which leaves 12 bytes of padding]
//     checksum is sum of header values [0-255]x512 with space[32] for cksum
// tar extensions ...
//    gnu tar - numeric fields can be bin not ascii [set msb of leading byte]
//    prefix[155] -> atime[12] ctime[12] offset[12] longname[4] pad[1]
//    4 sparse entries:{offset[12] bytes[12]} realsize[12] extended[1] pad[17]
//
// tar includes plenty of unused header space plus standard extension method
// - so can include everything in tar header
// => minimum size is 0.5k (300k per 600 histo) (1-2Mbyte per tigress run!)
// tigress eta files have 2600 spectra in 1150k (would be ~*2)
//
// will be trivial to convert these files to root files
// (using script which will work on any root version, not compiled program)
// ----------------------------------------------------------------------------
// diskfiles: just array of spectra[each with own header] (no file-header)

typedef struct file_header_struct { // should be 512 bytes
   char name[100];  char    mode[8];  char      uid[8];  char    gid[8];
   char  size[12];  char  mtime[12];  char    cksum[8];  char   type[1];
   char link[100];
   char  magic[6];  char version[2];  char   owner[32];  char group[32];
   char devmaj[8];  char  devmin[8];  char prefix[155];  char   pad[12];
} File_head;    // uid/gid:xbins/ybins
                // own/group:xmin_xmax  ymin_ymax

//typedef struct histo_header_struct { // soon will have to use this
//   int sbins; int xmin; int xmax;  int ybins; int ymin; int ymax;
//   int overflow;  int underflow;   int spare[8];
//} Histo_head;

static File_head file_head;
static char      file_body[FILE_BUFSIZ];

int delete_histograms(Config *cfg)
{
   TH1I *ptr;
   int i;
   for(i=0; i<cfg->nhistos; i++){
      ptr = (TH1I *)cfg->histo_list[i];
      free(ptr->data); //free(ptr);
   }
   cfg->nhistos = 0;
   return(0);
}

int write_histofile(Config *cfg, FILE *fp)
{
   int i, type, size, pad, cksum;

   if( fp == NULL ){
      fprintf(stderr,"*** No file open in write_histofile\n"); return(-1);
   }
   if( cfg == NULL ){ // not yet alloc'd live set
      fprintf(stderr,"No histos defined\n"); return(-1);
   }
   fseek(fp, 512, SEEK_SET);
   write_config(cfg, fp);
   size = (int)ftell(fp) - 512;
   pad = (512 - (size % 512)) % 512; // pad with 0 to multiple of 512 bytes
   memset(&file_head, 0, sizeof(file_head));
   sprintf(file_head.name,"config_file"); sprintf(file_head.mode,"0000755");
   sprintf(file_head.uid,"%07o", 0);      sprintf(file_head.gid,"%07o", 0);
   sprintf(file_head.mtime,"%011o",0);    memset(file_head.cksum,' ', 8);
   sprintf(file_head.owner ,"X");         sprintf(file_head.group ,"X");
   file_head.type[0] = 0;                 sprintf(file_head.magic,"ustar");
   file_head.version[0] = file_head.version[1] = '0';
   sprintf(file_head.size    ,"%011o", size );
   cksum = 0; for(i=0; i<512; i++){ cksum += file_head.name[i]; }
   sprintf(file_head.cksum   ,"%07o", cksum );
   fseek(fp, 0, SEEK_SET);
   fwrite( &file_head, sizeof(File_head), 1, fp);
   fseek(fp, size+pad+512, SEEK_SET);
   for(i=0; i<cfg->nhistos; i++){
      if( ((TH1I *)cfg->histo_list[i])->suppress ){ continue; }
      switch( type = ((TH1I *)cfg->histo_list[i])->type ){
      case INT32_1D: write_th1I(fp, cfg->histo_list[i] ); break;
      case INT32_2D: write_th1I(fp, cfg->histo_list[i] ); break;
      case INT32_2D_SYMM: write_th1I(fp, cfg->histo_list[i] ); break;
      }
   }
   return(0);
}

int close_histofile(Config *cfg)
{
   return( remove_config(cfg) );
}

///////////////////////////////////////////////////////////////////////////
#define COMPRESS_BUFSIZ 270000000
// allow for 8k^2 = 65M-chan * 4bytes = 256Mbytes [~270000000]
static char compress_buf[COMPRESS_BUFSIZ];
int compress_buffer(char *input, int size)
{
   z_stream strm;
   int status;

   strm.zalloc = Z_NULL;
   strm.zfree = Z_NULL;
   strm.opaque = Z_NULL;
   if( (status = deflateInit(&strm, Z_DEFAULT_COMPRESSION)) != Z_OK ){
      fprintf(stderr,"zlib compression error\n"); return(-1);
   }
   strm.avail_in = size;
   strm.next_in = (unsigned char *)input;

   strm.avail_out = COMPRESS_BUFSIZ; // should never be filled
   strm.next_out = (unsigned char *)compress_buf;
   deflate(&strm, Z_FINISH);
   status = COMPRESS_BUFSIZ - strm.avail_out;
   deflateEnd(&strm);
   return(status);
}
int decompress_buffer(char *input, int size)
{
   z_stream strm;
   int status;

   strm.zalloc = Z_NULL;
   strm.zfree = Z_NULL;
   strm.opaque = Z_NULL;
   strm.avail_in = 0;strm.next_in = Z_NULL;
   if( (status = inflateInit(&strm)) != Z_OK ){
      fprintf(stderr,"zlib decompression error\n"); return(-1);
   }
   strm.avail_in = size;
   strm.next_in = (unsigned char *)input;

   strm.avail_out = COMPRESS_BUFSIZ; // should never be filled
   strm.next_out = (unsigned char *)compress_buf;
   inflate(&strm, Z_NO_FLUSH);
   status = COMPRESS_BUFSIZ - strm.avail_out;
   inflateEnd(&strm);
   return(status);
}
int read_histo_data(Histogram *histo, FILE *fp)
{
   int bins = (histo->ybins != 0) ? histo->xbins*histo->ybins : histo->xbins;
   int size = histo->data_size;
   if( histo->data != NULL ){
      fprintf(stderr,"histo data already present:%s\n", histo->title );
      return(-1);
   }
   if( (histo->data = (int *)malloc(bins*sizeof(int))) == NULL){
      fprintf(stderr,"read_histo_data: data malloc failed\n");
      return(-1);
   }
   memset(histo->data, 0, bins*sizeof(int));
   if( fseek(fp, histo->file_data_offset, SEEK_SET) < 0 ){
      fprintf(stderr,"failed_seek histo:%s\n", histo->title );
      return(-1);
   }
   if( fread( &file_body, sizeof(char), size, fp) < size ){
      fprintf(stderr,"short read histo:%s[%d]\n", histo->title, size );
      return(-1);
   }
   // compressed data starts with signed bytes: 120,-100
   if( bins > SMALL_HISTO_BINS && file_body[0]==120 && file_body[1]==-100 ){
      size = decompress_buffer(file_body, size);
      memcpy(histo->data, compress_buf, size);
   } else {
      memcpy(histo->data, file_body, size);
   }
   return(0);
}
///////////////////////////////////////////////////////////////////////////
// alloc and read new histogram set (+config file)
Config *read_histofile(char *filename, int config_only)
{
   unsigned size, xbins, xmin, xmax, ybins, ymin, ymax, len;
   long file_offset = 0;
   int pad, err=0, bins;
   Config *cfg;
   TH1I *histo;
   char tmp[64];
   FILE *fp;

   if( (cfg=add_config(filename)) == NULL ){ return(NULL); }
   if( (fp=cfg->histo_fp=fopen(filename,"r")) == NULL){ // can't open
      fprintf(stderr,"can't open file:%s to read\n", filename );
      remove_config(cfg); return(NULL);
   }
   while( 1 ){
      if( fread( &file_head, sizeof(File_head), 1, fp) < 1 ){ break; }
      memcpy(tmp, file_head.size, 12); tmp[12]=0;
      if( sscanf(tmp, "%o", &size) < 1 ){
         fprintf(stderr,"can't read histo size from:%s\n", file_head.size);
         err=1; break;
      }
      file_offset += sizeof(File_head);
      pad = (512 - (size % 512)) % 512; // pad to multiple of 512 bytes
      if( cfg->nhistos == 0 && strcmp(file_head.name, "config_file") == 0 ){
         if( config_only ){
            if( fread( &file_body, sizeof(char), size+pad, fp) <  size+pad ){
               fprintf(stderr,"short read histo:%s[%d]\n",file_head.name,size);
               err=1; break;
            }
            load_config(cfg, NULL, file_body );
            return( cfg );
         } else { // do not try to load this as histogram - skip it
            if( fseek(fp, size+pad, SEEK_CUR) < 0 ){
               fprintf(stderr,"short seek histo:%s[%d]\n",file_head.name,size);
               err=1; break;
            }
            file_offset += size+pad;
            continue;
         }
      }
      if( config_only ){ // config file was not the first entry
         remove_config( cfg ); return( NULL );
      }
      memcpy(tmp, file_head.uid, 8); tmp[8]=0;
      if( sscanf(tmp, "%o", &xbins) < 1 ){
         fprintf(stderr,"can't read histo xbins from:%s\n", file_head.uid);
         err=1; break;
      }
      memcpy(tmp, file_head.gid, 8); tmp[8]=0;
      if( sscanf(tmp, "%o", &ybins) < 1 ){
         fprintf(stderr,"can't read histo ybins from:%s\n", file_head.gid);
         err=1; break;
      }
     if( sscanf(file_head.owner,"%d_%d", &xmin, &xmax) != 2 ){
         xmin = 0; xmax = xbins;
      }
      if( sscanf(file_head.group,"%d_%d", &ymin, &ymax) != 2 ){
         ymin = 0; ymax = ybins;
      }
      len = strlen(file_head.prefix);
      memcpy(cfg->current_path, file_head.prefix, len);
      cfg->current_path[len]=0;
      bins = ( ybins == 0 ) ? xbins : xbins*ybins;
      if( bins <= SMALL_HISTO_BINS*sizeof(int) ){
         if( fread( &file_body, sizeof(char), size+pad, fp) <  size+pad ){
            fprintf(stderr,"short read histo:%s[%d]\n",file_head.name,size);
            err=1; break;
         }
      } else {
         if( fseek(fp, size+pad, SEEK_CUR) < 0 ){
            fprintf(stderr,"short seek histo:%s[%d]\n", file_head.name, size);
            err=1; break;
         }
      }
      if( ybins == 0 ){
         histo = (TH1I *)H1_BOOK(cfg, file_head.name,file_head.link,xbins,0,xbins);
      } else {
        if(file_head.type[0] == 'C'){  ybins=SYMMETERIZE; } // set ybins to 0 for H2_BOOK to handle as symmetric
         histo = (TH1I *)H2_BOOK(cfg, file_head.name,file_head.link,xbins,0,xbins,ybins,0,ybins);
      }
      if( bins <= SMALL_HISTO_BINS ){
         memcpy(histo->data, file_body, size); // Do not include pad!
      } else {
         histo->file_data_offset = file_offset;  histo->data_size = size+pad;
      }
      file_offset += size+pad;
   }
   if( err ){ remove_config(cfg); return(NULL); }
   return(cfg);
}

int write_th1I(FILE *fp, void *ptr)
{
   int i, cksum, size, pad, count, mode, bins;
   TH1I *hist = (TH1I *)ptr;
   time_t filetime;
   time(&filetime);

   memset(&file_head, 0, sizeof(file_head));  file_head.cksum[7] = ' ';
   sprintf(file_head.name,"%s", hist->handle);
   //else { sprintf(file_head.name,"%s/%s", hist->path, hist->handle); }
   //sprintf(file_head.prefix,"%s", hist->title);
   sprintf(file_head.link,"%s", hist->title);
   sprintf(file_head.prefix,"%s", hist->path);
   sprintf(file_head.mode    ,"0000755"); // encode type?
   sprintf(file_head.uid     ,"%07o", hist->xbins);
   sprintf(file_head.gid     ,"%07o", (hist->type==INT_2D || hist->type==INT_2D_SYMM) ? hist->ybins : 0);
   //sprintf(file_head.size    ,"%011o", 0 ); // fill in proper size later
   sprintf(file_head.mtime   ,"%011lo", filetime );
   memset(file_head.cksum, ' ', 8);// cksum entry counted as 8 blanks, no null

   switch( hist->type ){
   case INT32_1D: file_head.type[0] = 'A'; break;
   case INT32_2D: file_head.type[0] = 'B'; break;
   case INT32_2D_SYMM: file_head.type[0] = 'C'; break;
   }

   // linkname[100] is from 157 to 256
   //sprintf(file_head.magic   ,"Grif1"    ); // ustar
   sprintf(file_head.magic   ,"ustar"    ); // ustar FOLLOWED BY NULL
   file_head.version[0] = file_head.version[1] = '0';
   sprintf(file_head.owner   ,"%d_%d", hist->xmin, hist->xmax );
   sprintf(file_head.group   ,"%d_%d",  hist->ymin, hist->ymax );
   sprintf(file_head.devmaj  ,""         );
   sprintf(file_head.devmin  ,""         );

   // binformat, #entries, compression format(use size for now)

   bins = (hist->type==INT_2D || hist->type==INT_2D_SYMM) ? hist->xbins*hist->ybins : hist->xbins;
   // check for empty (count non-zero at same time)
   count = 0; if( hist->data != NULL ){
      for(i=0; i<bins; i++){
         if( hist->data[i] != 0 ){
            ++count; if( hist->type==INT_2D || hist->type==INT_2D_SYMM){ break; }
         }
      }
   }
   if( count == 0 ){ size=0; mode=0; }
   //else if( count * 6 < ptr->xbins*sizeof(float) ){ size = count*6; mode=1; }
   else if( bins > 65536 ){
      size = compress_buffer((char *)(hist->data), bins*sizeof(int));
      mode = 3; // gzip compressed
   } else { size = bins*sizeof(int);  mode=2; } // just write data

   sprintf(file_head.size    ,"%011o", size );
   sprintf(file_head.owner   ,"%10d",  hist->underflow );
   sprintf(file_head.group   ,"%10d",  hist->overflow  );

   cksum = 0; for(i=0; i<512; i++){ cksum += file_head.name[i]; }
   sprintf(file_head.cksum   ,"%07o", cksum );
   fwrite( &file_head, sizeof(File_head), 1, fp);

   pad = (512 - (size % 512)) % 512; // pad with 0 to multiple of 512 bytes
   memset(file_body+size, 0, pad);
   switch(mode){
   case  0: break;
   case  1: size = pad = 0; break;
   case  2: memcpy( file_body, hist->data, size ); break;
   case  3: memcpy( file_body, compress_buf, size ); break;
   default: break;
   }
   fwrite( &file_body, sizeof(char), size+pad, fp);

   // run test and make file, do root script to view spectra

   return(0);
}

int sum_th1I(Config *dst_cfg, Config *src_cfg, TH1I *src)
{
   int i, bins, ybins;
   TH1I *dst;

   if( (dst=find_histo(dst_cfg, src->handle)) == NULL ){
      memcpy(dst_cfg->current_path, src->path, HISTO_FOLDER_LENGTH);
      if( src->type == INT_1D ){
         dst = H1_BOOK(dst_cfg, src->handle, src->title, src->xbins, src->xmin, src->xmax);
      } else {
        // Handle symmetrized and non-symmetrized matrices
        if( src->symm == 1 ){
          ybins = SYMMETERIZE;
        }else{
          ybins = src->ybins;
        }
         dst = (TH1I *)H2_BOOK(dst_cfg, src->handle, src->title, src->xbins, src->xmin, src->xmax, ybins, src->ymin, src->ymax);
      }
      if( dst == NULL ){ return(-1); }
      if( dst->data == NULL ){
         bins = (dst->ybins != 0) ? dst->xbins*dst->ybins : dst->xbins;
         if( (dst->data = (int *)malloc(bins*sizeof(int))) == NULL){
            fprintf(stderr,"sum_TH1I: data malloc failed\n");
            return(-1);
         }
      }
      memcpy(dst->data, src->data, dst->valid_bins*sizeof(int) );
      return(0);
   }
   for(i=0; i<dst->valid_bins && i<src->valid_bins; i++){ dst->data[i] += src->data[i]; }
   return(0);
}

int old_sum_th1I(Config *dst_cfg, TH1I *dst, TH1I *src)
{
   int i, bins, ybins;
   if( dst == NULL ){
      memcpy(dst_cfg->current_path, src->path, HISTO_FOLDER_LENGTH);
      if( src->type == INT_1D ){
         dst = H1_BOOK(dst_cfg, src->handle, src->title, src->xbins, src->xmin, src->xmax);
      } else {
        // Handle symmetrized and non-symmetrized matrices
        if( src->symm == 1 ){
          ybins = SYMMETERIZE;
        }else{
          ybins = src->ybins;
        }
         dst = (TH1I *)H2_BOOK(dst_cfg, src->handle, src->title, src->xbins, src->xmin, src->xmax, ybins, src->ymin, src->ymax);
      }
      if( dst == NULL ){ return(-1); }
      if( dst->data == NULL ){
         bins = (dst->ybins != 0) ? dst->xbins*dst->ybins : dst->xbins;
         if( (dst->data = (int *)malloc(bins*sizeof(int))) == NULL){
            fprintf(stderr,"sum_TH1I: data malloc failed\n");
            return(-1);
         }
      }
      memcmp(dst->data, src->data, dst->valid_bins*sizeof(int) );
      return(0);
   }
   for(i=0; i<dst->valid_bins && i<src->valid_bins; i++){ dst->data[i] += src->data[i]; }
   return(0);
}

////////////////////////////////////////////////////////////////////////////
///////////////////////    FOLDER HANDLING  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// max length ~256, no odd characters or /
int check_folder(char *folder) // return valid length (or -1) if invalid
{
   int i, len;

   for(i=0; i<HISTO_FOLDER_LENGTH; i++){
      if( folder[i] == '\0' ){ break; }
   }
   if( (len=i) >= HISTO_FOLDER_LENGTH ){
      fprintf(stderr,"folder longer than %d\n", HISTO_FOLDER_LENGTH );
      return(-1);
   }
   if( folder[len-1] == '/' ){ --len; } // remove trailing /
   for(i=0; i<len; i++){ // check input - no '/'
      if( folder[i] == '/' ){ break; }
      if( folder[i]  < ' ' ){ break; } // space = 32
      if( folder[i]  > '~' ){ break; } //     ~ = 127
   }
   if( i != len ){
      fprintf(stderr,"invalid characters in folder:%s\n", folder );
      return(-1);
   }
   return(len);
}

// could be a little simpler to do this as histograms/folders are created
//    [NOTE going backwards (close_folder) is awkward in singly-linked list]
//    [     additional variables are required to store reverse path ...    ]
//    [  e.g.: list of pointers (to current folder and each higer folder)  ]
// BUT, doing after histos are all created is also simple enough
// ---------------- Example Tree ----------------------------
//     TOP [--NO]                 Top has no Name/Next-Folder
//      |
//      A------B--C
//      |      |  |
//      D--E   F  G--H
//         |   |  |  |
//         I   J  K  L
// ----------------------------------------------------------
// generate nested list from unsorted histogram list
// depth first, otherwise keep order same as when created
int create_histo_tree(Config *cfg)
{
   char name[HISTO_FOLDER_LENGTH];
   char *fldstart, *p;
   Folder *curr_folder, *f;
   TH1I *curr_histo;
   int i, len;

   if( cfg->folders_valid ){ return(0); } // already up to date
   else {
      delete_histo_tree(cfg);             // start over
      cfg->folders_valid = 1;
   }
   for(i=0; i<cfg->nhistos; i++){
      if( ((TH1I *)cfg->histo_list[i])->suppress ){ continue; }
      // split path into folders, find position in folder tree, and add histo
      p = ((TH1I *)cfg->histo_list[i])->path;
      curr_folder = &cfg->first_folder;
      while(1){ // loop over path ..
	 if( *p == '/' ){ ++p; } // skip over any folder separator
	 fldstart = p;
         while( *p != '/' && *p != 0 ){ ++p; }
         if( (len = p - fldstart) == 0 ){ // NULL name => reached EndOfPath
            if( (curr_histo = curr_folder->first_histo) == NULL ){
	       curr_folder->first_histo = cfg->histo_list[i]; break;//DONE
   	    }
            while( curr_histo->next != NULL ){ curr_histo = curr_histo->next; }
            curr_histo->next = cfg->histo_list[i]; break;           //DONE
         } else {                             // descend into this folder
	    memcpy(name, fldstart, len); name[len] = 0;
            if( curr_folder->next_subfolder == NULL ){  // no subfolders exist
	       if((f=curr_folder->next_subfolder=malloc(sizeof(Folder)))==NULL){
                  fprintf(stderr,"folder_tree: structure malloc failed\n");
	          return(-1);
	       }
               memcpy(f->name, name, len+1);
               f->next_subfolder = f->next_folder = NULL;
               f->first_histo = NULL;
               curr_folder = f;
	    } else {                       // subfolders do exist - find or add
	       curr_folder = curr_folder->next_subfolder;
               while( strncmp(name, curr_folder->name, len) !=   0 ||
		                  strlen(curr_folder->name) != len ){
	          if( curr_folder->next_folder != NULL ){
 		     curr_folder = curr_folder->next_folder; continue;
		  }
                  // did not find - add new "next_folder"
                  if((f=curr_folder->next_folder=malloc(sizeof(Folder)))==NULL){
                     fprintf(stderr,"folder_tree: structure malloc failed\n");
	             return(-1);
	          }
                  memcpy(f->name, name, len+1);
                  f->next_subfolder = f->next_folder = NULL;
                  f->first_histo = NULL;
                  curr_folder = f;
		  break;
               }
	    }
	 } // done descent - continue along path
      } // done current histogram path
   } // done all histograms
   return(0);
}

// depth first - reuse current_path for deletion
// depth=0 is for "first_folder" with no name/nxtfolder, and part of set
// (do not do any freeing on depth=0)
int delete_histo_tree(Config *cfg)
{
   Folder *tmp, *folder = cfg->treepath[0] = &cfg->first_folder;
   int depth = 0;

   while( 1 ){
      if( folder->next_subfolder != NULL ){
         folder = cfg->treepath[++depth] = folder->next_subfolder;
      }
      // no next_subfolder from here - free current, and go to nxt folder
      if( depth == 0 ){ break; }
      tmp = folder->next_folder;
      free(folder->name); free(folder);
      if( (folder = cfg->treepath[depth] = tmp) == NULL ){ // no nxt - go up
         folder = cfg->treepath[--depth];
      }
   }
   return(0);
}

// output histogram tree item-by-item
// types ... 0:EnterFolder    [return folder name]
//           1:AddHisto       [return histo title]
//           2:ExitFolder     [return zero-length string]
// Null return => done [type=2 - exit top folder]
// first call [enter top folder, with empty name] => return empty string, type0
#define MAX_DEPTH HISTO_FOLDER_LENGTH
static int curr_depth;
char *next_histotree_item(Config *cfg, int reset, int *type, int *ascend)
{
   static char name[HISTO_FOLDER_LENGTH];
   static Folder* folder;
   static TH1I *histo;

   *ascend = 0;
   if( reset ){
      create_histo_tree(cfg);
      folder = cfg->treepath[0] = &cfg->first_folder;
      cfg->current_depth = 0;
      histo = NULL;
      name[0] = 0;
      *type = 0; return(name);
   }
   // depth first, continue till next_subfolder is null
   if( folder->next_subfolder != NULL ){
      folder = cfg->treepath[++cfg->current_depth] = folder->next_subfolder;
      histo = NULL;
      memcpy(name,  folder->name, strlen(folder->name)+1);
   // To use restricted format where 1st element of array is final folder name
   // Need to check next_subfolder, and if null, set type = 3 instead of 0
   // This START AN ARRAY AND make the folder name an array element
   //                                               not a sub-object
      if( folder->next_subfolder == NULL ){
         *type = 3; return(name);
      }
      *type = 0; return(name);
   }
   // histo starts NULL, list ends when curr_histo->next_histo is NULL
   if( histo == NULL ){
      histo = folder->first_histo;
   } else {
      histo = histo->next;
   }
   if( histo == NULL ){ // end of histos in this folder - exit folder
      // **before going up a level, need to go to next folder on this level
      if( folder->next_folder != NULL ){
         folder = cfg->treepath[cfg->current_depth] = folder->next_folder;
         histo = NULL;
         memcpy(name,  folder->name, strlen(folder->name)+1);
         if( folder->next_subfolder == NULL ){
            *ascend = -1; *type = 1; return(name);
         }
         *ascend = 1; *type = 0; return(name);
      } else {
      // **when finally do go up a level, go to *next* folder on that level
      //                           (or will repeat what has just been done)
	 *ascend = 1;
 	 while( --cfg->current_depth > 0 ){ ++*ascend;
	    if( (folder = cfg->treepath[cfg->current_depth]->next_folder) != NULL ){
               cfg->treepath[cfg->current_depth] = folder;
               histo = NULL;
               memcpy(name,  folder->name, strlen(folder->name)+1);
               if( folder->next_subfolder == NULL ){
                  *type = 3; return(name);
               }
               *type = 0; return(name);
	    }
	 }
         *type = 2; return(NULL);  // all done (curr_depth = 0)
      }
   }
   if( histo->type == INT_2D  || histo->type==INT_2D_SYMM){
      sprintf(name, "%s:2d", histo->title );
   } else {
      sprintf(name, "%s", histo->title );
   }
   *type = 1; return(name);
}
// ----------------------------------------------------------------------------
// json format ...
//
// {} = object start/end, comma-sep list of zero or more "string":value
// [] = array start/end,  comma-sep list of zero or more          value
//
//     value: null, string, number[int/float], bool, array or object
//
// *NOTE* spec says NO-TRAILING-COMMAS (although they are usually accepted)
// ----------------------------------------------------------------------------

int dump_histo_tree(Folder *folder)
{
   while( 1 ){
      //printf("\"%s\"\n", (folder->name == NULL) ? "NULL" : folder->name); // warning: comparison of array 'folder->name' equal to a null pointer is always false
      printf("\"%s\"\n", folder->name);
      if( folder->next_subfolder != NULL ){
 	 printf(" / \n");
         dump_histo_tree(folder->next_subfolder);
      }
      if( folder->next_folder != NULL ){
	 folder = folder->next_folder;
      } else {
         break;
      }
   }
   printf(" .. \n");
   return(0);
}
