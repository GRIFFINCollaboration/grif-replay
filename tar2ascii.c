// gcc -g -o tar2ascii tar2ascii.c -lz

#include <stdio.h>
#include <stdlib.h>   // malloc
#include <string.h>   // memset
#include <unistd.h>   // chdir
#include <sys/stat.h>  // mkdir
#include <sys/types.h> // chdir
#include <errno.h>     // errno
#include <zlib.h>      //

#define DEFAULT_OUTDIR "tar2ascii_dir" // do not use program name "tar2ascii"
#define SMALL_HISTO_BINS 65536    // 64k bins
#define FILE_BUFSIZ     (1024*1024*1024) // allows matrices up to 16k*16k*4

typedef struct file_header_struct { // should be 512 bytes
   char name[100];  char    mode[8];  char      uid[8];  char    gid[8];
   char  size[12];  char  mtime[12];  char    cksum[8];  char   type[1];
   char link[100];
   char  magic[6];  char version[2];  char   owner[32];  char group[32];
   char devmaj[8];  char  devmin[8];  char prefix[155];  char   pad[12];
} File_head;    // uid/gid:xbins/ybins
                // own/group:xmin_xmax  ymin_ymax

static File_head file_head;
static char      file_body[FILE_BUFSIZ];
static char      path[2048];
static int      *histo_data;
static FILE     *histo_fp;

int read_histofile(FILE *fp);
int check_path(char *path);
int write_histo(char *path, int xmin, int xmax, int xbins, int ymin, int ymax, int ybins);

#define COMPRESS_BUFSIZ 270000000
// allow for 8k^2 = 65M-chan * 4bytes = 256Mbytes [~270000000]
static char compress_buf[COMPRESS_BUFSIZ];
int decompress_buffer(char *input, int size);

int main(int argc, char *argv[])
{
   struct stat info;
   int len;
   
   if( argc >= 2 ){
      if( (histo_fp=fopen(argv[1],"r")) == NULL){ // can't open
         fprintf(stderr,"can't open file:%s to read\n", argv[1] );
         exit(-1);
      }
      len = strlen(argv[1]);
      if( strncmp(argv[1]+len-4, ".tar", 4) == 0 ){ len -=4; }
      memcpy(path, argv[1], len); path[len] = 0; 
   } else { histo_fp = stdin; sprintf(path, "%s", DEFAULT_OUTDIR); }

   if( !(stat(path, &info) == 0 && S_ISDIR(info.st_mode)) ){
       if( mkdir(path, 0777) != 0 ){
          fprintf(stderr,"can't create directory %s %s\n", path, strerror(errno) );
          exit(-1);
       }
   }
   if( chdir(path) ){
      fprintf(stderr,"can't change directory %s %s\n", path, strerror(errno) );
      exit(-1);
   }
   read_histofile(histo_fp);
   exit(0);
}

int read_histofile(FILE *fp)
{
   unsigned size, xbins, xmin, xmax, ybins, ymin, ymax, len, data_size;
   long file_offset = 0;
   int pad, err=0, bins;
   char tmp[64];

   while( 1 ){
      if( fread( &file_head, sizeof(File_head), 1, fp) < 1 ){ break; }
      memcpy(tmp, file_head.size, 12); tmp[12]=0;
      if( sscanf(tmp, "%o", &size) < 1 ){
         fprintf(stderr,"can't read histo size from:%s\n", file_head.size);
         err=1; break;
      }
      file_offset += sizeof(File_head);
      pad = (512 - (size % 512)) % 512; // pad to multiple of 512 bytes
      if( strcmp(file_head.name, "config_file") == 0 ){
         if( fseek(fp, size+pad, SEEK_CUR) < 0 ){
            fprintf(stderr,"short seek histo:%s[%d]\n",file_head.name,size);
            err=1; break;
         }
         file_offset += size+pad;
         continue;
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
      memcpy(path, file_head.prefix, len); path[len]=0;
      check_path(path);
      path[len]='/';
      sprintf(path+len+1, "%s.txt", file_head.name);
      if( fread( &file_body, sizeof(char), size+pad, fp) <  size+pad ){
         fprintf(stderr,"short read histo:%s[%d]\n",file_head.name,size);
         err=1; break;
      }
      data_size = size;
      bins = ( ybins == 0 ) ? xbins : xbins*ybins;
      if( (histo_data = (int *)malloc(bins*sizeof(int))) == NULL){
         fprintf(stderr,"read_histo_data: data malloc failed\n");
         return(-1);
      }
      memset(histo_data, 0, bins*sizeof(int));
      // compressed data starts with signed bytes: 120,-100
      if( bins > SMALL_HISTO_BINS && file_body[0]==120 && file_body[1]==-100 ){
         data_size = decompress_buffer(file_body, data_size);
         memcpy(histo_data, compress_buf, data_size);
      } else {
         memcpy(histo_data, file_body, data_size);
      }
      write_histo(path, xmin, xmax, xbins, ymin, ymax, ybins);
   }
   if( err ){ return(-1); }
   return(0);
}

int check_path(char *path)
{
   struct stat info;
   char *ptr = path;

   while(1){
      while( *ptr != '/' && *ptr != 0 ){ ++ptr; }
      if( *ptr == '/' ){
         *ptr = 0;
         if( !(stat(path, &info) == 0 && S_ISDIR(info.st_mode)) ){
            if( mkdir(path, 0777) != 0 ){
               fprintf(stderr,"can't create directory %s %s\n", path, strerror(errno) );
               return(-1);
            }
         }
         *ptr = '/'; ++ptr; continue;
      }
      if( !(stat(path, &info) == 0 && S_ISDIR(info.st_mode)) ){
         if( mkdir( path, 0777) != 0 ){
            fprintf(stderr,"can't create directory %s %s\n", path, strerror(errno) );
            return(-1);
         }
      }
      break;
   }
   return(0);
}

int write_histo(char *path, int xmin, int xmax, int xbins, int ymin, int ymax, int ybins)
{
   FILE *fp;
   int i, j, k;
   if( (fp=fopen(path,"w")) == NULL){ // can't open
      fprintf(stderr,"can't open file:%s to write\n", path );
      return(-1);
   }
   if( ybins == 0 ){
      for(i=0; i<xbins; i++){
         if( histo_data[i] == 0 ){ continue; }
         fprintf(fp, "%d %d\n", i, histo_data[i]);
      }
   } else {
      //if(file_head.type[0] == 'C'){  ybins=SYMMETERIZE; }
      for(j=0; j<ybins; j++){
      for(i=0; i<xbins; i++){
         if( (k=histo_data[i+j*xbins]) == 0 ){ continue; }
         fprintf(fp, "%d %d %d\n", i, j, k );
      }
      }
   }
   fclose(fp);
   return(0);
}

///////////////////////////////////////////////////////////////////////////
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
