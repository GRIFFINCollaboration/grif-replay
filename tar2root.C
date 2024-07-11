// .x tar2root.C+("test.tar","test.root")
//
// radware matrices were always 4k and 2bytes[.mat]  [later 4 bytes for .m4b]

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#define FILE_BUFSIZ    1048576 // 64k * 4 * 4 [1M]
#define MAX_SPECTRA      65535
#define MAX_2D            1024

typedef struct tar_file_header_struct { // should be 512 bytes
   char name[100];  char    mode[8];  char      uid[8];  char    gid[8];
   char  size[12];  char  mtime[12];  char    cksum[8];  char   type[1];
   char link[100];
   char  magic[6];  char version[2];  char   owner[32];  char group[32];
   char devmaj[8];  char  devmin[8];  char prefix[155];  char   pad[12];
} Tar_file_head;    // uid/gid:xbins/ybins
                // own/group:xmin_xmax  ymin_ymax

Tar_file_head file_head;
int           file_body[128*FILE_BUFSIZ];

TH1F *spectrum[MAX_SPECTRA];
TH2F   *spec2d[MAX_2D];

int tar2root(char *tarfile, char *rootfile)
{
   int i, j, pad, count1, count2, size, xbins, ybins, len;
   char tmp[256], handle[128], title[128], fname[64];
   FILE *fp, *tmp_fp;

   TFile *outf = new TFile(rootfile, "RECREATE","Root Histograms");
   if( outf == 0 ){
      printf("Cannot create file for histograms: %s\n", rootfile);
      return(-1);
   }
   if( (fp=fopen(tarfile,"r")) == NULL){ // can't open
      fprintf(stderr,"can't open file:%s to read\n", tarfile );
      return(-1);
   }
   if( (tmp_fp=fopen("dbgout","w")) == NULL){ // can't open
      fprintf(stderr,"can't open debug file to write\n");
      return(-1);
   }   
   count1 = count2 = 0;
   while( 1 ){
      if( fread( &file_head, 512, 1, fp) < 1 ){ break; }
      memcpy(tmp, file_head.size, 12); tmp[12]=0;
      printf("size:%s\n", tmp);
      if( sscanf(tmp, "%o", &size) < 1 ){
         fprintf(stderr,"can't read histo size from:%s\n", file_head.size);
      }
      pad = (512 - (size % 512)) % 512; // pad to multiple of 512 bytes
      printf("size+pad:%d\n", size+pad);
      if( fread( &file_body, sizeof(char), size+pad, fp) <  size+pad ){
         fprintf(stderr,"short read on histo:%s[%d]\n", file_head.name, size);
         break;
      }
      if( strcmp(file_head.name, "config_file") == 0 ){ continue; }
      memcpy(tmp, file_head.uid, 8); tmp[8]=0;
      if( sscanf(tmp, "%o", &xbins) < 1 ){
         fprintf(stderr,"can't read histo xbins from:%s\n", file_head.uid);
      }
      memcpy(tmp, file_head.gid, 8); tmp[8]=0;
      if( sscanf(tmp, "%o", &ybins) < 1 ){
         fprintf(stderr,"can't read histo ybins from:%s\n", file_head.gid);
      }
      printf("HISTO %s %s %d %d\n", file_head.name, file_head.link, xbins, ybins);
      if( ybins == 0 ){
         spectrum[count1]= new TH1F(file_head.name,file_head.link,xbins,0,xbins);
         if( size+pad > 0 ){
            for(i=0; i<xbins; i++){
               spectrum[count1]->SetBinContent(i, file_body[i]);
            }
         }
         spectrum[count1++]->Write();
      } else {
         spec2d[count2]= new TH2F(file_head.name,file_head.link,xbins,0,xbins,ybins,0,ybins);
         if( size+pad > 0 ){
            for(j=0; j<ybins; j++){ for(i=0; i<xbins; i++){
                  spec2d[count2]->SetBinContent(i,j,file_body[i+j*xbins]);
            }}
         }
         spec2d[count2++]->Write();
      }
   }
   fclose(fp);
   outf->Close();
   delete outf;
   return(0);
}
