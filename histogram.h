#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "config.h"

#define FLOAT_1D             1
#define DOUBLE_1D            2
#define INT32_1D             3
#define INT16_1D             4
#define INT64_1D             5
#define FLOAT_2D            11
#define DOUBLE_2D           12
#define INT32_2D            13
#define INT32_2D_SYMM       14
#define INT64_2D            15

#define INT_1D        INT32_1D
#define INT_2D        INT32_2D
#define INT_2D_SYMM   INT32_2D_SYMM

// to allow flexible histogram ranges and scaling ...
// should provide BOTH variable scaling var[gain+offset]->bin#
//                     and axis scaling bin[gain+offset]->value
//                       (possibly even log option for both)
// -------------
// Do not need bin contents to be floating point (float ~24bits int precision)
// as can also add bin content gain/offset
//    => can keep all histogram bin types as int32
//  **UNLESS need histograms containing > 4.2billion counts in a single bin **
//  in this case, need to change data to void *, and check type when accessing
// -------------
//  3d histograms are easily doable [SPARINGLY] - add 3d stuff later
//     3d histos probably require zero-suppression (already have this code)
// [2d probably now do not need this - unlike 30 years ago]
//  if this was wanted - add controls for disk flushing - op count and/or time
//     i.e. every few thousand ops or every few seconds - update disk
//          possibly better to use mmap - so done automatically

extern int open_folder(Config *cfg, char* path);
extern int close_folder(Config *cfg);
extern int Zero_Histograms(Config *cfg);
extern int delete_histograms(Config *cfg);
extern TH1I *hist_querytitle(Config *cfg, char *name);
extern TH1I *hist_queryhandle(Config *cfg, char *name);
extern int write_histofile(Config *cfg, FILE *fp);
extern Config *read_histofile(char *filename, int config_only);
extern char *next_histotitle(Config *cfg, int reset);
extern char *next_histotree_item(Config *cfg, int rst, int *type, int *asc);
extern int delete_histo_tree(Config *cfg);
extern Config *add_histoset();
extern TH1I *H1_BOOK(Config *cfg, char *name, char *title, int xbins, int xmin, int xmax);
extern int TH1I_Reset(TH1I *);
extern int TH1I_Fill(TH1I *, int bin, int count);
extern int TH1I_SetBinContent(TH1I *, int bin, int value);
extern int TH1I_GetBinContent(TH1I *, int bin);
extern int TH1I_SetValidLen(TH1I *, int bins);
extern TH2I *H2_BOOK(Config *cfg, char *name, char *title, int xbins, int xmin, int xmax, int ybins, int ymin, int ymax);
extern int TH2I_Reset(TH2I *);
extern int TH2I_Fill(TH2I *, int xbin, int ybin, int count);
extern int TH2I_SetBinContent(TH2I *, int xbin, int ybin, int value);
extern int TH2I_GetBinContent(TH2I *, int xbin, int ybin);
extern int TH2I_SetValidLen(TH2I *, int bins);
extern int write_th1I(FILE *fp, void *ptr);

#endif
