
#define HEAD_EVENT 1
#define TAIL_EVENT 2

#define  V792_MAXCHAN 32  // number of channels on module
#define V1190_MAXCHAN 64  //                "

typedef struct io32_event_struct { // 6 words (rest are duplicates)
   int fw_version; // IO32 firmware version number
   int trig_count; // since BOR
   int trig_ts; int readout_begin; int readout_end; int trig_bitmask;
} Io32_event;

#define MAX_IO32_FIFOWORDS 4
typedef struct io32_fifo_struct { // 6+4 = 10 words
   int fw_version;  int timestamp;  int ts_route;
   int ts_ctrl;  int rollover; int nwords; int data[MAX_IO32_FIFOWORDS];
} Io32_fifo;

#define MAX_V1190_WORDS 40
typedef struct v1190_event_struct { // 5+16 = 21
    int ev_count;  short bunch_id; short status; 
   int time_tag;  int d_count;     int data[MAX_V1190_WORDS];
} V1190_event;

#define MAX_V792_WORDS 40 // sometimes do get all 32 channels
typedef struct v792_event_struct { // 32+2 = 34
   int ev_count;
   int d_count;     int data[MAX_V792_WORDS];
} V792_event;

// sort into these structures - [keep original data - for ungainmatched etc]
#define BGO_MAXCHAN  32 // currently 30 are used
#define SB_MAXCHAN    2
#define NAI_MAXCHAN   2
#define DSSD_MAXCHAN 32
#define DSSD_TDCCHAN  2  // front, back
#define IC_MAXCHAN    5
#define MCP_MAXCHAN   4
#define MCP_TDCCHAN   2  // tac, tdc0, tdc1

#define SUBSYS_BGO    1 // subsys_name[][] in default_sort,c
#define SUBSYS_SB     2 //   should match the ordering here
#define SUBSYS_NAI    3
#define SUBSYS_DSSD   4
#define SUBSYS_IC     5
#define SUBSYS_MCP    6
#define SUBSYS_MCPTAC 7
#define SUBSYS_GE     8
#define SUBSYS_XTDC   9
#define SUBSYS_RFTDC 10
#define SUBSYS_TDC0  11

typedef struct head_data {
   int  xtdc; int rftdc; int tdc0; 
   float ge_energy;
   float bgo_energy[BGO_MAXCHAN];
   int   bgo_time  [BGO_MAXCHAN];
   float sb_energy [SB_MAXCHAN];
   float nai_energy[NAI_MAXCHAN];
} Head_data;

typedef struct tail_data {
   int  xtdc; int rftdc; int tdc0; 
   float dssd_energy[DSSD_MAXCHAN];
   int   dssd_time  [DSSD_TDCCHAN]; // front, back
   float ic_energy  [IC_MAXCHAN];
   int   ic_time    [IC_MAXCHAN];
   float mcp_energy [MCP_MAXCHAN];
   float mcptac_energy;
   int   mcp_time   [MCP_TDCCHAN];  // tac, tdc0, tdc1
} Tail_data;

// can save half the memory space using union for head/tail data
typedef union head_tail_data_union {
   Head_data head_data;
   Tail_data tail_data;
} Head_tail_data;

#define DRAGON_EVENT_MARK 0xfedccdef // help recover resync if lost for some reason 
typedef struct dragon_event_struct {  // raw:110 words - 10x grifevents size 
   int  begin_marker;                 // full:~190       19x grifevents size 
   int  type;         // head or tail
   long ts;           // trigger timestamp from io32_fifo
   Io32_event   io32; //  6 words
   Io32_fifo io32_ts; // 10 words
   V792_event  v792a; // 34 words
   V792_event  v792b; // 34 words  // v792b is tail only
   V1190_event v1190; // 21 words
   Head_tail_data head_tail_data;
   //Head_data head_data;
   //Tail_data tail_data;
} Dragon_event;

#define DRAGON_EVENTWORDS (sizeof(Dragon_event)/sizeof(int))

#define EVT_BUFSIZE         16384

extern int   ge_adc_chan, ge_adc_module, head_xtdc_chan, head_rftdc_chan, head_tdc0_chan;
extern float ge_adc_slope, ge_adc_offset, head_xtdc_slope, head_xtdc_offset;
extern float head_rftdc_slope, head_rftdc_offset, head_tdc0_slope, head_tdc0_offset;
extern int   bgo_adc_chan    [BGO_MAXCHAN];
extern int   bgo_adc_module  [BGO_MAXCHAN];
extern float bgo_adc_slope   [BGO_MAXCHAN];
extern float bgo_adc_offset  [BGO_MAXCHAN];
extern int   bgo_adc_pedestal[BGO_MAXCHAN];
extern int   bgo_tdc_chan    [BGO_MAXCHAN];
extern float bgo_tdc_slope   [BGO_MAXCHAN];
extern float bgo_tdc_offset  [BGO_MAXCHAN];
extern float bgo_tdc_xposn   [BGO_MAXCHAN];
extern float bgo_tdc_yposn   [BGO_MAXCHAN];
extern float bgo_tdc_zposn   [BGO_MAXCHAN];
// also bgo/hv/channel
extern int   sb_adc_chan     [SB_MAXCHAN];
extern int   sb_adc_module   [SB_MAXCHAN];
extern float sb_adc_slope    [SB_MAXCHAN];
extern float sb_adc_offset   [SB_MAXCHAN];
extern int   nai_adc_chan    [NAI_MAXCHAN];
extern int   nai_adc_module  [NAI_MAXCHAN];
extern float nai_adc_slope   [NAI_MAXCHAN];
extern float nai_adc_offset  [NAI_MAXCHAN];
// tail detectors ...
extern int   tail_xtdc_chan, tail_rftdc_chan, tail_tdc0_chan;
extern float tail_xtdc_slope, tail_xtdc_offset, tail_rftdc_slope, tail_rftdc_offset;
extern float tail_tdc0_slope, tail_tdc0_offset; 
extern int   dssd_adc_chan    [DSSD_MAXCHAN];
extern int   dssd_adc_module  [DSSD_MAXCHAN]; // 32 dssd chan are in two 16-chan adcs
extern float dssd_adc_slope   [DSSD_MAXCHAN];
extern float dssd_adc_offset  [DSSD_MAXCHAN];
extern int   dssd_tdc_chan    [DSSD_TDCCHAN]; // was called front/back
extern float dssd_tdc_slope   [DSSD_TDCCHAN];
extern float dssd_tdc_offset  [DSSD_TDCCHAN];
extern int   ic_adc_chan      [IC_MAXCHAN];
extern int   ic_adc_module    [IC_MAXCHAN];
extern float ic_adc_slope     [IC_MAXCHAN];
extern float ic_adc_offset    [IC_MAXCHAN];
extern int   ic_tdc_chan      [IC_MAXCHAN];
extern float ic_tdc_slope     [IC_MAXCHAN];
extern float ic_tdc_offset    [IC_MAXCHAN];
extern int   mcp_adc_chan     [MCP_MAXCHAN];
extern int   mcp_adc_module   [MCP_MAXCHAN];
extern float mcp_adc_slope    [MCP_MAXCHAN];
extern float mcp_adc_offset   [MCP_MAXCHAN];
extern int   mcptac_adc_chan;
extern int   mcptac_adc_module;
extern float mcptac_adc_slope;
extern float mcptac_adc_offset;
extern int   mcp_tdc_chan     [MCP_TDCCHAN]; // was called tdc0,tdc1
extern float mcp_tdc_slope    [MCP_TDCCHAN];
extern float mcp_tdc_offset   [MCP_TDCCHAN];
//derived tables
extern int head_adc_dettype[2*V792_MAXCHAN]; extern int head_adc_dstchan[2*V792_MAXCHAN]; 
extern int head_tdc_dettype[ V1190_MAXCHAN]; extern int head_tdc_dstchan[ V1190_MAXCHAN]; 
extern int tail_adc_dettype[2*V792_MAXCHAN]; extern int tail_adc_dstchan[2*V792_MAXCHAN]; 
extern int tail_tdc_dettype[ V1190_MAXCHAN]; extern int tail_tdc_dstchan[ V1190_MAXCHAN]; 

int unpack_io32_bank(Dragon_event *evt, int *ptr, int len, int type);
int unpack_io32_fifobank(Dragon_event *evt, int *ptr, int len, int type);
int unpack_v1190_bank(Dragon_event *evt, int *ptr, int len, int type);
int unpack_v792_bank(Dragon_event *evt, int *ptr, int len, int id, int type);
