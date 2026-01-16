
#define HEAD_EVENT 1
#define TAIL_EVENT 2

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

#define MAX_V1190_WORDS 16
typedef struct v1190_event_struct { // 5+16 = 21
    int ev_count;  short bunch_id; short status; 
   int time_tag;  int d_count;     int data[MAX_V1190_WORDS];
} V1190_event;

#define MAX_V792_WORDS 32 // sometimes do get all 32 channels
typedef struct v792_event_struct { // 32+2 = 34
   int ev_count;
   int d_count;     int data[MAX_V792_WORDS];
} V792_event;

//typedef struct bgo_event_struct {
//} BGO_event;
//typedef struct mcp_event_struct {
//} Mcp_event;
//typedef struct mcp_event_struct {
//} Dssd_event;

#define DRAGON_EVENT_MARK 0xfedccdef // help recover resync if lost for some reason 
typedef struct dragon_event_struct {  // 110 words - 10 times as large as grif events 
   int  begin_marker;
   int  type;         // head or tail
   long ts;           // trigger timestamp from io32_fifo
   Io32_event   io32; //  6 words
   Io32_fifo io32_ts; // 10 words
   V792_event  v792a; // 34 words
   V792_event  v792b; // 34 words  // v792b is tail only
   V1190_event v1190; // 21 words
} Dragon_event;

#define DRAGON_EVENTWORDS (sizeof(Dragon_event)/sizeof(int))

#define EVT_BUFSIZE         16384

int unpack_io32_bank(Dragon_event *evt, int *ptr, int len, int type);
int unpack_io32_fifobank(Dragon_event *evt, int *ptr, int len, int type);
int unpack_v1190_bank(Dragon_event *evt, int *ptr, int len, int type);
int unpack_v792_bank(Dragon_event *evt, int *ptr, int len, int id, int type);

