//#define RECORDSIZE    (1<<24)  // use 16Mbyte records (any size will do)
#define RECORDSIZE    (1<<20)  // use 1Mbyte records (any size will do)
#define EVBUFSIZE     (1<<26)  // allow up to 64Mbyte events
#define MAX_BANK_SIZE (1<<18)  // allow up to 256kbyte Banks

typedef struct midas_event_header_struct {
   short     event_id;   short trigger_mask;   int     serial_num;
   int      timestamp;   int      data_size;
} Midas_event_header;

typedef struct midas_allbank_header_struct {
   int    allbanksize;   int          flags;
} Midas_allbank_header;

typedef struct midas_bank16_header_struct { /* old 16 bit version of header */
   char       name[4];   short    data_type;   unsigned short data_size;
} Midas_bank16_header;

typedef struct midas_bank_header_struct {
   char       name[4];   int      data_type;   int      data_size;
} Midas_bank_header;

//#define BANK_BUFSIZE (64*1024*1024)
#define BANK_BUFSIZE (4*1024*1024) // this is 32bit-words
extern unsigned bankbuf[BANK_BUFSIZE];
extern volatile unsigned long bankbuf_wrpos;
extern volatile unsigned long bankbuf_rdpos;
// without brackets, and using wrpos = bankbuf_wrpos % BANK_BUFSIZE;
//    when bankbuf_wrpos = 34104, wrpos = 58720256

extern void midas_main(Sort_status *arg);
extern int next_record(Sort_status *arg);
extern int next_event(Sort_status *arg);
extern int next_bank(Sort_status *arg, char **bank_name);
extern int copy_bank(unsigned *ptr, int size);

#define REORDER_BUFSIZE (1024*1024)
#define EVENTBUFSIZE (8*REORDER_BUFSIZE)
extern unsigned event_buffer[EVENTBUFSIZE+128]; // add space at end for debug
extern volatile unsigned long eventbuf_rdpos;
extern volatile unsigned long eventbuf_wrpos;
