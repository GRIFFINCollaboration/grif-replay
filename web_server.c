// server loop-> accept->handle->get_req->get_line(buf,urllen)
//   ->parse_line    first line: RmvTrailSpc, SkipGet/Head, SkipLeadingSpace
//                    copy buf->url up to first space, then skip spaces
//                   IF next 5 chars are HTTP/ return more_flws else return 0
//         SUBSEQUENT LINES - Ignore contents, return more_flws if non-empty
/////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <sys/select.h> /* select(), FD_SET() etc */
#include <netinet/in.h> /* sockaddr_in, IPPROTO_TCP, INADDR_ANY */
#include <signal.h>     /* signal() */
#include <stdlib.h>     /* exit() */
#include <string.h>     /* memset() */
#include <unistd.h>     // write(), read(), close()
#include <ctype.h>      // isxdigit(), isspace()
#include <signal.h>
#include "histogram.h"
#include "web_server.h"

#define MAX_QUEUE_LEN    6
#define REQUEST_TIMEOUT 30 // 30 seconds (was 10 seconds)
#define URLLEN        2048 //      maximum url and other string lengths
#define WEBPORT       9093 //    http standard recommends ~8000bytes BUT
                           // windows browsers will not handle more than 2000

// The following is required on MacOS
 #ifndef SOCK_NONBLOCK
 #include <fcntl.h>
 #define SOCK_NONBLOCK 0
 #endif

#if defined(__APPLE__) &&  defined(__MACH__)
// Operating system is Mac OS
#define PLATFORM_IS_MACOS 1
#else
#define PLATFORM_IS_MACOS 0
#endif

int handle_connection(int fd);

void web_main(int *arg)
{
   struct sockaddr_in sock_addr;
   int sock_fd, client_fd;
   int sockopt=1; // Re-use the socket

   signal(SIGPIPE, SIG_IGN);

    if(PLATFORM_IS_MACOS){ // MAC Os
       if( (sock_fd=socket(AF_INET,SOCK_STREAM|SOCK_NONBLOCK,0)) == -1){
          perror("create socket failed");  exit(-1);
       }
    }else{ // linux etc
      if( (sock_fd=socket(PF_INET,SOCK_STREAM|SOCK_NONBLOCK,IPPROTO_TCP)) == -1){
            perror("create socket failed");  exit(-1);
         }
    }

   setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, &sockopt, sizeof(int));
   memset(&sock_addr, 0, sizeof(sock_addr));
   sock_addr.sin_family = AF_INET;
   sock_addr.sin_port = htons(WEBPORT);
   sock_addr.sin_addr.s_addr = INADDR_ANY;

   if( bind(sock_fd,(struct sockaddr *)&sock_addr, sizeof(sock_addr)) == -1){
      perror("bind failed"); close(sock_fd); return;
   }
   if( listen(sock_fd, MAX_QUEUE_LEN) == -1 ){
      perror("listen failed"); close(sock_fd); return;
   }
   fprintf(stdout,"Launch data server ...\n");
   if( init_config() ){ close(sock_fd); return; }
   while(1){
      extern volatile int shutdown_server;
      // use select and non-blocking accept
      // on select timeouts (~100ms) check shutdown flag
      int num_fd;  fd_set read_fds;  struct timeval tv;
      tv.tv_sec = 0; tv.tv_usec = 100000; // 100ms
      num_fd = 1; FD_ZERO(&read_fds); FD_SET(sock_fd, &read_fds);
      if( select(sock_fd+1, &read_fds, NULL, NULL, &tv) > 0 ){
	 //fprintf(stdout,"calling accept ...\n");
         if( (client_fd=accept(sock_fd, NULL, NULL)) < 0 ){
            perror("accept failed"); // close(sock_fd);  exit(-1);
         } else {
            handle_connection(client_fd); close(client_fd);
	 }
      }
      if( shutdown_server != 0 ){  break; }
   }
   fprintf(stdout,"shutting down data server ...\n");
   close(sock_fd);
   return;
}

///////////////////////////////////////////////////////////////////////////
//////////////         server request handling          ///////////////////
///////////////////////////////////////////////////////////////////////////

#define MORE_FOLLOWS 1
#define OPTIONS_REQUEST 2
#define MAXURLARGS 128
#define ROOTDIR "/home/griffin/daq/analyser"
static char buf[URLLEN], url[URLLEN], filename[URLLEN];
static char url_args[MAXURLARGS][STRING_LEN];
char *histo_list[MAXURLARGS];

int get_request(int fd);
int parse_url(int fd, int *cmd, int *narg);
int send_file(char *filename, int fd);

int handle_connection(int fd)
{
   extern int handle_command(int fd, int narg, char url_args[][STRING_LEN]);
   int request_count, content_type, command, arg;
   char *arg2 = url_args[1];
   char *ptr = url;
   int narg;

   if( (request_count=get_request(fd)) /* #lines */   < 0 ){return(-1);}
   while( *ptr != '\0' ){ if(*(ptr++) == '?'){break;} } // skip any "GET /?"
   if( strncmp(ptr, "cmd=", 4)     != 0 ){ return(-1); }
   if( (narg=split_cmdurl(ptr+4))  <= 0 ){ return(-1); }
   //if( send_header(fd, APP_JSON)    < 0 ){ return(-1); }
  // if( send_http_error_response(fd, 1,(char*)"My custom message here")    < 0 ){ return(-1); } // Test this here but will move
   fprintf(stderr,"URL:%s\n", ptr);
   return( handle_command(fd, narg, url_args) );
}

int split_cmdurl(char *url)  // split ?cmd=XXX&arg1=XXX?...
{
   int i=0, j=0, err=0, name_val=0;
   char *ptr = url;
   if( strlen(url) == 0 ){ return(-1); }
   while( *ptr != '\0' && ptr-url < URLLEN ){
      if( strncmp(ptr,"%20",3) == 0 ){
         url_args[i][j++] = ' '; ptr+=3; continue;
      }
      url_args[i][j++] = *ptr++;
      if( j >= STRING_LEN ){ --j; err=1; }
      if( *ptr != '&' && *ptr != '=' ){ continue; }
      if( err ){ fprintf(stderr,"arg %d too long in %s\n", i, url); }
      url_args[i][j] = 0;
      if( *ptr == '&' ){
         i += ( name_val == 0 ) ? 2 : 1;
         j=0; err=0; name_val=0; ++ptr;
      } else {
         if(  name_val != 0 ){
            fprintf(stderr,"name: %s in %s >1 value\n", url_args[i-1], url);
            --i; // ignore all but last
         }
         ++i; j=0; err=0; name_val=1; ++ptr;
      }
      if( i >= MAXURLARGS ){
         i = MAXURLARGS-1;  fprintf(stderr,"too many args in %s\n", url);
      }
   }
   url_args[i][j] = 0;
   for(j=i+1; j<MAXURLARGS; j++){ url_args[j][i] = 0; } // clear other args
   return(i+1);
}

///////////////////////////////////////////////////////////////////////////
////////////////////      http-specific stuff      ////////////////////////
///////////////////////////////////////////////////////////////////////////
int parse_line(char *buf, int first);
int remove_trailing_space(char *buf);
void decodeurl(char *dst, const char *src); // currently unused

#define TESTA "<HTML>\n<HEAD>\n<TITLE>Served File</TITLE>\n</HEAD>\n\n<BODY>\n"
#define TESTB "<H1>Heading</H1>\nHardcoded test page.</BODY>\n</HTML>\n"
// put_line(fd,TESTA,63-6); put_line(fd,TESTB,56-3);

#define HDRA "HTTP/1.0 200 OK\r\nAccess-Control-Allow-Origin: *\r\n"
#define HDRB "Server: javaspec v0\r\nContent-Type: "

static char content_types[][32]={
   "text/html","text/css","text/javascript","application/json"
};
int send_header(int fd, int type)
{
   put_line(fd, HDRA, 49); put_line(fd, HDRB, 35);
   put_line(fd, content_types[type], strlen(content_types[type]) );
   put_line(fd, (char*)"\r\n\r\n", 4);
   return 0;
}

int send_http_error_response(int fd, int type, char *error_message)
{
   char temp[512], tmp[32];
   sprintf(temp,"\"%s\"\r\n", error_content_bodys[type]);
   put_line(fd, (char*)"HTTP/1.0 ", 9); put_line(fd, error_content_titles[type], strlen(error_content_titles[type]) );
   put_line(fd, (char*)"\r\nAccess-Control-Allow-Origin: *\r\n", 34);
   put_line(fd, (char*)"Content-Type: ", 13); put_line(fd, (char*)"text/html;", 10 ); // Error responses are always just text
   sprintf(tmp, "\r\nContent-Length: %lu", strlen(error_message) );
   put_line(fd, tmp, strlen(tmp) );
   put_line(fd, (char*)"\r\n\r\n", 4); // Empty line to separate header and body
   put_line(fd, error_message, strlen(error_message) );
   put_line(fd, (char*)"\r\n", 2);
   return(-1); // Always quit after an error response
}

int send_file(char *filename, int fd)
{
   char temp[256];
   FILE *fp;
   fprintf(stdout,"sending file: %s\n", filename);
   if( (fp=fopen(filename,"r")) == NULL){
      perror("can't open file"); return(-1);
   }
   while( fgets(temp, 256, fp) != NULL ){
      put_line(fd, temp, strlen(temp) );
   }
   fclose(fp);
   return 0;
}

/* HTTP/1.0: GET HEAD POST - GET and HEAD must be implemented */
/* HTTP/1.0: OPTIONS PUT DELETE TRACE CONNECT */
/* GET URL[ HTTP/Version] (no spaces in URL) */
/* GET - just return data, others return header then data  */
int parse_line(char *buf, int first)
{
   char *p=url;
   //fprintf(stdout,"Rcv:%s", buf);
   remove_trailing_space(buf);
   if(strlen(buf) == 0 ){ return(0); }
   if( first ){
      if(        ! strncmp(buf, "GET ",  4) ){ buf += 4;
      } else if( ! strncmp(buf, "HEAD ", 5) ){ buf += 5;
    //  } else if( ! strncmp(buf, "OPTIONS ", 8) ){ return(OPTIONS_REQUEST);
    } else if( ! strncmp(buf, "OPTIONS ", 8) ){ buf += 8; return(OPTIONS_REQUEST);
      } else {
         fprintf(stdout,"Unimplemented"); return(-1);
      }
      while( *buf == ' ' ){ ++buf; } /* skip space */
      while( *buf != ' ' && *buf != '\0' ){ *p++ = *buf++; } /* copy url */
      while( *buf == ' ' ){ ++buf; } /* skip space */
      *p = '\0';
      if( ! strncmp(buf, "HTTP/",  5) ){
         return(MORE_FOLLOWS);
      }
      return(0);
   }
   return(MORE_FOLLOWS);
}

int remove_trailing_space(char *buf)
{
   char *p = buf + strlen(buf)-1;
   while( p >= buf ){ // if( !isspace(*p) ){ break; }
      if( *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r' ){ break; }
      *(p--) = '\0';
   }
   return(0);
}

// decode ascii hex strings ['0'=48 'A'=65 'a'=97]
//  - convert to uppercase, then subtract either (65-10):letter or 48:digit
void decodeurl(char *dst, const char *src)
{
   char a, b;
   while(*src){ // != '\0'
      if( (*src=='%') && ((a = src[1]) && (b = src[2])) && // a,b != '\0'
                         ( isxdigit(a) && isxdigit(b) ) ){ // a,b xdigits
	 if( a >= 'a' ){ a -=  'A' - 'a'; } //->uppercase  //     [0-9a-fA-F]
         if( a >= 'A' ){ a -= ('A' - 10); } else { a -= '0'; }
         if( b >= 'a' ){ b -=  'A' - 'a'; } //->uppercase
         if( b >= 'A' ){ b -= ('A' - 10); } else { b -= '0'; }
         *dst++ = 16*a+b;  src+=3;
      } else {
         *dst++ = *src++;
      }
   }
   *dst++ = '\0';
}

///////////////////////////////////////////////////////////////////////////
////////////////////           socket I/O          ////////////////////////
///////////////////////////////////////////////////////////////////////////

/* get next request - should be within timeout of connecting             */
/* returns number of lines contained in request                          */
/* Example request:                                                      */
/*    [GET /?cmd=callspechandler&spectrum0=HITPATTERN_Energy HTTP/1.1\n] */
/*    [Host: grifstore1.triumf.ca:9093\n]                                */
/*    [User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:52.0)              */
/*                                        Gecko/20100101 Firefox/52.0\n] */
/*    [Accept: x/x\n]                                                    */
/*    [Accept-Language: en-US,en;q=0.5\n]                                */
/*    [Accept-Encoding: gzip, deflate\n]                                 */
/*    [Referer: http://griffincollaboration.github.io/SpVwr/spVwr.html\n]*/
/*    [Origin: http://griffincollaboration.github.io\n]                  */
/*    [DNT: 1\n]                                                         */
/*    [Connection: keep-alive\n]                                         */
/*    [\n]                                                               */
/*                                                                       */
int get_request(int fd)
{
   int status, line_count=0;
   struct timeval tv;
   fd_set read_fd;

   while(1){
      tv.tv_sec = REQUEST_TIMEOUT;  tv.tv_usec = 0;
      FD_ZERO(&read_fd);  FD_SET(fd, &read_fd);

      /* select first arg is highest fd in any set + 1 */
      /* args 2-4 are 3 separate sets for read, write, exception */
      if( (status = select (fd+1, &read_fd, NULL, NULL, &tv )) < 0 ){
         return(-1); /* error */
      }
      if( status == 0 ){ return(-1); } /* timeout */
      if( get_line(fd, buf, URLLEN) <= 0 ){
         fprintf(stderr,"empty get_line at line%d\n", line_count);  return(-1);
      }
      if( (status = parse_line(buf, (line_count==0) )) < 0 ){
         fprintf(stderr,"parse error [%s]\n", buf);  return(-1);
      }

      if( status == OPTIONS_REQUEST ){
       // This is a preflight CORS request. We just need to send the success header.
       // Send special header
          put_line(fd, "HTTP/1.0 200 OK\r\nAccess-Control-Allow-Origin: *\r\nAccess-Control-Allow-Methods: GET, OPTIONS\r\n", 93);
          put_line(fd, "Access-Control-Max-Age: 86400\r\n", 31); // Allows the response to be cached for 24 hrs.
          put_line(fd, HDRB, 35);
          put_line(fd, content_types[APP_JSON], strlen(content_types[APP_JSON]) );
          put_line(fd, (char*)"\r\n\r\n", 4);
       return(-1); // Do not allow this request to go any further
      }

      ++line_count;
      if( status != MORE_FOLLOWS ){
        break;
      }
   }
   return(line_count);
}

/* read up to (+including) a newline, 1 character at a time */
int get_line(int fd, char *buf, int maxlen)
{
   int i, status;
   for(i=0; i<(maxlen-1); i++){
      if( (status=read(fd, buf, 1)) < 0 ){     /* error  -  could probably */
         perror("get_line failed"); return(-1);/*  continue on some errors */
      }
      if( status == 0  ){ *(buf)  ='\0'; return(i);   } /* EOF */
      if( *buf == '\n' ){ *(buf+1)='\0'; return(i+1); } /* End of Line */
      ++buf;
   }
   return(0);
}

/* write complete line as many characters at a time as we can */
int put_line(int fd, char *buf, int length)
{
   int sent;
   while( length ){
      if( (sent = write(fd, buf, length)) <= 0 ){
         perror("put_line failed"); return(-1);
      }
      length -= sent; buf += sent;
   }
   return(0);
}
