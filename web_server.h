// web_server.h

extern int get_line(int fd, char *buf, int maxlen);
extern int put_line(int fd, char *buf, int length);
int split_cmdurl(char *url);


///////////////////////////////////////////////////////////////////////////
//////////////         server request handling          ///////////////////
///////////////////////////////////////////////////////////////////////////

#define TEXT_HTML 0
#define TEXT_CSS  1
#define TEXT_JS   2
#define APP_JSON  3

extern int send_header(int fd, int type);          // "HTTP 200\nSvr\nCntntTypeX\n\n"
extern int send_http_error_response(int fd, int type, char *error_message);



///////////////////////////////////////////////////////////////////////////
//////////////     server error codes and responses     ///////////////////
///////////////////////////////////////////////////////////////////////////

#define ERR_HDRA "HTTP/1.0 500 \r\nAccess-Control-Allow-Origin: *\r\n"
// HTTP STATUS CODES...
// 4xx: Client Error - The request contains bad syntax or cannot be fulfilled
// 5xx: Server Error - The server failed to fulfill an apparently valid request
#define STATUS_CODE_400 0 // 400 Bad Request. The server cannot or will not process the request due to client error
#define STATUS_CODE_404 1 // 404 Not Found. The server cannot find the requested resource.
#define STATUS_CODE_422 2 // 422 Unprocessable Content. Well-formed reuqest but was unable to be followed due to semantic errors.
#define STATUS_CODE_500 3 // 500 Generic Server Error
static char error_content_titles[][32]={
   "400 Bad Request",
   "404 Not Found",
   "422 Unprocessable Content",
   "500 Internal Server Error"
};
static char error_content_bodys[][128]={
   "400 Bad Request",
   "404 Not Found",
   "422 Unprocessable Content",
   "500 Internal Server Error"
};
