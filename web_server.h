// web_server.h

extern int get_line(int fd, char *buf, int maxlen);
extern int put_line(int fd, char *buf, int length);
int split_cmdurl(char *url);
