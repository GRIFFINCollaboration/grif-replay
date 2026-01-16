OBJECTS = config.o grif-replay.o midas-format.o grif-format.o histogram.o \
          web_server.o reorder.o user_sort.o default_sort.o test_config.o

CFLAGS  = -g -O3 -fPIC

grif-replay: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) -rdynamic -lz -ldl -lm -lpthread

midas: midas_module.so

midas_module.so: midas_module.c libmidas.a
	$(CC) $(CFLAGS) -rdynamic -shared -o $@ $^ -lrt -lz -lutil -lnsl -lpthread

.c.o:
	$(CC) -c $(CFLAGS) $<

config.o:       config.c config.h grif-replay.h histogram.h
grif-format.o:  grif-format.c grif-replay.h grif-format.h midas-format.h
midas-format.o: midas-format.c grif-replay.h midas-format.h
default_sort.o: default_sort.c config.h grif-format.h histogram.h
grif-replay.o:  grif-replay.c config.h grif-format.h midas-format.h
midas_module.o: midas_module.c grif-format.h config.h histogram.h midas-format.h midas.h
user_sort.o:    user_sort.c grif-format.h config.h histogram.h
histogram.o:    histogram.c config.h histogram.h
reorder.o:      reorder.c grif-replay.h midas-format.h
web_server.o:   web_server.c histogram.h

#SOURCES = config.c grif-replay.c midas-format.c grif-format.c histogram.c \
           web_server.c reorder.c user_sort.c default_sort.c test_config.c

clean:
	rm -f *.o grif-replay midas_module.so

# if there is a file called "clean", above will fail without this ...
.PHONY: clean
