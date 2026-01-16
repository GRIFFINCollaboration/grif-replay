GRIF_OBJECTS = config.o grif-replay.o midas-format.o grif-format.o histogram.o \
          web_server.o grif-reorder.o user_sort.o default_sort.o test_config.o  odb.o

DRAG_OBJECTS = config.o dragon-replay.o midas-format.o histogram.o web_server.o \
          dragon-reorder.o user_sort.o dragon_sort.o test_config.o \
          dragon-format.o odb.o

CFLAGS  = -g -O3 -fPIC

grif-replay:   SYS = -DGRIFFIN_SORT
dragon-replay: SYS = -DDRAGON_SORT

griffin-replay: grif-replay

dragon-replay: $(DRAG_OBJECTS)
	$(CC) $(CFLAGS) $(SYS) -o $@ $^ -rdynamic -lz -ldl -lm -lpthread

grif-replay: $(GRIF_OBJECTS)
	$(CC) $(CFLAGS) $(SYS) -o $@ $^ -rdynamic -lz -ldl -lm -lpthread

midas: midas_module.so

midas_module.so: midas_module.c libmidas.a
	$(CC) $(CFLAGS) $(SYS) -rdynamic -shared -o $@ $^ -lrt -lz -lutil -lnsl -lpthread

.c.o:
	$(CC) -c $(CFLAGS) $(SYS) $<

odb.o:            odb.c odb.h
config.o:         config.c config.h grif-replay.h histogram.h
grif-format.o:    grif-format.c grif-replay.h grif-format.h midas-format.h
midas-format.o:   .FORCE midas-format.c grif-replay.h midas-format.h dragon-format.h
default_sort.o:   default_sort.c config.h grif-format.h histogram.h
dragon_sort.o:    default_sort.c config.h dragon-format.h histogram.h
grif-replay.o:    grif-replay.c config.h grif-format.h midas-format.h
dragon-replay.o:  grif-replay.c config.h dragon-format.h midas-format.h
dragon-format.o:  dragon-format.c config.h dragon-format.h midas-format.h
midas_module.o:   midas_module.c dragon-format.h config.h histogram.h midas-format.h midas.h
user_sort.o:      user_sort.c dragon-format.h config.h histogram.h
histogram.o:      histogram.c config.h histogram.h
grif-reorder.o:   grif-reorder.c grif-replay.h midas-format.h
dragon-reorder.o: dragon-reorder.c grif-replay.h midas-format.h dragon-format.h
web_server.o:     web_server.c histogram.h

#SOURCES = config.c grif-replay.c midas-format.c grif-format.c histogram.c \
           web_server.c reorder.c user_sort.c default_sort.c test_config.c

clean:
	rm -f *.o grif-replay midas_module.so

.FORCE:

# if there is a file called "clean", above will fail without this ...
.PHONY: clean .FORCE
