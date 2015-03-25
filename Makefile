OBJECTS= draw.o display.o matrix.o utils.o parser.o
DEBUG= -DDEBUG -g
WARNINGS_QUIET= -Wall -Wno-unused-variable -Wno-unused-function
WARNINGS_ALL= -Wall -Wpadded
CFLAGS=$(WARNINGS_QUIET)
LIBS= -lm 
CC= gcc

all: build clean

cli: build-repl clean

build: $(OBJECTS) main.o
	$(CC) $(DEBUG) -o main $(OBJECTS) main.o $(LIBS)

build-repl: $(OBJECTS) repl.o
	$(CC) $(DEBUG) -o cli $(OBJECTS) repl.o $(LIBS) -lreadline

main.o: main.c main.h
	$(CC) $(DEBUG) $(CFLAGS) -c main.c

repl.o: repl.c repl.h
	$(CC) $(DEBUG) $(CFLAGS) -c repl.c

draw.o: draw.c draw.h
	$(CC) $(DEBUG) $(CFLAGS) -c draw.c

display.o: display.c display.h
	$(CC) $(DEBUG) $(CFLAGS) -c display.c

matrix.o: matrix.c matrix.h
	$(CC) $(DEBUG) $(CFLAGS) -c matrix.c

utils.o: utils.c utils.h
	$(CC) $(DEBUG) $(CFLAGS) -c utils.c

parser.o: parser.c parser.h
	$(CC) $(DEBUG) $(CFLAGS) -c parser.c

clean:
	rm -f *.o
