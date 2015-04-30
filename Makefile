OBJECTS= draw.o display.o matrix.o utils.o parser-old.o gmath.o symtab.o print_pcode.o stack.o
DEBUG= -DDEBUG -g
WARNINGS_QUIET= -Wall -Wno-unused-variable -Wno-unused-function
WARNINGS_ALL= -Wall -Wpadded
CFLAGS=$(WARNINGS_QUIET)
LIBS= -lm 
CC= gcc

all: build clean

cli: build-repl clean

build: $(OBJECTS) lex.yy.c y.tab.c y.tab.h
	$(CC) $(DEBUG) -o mdl $(OBJECTS) lex.yy.c y.tab.c y.tab.h $(LIBS)

build-parser: $(OBJECTS) main.o
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

parser-old.o: parser-old.c parser-old.h
	$(CC) $(DEBUG) $(CFLAGS) -c parser-old.c

gmath.o: gmath.c gmath.h
	$(CC) $(DEBUG) $(CFLAGS) -c gmath.c

symtab.o: symtab.c symtab.h
	$(CC) $(CFLAGS) -c symtab.c

stack.o: stack.c stack.h
	$(CC) $(CFLAGS) -c stack.c 

print_pcode.o: print_pcode.c print_pcode.h y.tab.h
	$(CC) $(CFLAGS) -c print_pcode.c

lex.yy.c: mdl.l y.tab.h 
	flex -I mdl.l

y.tab.c: mdl.y symtab.h parser.h
	bison -d -y mdl.y

y.tab.h: mdl.y 
	bison -d -y mdl.y

clean:
	rm -f *.o
	rm -f y.tab.c y.tab.h
	rm -f lex.yy.c
