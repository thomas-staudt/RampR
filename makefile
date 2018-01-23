CC = gcc
CFLAGS = -Wall -Wfatal-errors -Iinclude -O3 -march=native -std=c99
CFLAGS_SHARED = -shared -fPIC
CFILES = src/rampr.c src/api.c src/routine.c src/step.c src/evolution.c
HFILES = include/rampr.h


all: lib

lib: $(HFILES) $(CFILES)	
	mkdir -p lib
	$(CC) $(CFLAGS) $(CFLAGS_SHARED) $(CFILES) -o lib/librampr.so -lm

clean:
	rm -f pyrampr/*.pyc 
	rm -rf lib
