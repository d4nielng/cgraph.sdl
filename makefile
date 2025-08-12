CC=g++
CFLAGS=-O2
LFLAGS=
SDL_CONFIG=$(shell sdl2-config --cflags --libs)

#.phony: all clean

ifeq ($(OS),Windows_NT)
	RM=del
	EXE=.exe
	LFLAGS+=-s
else
	UNAME_S := $(shell uname -s)
	EXE=
	RM=rm -f
endif

all: demo sierpinski cantor mandel viewbmp nqueens tree hilbert tetris

.cpp:
	$(CC) $(SDL_CONFIG) $(CFLAGS) $< -o $@ $(LFLAGS)

#.cpp.exe:
#	$(CC) $(CFLAGS) $< -o $@ $(LFLAGS)

clean:
	$(RM) demo$(EXE)
	$(RM) sierpinski$(EXE)
	$(RM) viewbmp$(EXE)
	$(RM) cantor$(EXE)
	$(RM) mandel$(EXE)
	$(RM) nqueens$(EXE)
	$(RM) hilbert$(EXE)
	$(RM) tree$(EXE)
