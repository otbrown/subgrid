# directory paths
D3Q15DIR = ../libd3q15
VPATH = $(D3Q15DIR)

# compiler
CC = gcc
CFLAGS = -Wall -O3

# dependencies
HEADERS = d3q15.h
LIBS = -ld3q15 -lm

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -I$(D3Q15DIR)

pulse: pulse.o $(LIBS)
	$(CC) $< -o $@ -L$(D3Q15DIR) $(LIBS)

all: pulse

.PHONY: clean

clean:
	rm -f *.o pulse
