LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include

CC = gcc
CFLAGS = -Wall -O3 -std=c99
AR = ar rc

LIB = libd3q15.a

libobjs = d3q15.o checks.o noise.o \
	init_uniform.o init_point.o \
	bc_pbc.o bc_noslip.o bc_freeslip.o bc_wall.o \
	force_none.o force_const.o

headers = d3q15.h Lattice.h noise.h eigenvectors.h

$(LIB) : $(libobjs)
	$(AR) $@ $^
	ranlib $@

%.o : %.c $(headers)
	$(CC) $(CFLAGS) -c $<

eigenvectors.h : eigenvectors.py
	python $< C > $@

.PHONY : clean
clean :
	-rm libd3q15.a $(libobjs) eigenvectors.h
