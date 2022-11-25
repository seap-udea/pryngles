FC=gfortran
FFLAGS=
FLFLAGS=
CC=gcc
CFLAGS=-I. -I/usr/include/python3.10
CLFLAGS=-lm

clean:
	rm -rf *.o *.out *.fo *.fout *~

#############################################################
# FORTRAN COMPILING RULES
#############################################################
%.fo:%.f
	$(FC) $(FFLAGS) -c $^ -o $@ 

%.fout:%.fo pixx.fo 
	$(FC) $(FLFLAGS) $^ -o $@

#############################################################
# C COMPILING RULES
#############################################################
%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@ 

%.out:%.o
	$(CC) $(CLFLAGS) $^ -o $@
