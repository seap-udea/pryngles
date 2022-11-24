FC=gfortran
FFLAGS=
FLFLAGS=
CC=gcc
CFLAGS=
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

%.out:%.o cpixx.o 
	$(CC) $(CLFLAGS) $^ -o $@

