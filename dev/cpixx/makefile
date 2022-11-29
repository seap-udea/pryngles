#PYTHONDIR=/usr/include/python3.10
PYTHONDIR=/Users/jorgezuluagacallejas/opt/anaconda3/include/python3.9/
FC=gfortran
FFLAGS=
FLFLAGS=
CC=gcc
CFLAGS=-I. -I$(PYTHONDIR)
CLFLAGS=-lm
CLIBFLAGS=-fPIC -shared $(CFLAGS) 

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

%.so:%.c
	$(CC) $(CLIBFLAGS) -o $@ $^
