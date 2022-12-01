#PYTHONDIR=/usr/include/python3.10
PYTHON=/usr/local/bin/python3.11
PYTHONDIR=/Library/Frameworks/Python.framework/Versions/3.11/include/python3.11/
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

#############################################################
# TEST
#############################################################
test_dev:
	make cpixx_dev.so
#	$(PYTHON) test_cpixx_dev.py 
	$(PYTHON) test_cpixx.py 

test_cpixx:
	make cpixx.so
	$(PYTHON) test_cpixx.py 
