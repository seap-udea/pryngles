#PYTHONDIR=/usr/include/python3.10
PYTHON=/usr/local/bin/python3.11
PYTHONDIR=/Library/Frameworks/Python.framework/Versions/3.11/include/python3.11/
CC=gcc
CFLAGS=-I. -I$(PYTHONDIR)
CLFLAGS=-lm
CLIBFLAGS=-fPIC -shared $(CFLAGS)

clean:
	rm -rf *.o *.out *.fo *.fout *~

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
test_cpixx:
	make cpixx.so
	$(PYTHON) test_cpixx.py 

test_mul:
	make cpixx.so
	$(PYTHON) test_mularray.py 
