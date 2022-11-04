# Makefile for the Fourier geometries program
#
# Author: D. M. Stam
# Date: July 2013
# 

FC= gfortran
FFLAGS=  -g -O -fbounds-check -Wall

OBJECTS = main.o \
          bracks.o \
          splint.o \
          rdfous.o \
          getgeos.o \
          spline.o 

pixx: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o pixx

main.o: main.f max_incl
rdfous.o: rdfous.f max_incl
getgeos.o: getgeos.f max_incl

clean: 
	rm -f pixx
	rm -f *.[oc]
