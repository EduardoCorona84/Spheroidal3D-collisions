FFLAGS=-g -O2 -fPIC -Wargument-mismatch -fallow-argument-mismatch

TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
MEX=   mex -largeArrayDims $(OOFLAG) 
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
OBJSUF=o
MODSUF=mod

MWFLAGS=-c99complex -mex     
MEX_LINK=
MEX_FLAGS=
# ifeq (,$(findstring openmp,$(TARGET)))
# openmp is NOT in TARGET    
# MEX_LIBS=-lgfortran -lm
# else 
# openmp is in TARGET
MEX_LIBS=-lgfortran -lgomp -lm
# endif



