FFLAGS=-g -O2 -fPIC -fallow-argument-mismatch 

TESTC99COMPLEX=test_c99_complex
FC=gfortran
CC=gcc
CXX=g++

MEX=mex -largeArrayDims $(OOFLAG)
MEX_FC=$(FC)
MEX_CC=$(CXX)
MW=mwrap
OBJSUF=o
MODSUF=mod

MWFLAGS=-c99complex -mex     
MEX_LINK=
MEX_FLAGS=

## Mac
MEX_LIBS=-lgfortran -lm
MEX=mex $(OOFLAG)


## Linux
# MEX=mex -largeArrayDims $(OOFLAG)
## no openmp
# MEX_LIBS=-lgfortran -lm
## openmp
# MEX_LIBS=-lgfortran -lgomp -lm



