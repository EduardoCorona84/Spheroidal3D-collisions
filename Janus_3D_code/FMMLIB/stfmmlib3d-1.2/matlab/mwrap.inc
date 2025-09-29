## Edited in 2025 by Nic Rummel to run on BLANCA and ALPINE 
# Set compile cmds
MEX=mex -largeArrayDims 
MW=mwrap
MEX_FC=gfortran
MEX_CC=g++
# Set compile flags
FFLAGS=-g -O3 --openmp -fPIC -fallow-argument-mismatch
CFLAGS= -fPIC -O3 -march=native -funroll-loops -std=gnu17
# Idk if this is necessary
TESTC99COMPLEX=test_c99_complex
MWFLAGS=-c99complex -mex	
MEX_LIBS=-lgomp -lgfortran -lm
MEX_LINK=-output $(MEX_GATEWAY).mexa64
MEX_FLAGS= -D_OPENMP
OBJSUF=o
MODSUF=mod
