all: linux clean

linux: linux-openmp

linux:
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64 clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64 -j

linux-openmp:
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64-openmp -j

clean:
	rm -f *.o