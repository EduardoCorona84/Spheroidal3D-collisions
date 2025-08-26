all: linux clean

mac: macos-applesilicon

linux: linux-openmp clean 

linux:
	$(MAKE) -f makefile.mwrap TARGET=linux clean
	$(MAKE) -f makefile.mwrap TARGET=linux -j

linux-openmp:
	$(MAKE) -f makefile.mwrap TARGET=linux-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=linux-openmp -j

macos-applesilicon:
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64 clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64 -j

clean:
	rm -f *.o