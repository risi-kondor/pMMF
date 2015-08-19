ROOTDIR= ..
UTILITYDIR= $(ROOTDIR)/utility
MATRIXDIR= $(ROOTDIR)/matrices
FILETYPESDIR= $(ROOTDIR)/filetypes
BLOCKEDDIR= $(ROOTDIR)/blocked
include $(ROOTDIR)/Makefile.options

OBJFILES= $(UTILITYDIR)/*.o $(FILETYPESDIR)/*.o $(MATRIXDIR)/*.o $(BLOCKEDDIR)/*.o ./*.o  
CFLAGS= -I $(ROOTDIR) -I $(UTILITYDIR) -I $(FILETYPESDIR) -I $(MATRIXDIR) -I $(BLOCKEDDIR) -I $(EIGENDIR) -I$(MATIO)/include -std=c++11 
LIBS= -lstdc++ -lm -lpthread -L$(MATIO)/lib -lmatio


Wtransform.o: Wtransform.hpp Wtransform.cpp
	$(CC) -c Wtransform.cpp $(CFLAGS)

MMFremap.o: MMFremap.hpp MMFremap.cpp
	$(CC) -c MMFremap.cpp $(CFLAGS)

MMFprocess.o: MMFprocess.hpp MMFprocess.cpp
	$(CC) -c MMFprocess.cpp $(CFLAGS)

MMFmatrix.o: MMFmatrix.hpp MMFmatrix.cpp MMFprocess.o MMFremap.o 
	$(CC) -c MMFmatrix.cpp $(CFLAGS)

MMFchannel.o: MMFchannel.hpp MMFchannel.cpp
	$(CC) -c MMFchannel.cpp $(CFLAGS)

MMFstage.o: MMFstage.hpp MMFstage.cpp MMFchannel.o MMFremap.o 
	$(CC) -c MMFstage.cpp $(CFLAGS)

MMFparams.o: MMFparams.hpp MMFparams.cpp 
	$(CC) -c MMFparams.cpp $(CFLAGS)

MMF.o: MMF.hpp MMF.cpp MMFmatrix.o MMFstage.o  Wtransform.o MMFparams.o  
	$(CC) -c MMF.cpp $(CFLAGS)



testElimination: testElimination.cpp MMFprocess.o  
	$(CC) -o testElimination testElimination.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

testCompressionMap: testCompressionMap.cpp MMFmatrix.o 
	$(CC) -o testCompressionMap testCompressionMap.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

testCumulativeError: testCumulativeError.cpp MMF.o 
	$(CC) -o testCumulativeError testCumulativeError.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

testMMF: testMMF.cpp MMF.o 
	$(CC) -o testMMF testMMF.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

testMMFprod: testMMFprod.cpp MMF.o 
	$(CC) -o testMMFprod testMMFprod.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

testMMFgraph: testMMFgraph.cpp MMF.o 
	$(CC) -o testMMFgraph testMMFgraph.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 



MMFcompress: MMFcompress.cpp MMF.o 
	$(CC) -o MMFcompress MMFcompress.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

MMFsave: MMFsave.cpp MMF.o 
	$(CC) -o MMFsave MMFsave.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

SpectralNorms: SpectralNorms.cpp MMF.o 
	$(CC) -o SpectralNorms SpectralNorms.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 

timeMMFcompress: timeMMFcompress.cpp MMF.o 
	$(CC) -o timeMMFcompress timeMMFcompress.cpp $(OBJFILES) $(CFLAGS) $(LIBS) 



objects: Wtransform.o MMFremap.o MMFprocess.o MMFmatrix.o  MMFchannel.o MMFstage.o MMFparams.o MMF.o 

tests: testCompressionMap testCumulativeError testMMF testElimination testMMFprod testMMFgraph timeMMFcompress

execs: MMFcompress MMFsave 


all: objects tests execs

clean: 
	rm -f *.o testCompressionMap testCumulativeError testMMF testElimination testMMFprod testMMFgraph MMFcompress MMFsave timeMMFcompress 

anew: clean all 
