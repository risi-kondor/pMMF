ROOTDIR= ..
UTILITYDIR= $(ROOTDIR)/utility
MATRIXDIR= $(ROOTDIR)/matrices
FILETYPESDIR= $(ROOTDIR)/filetypes
BLOCKEDDIR= $(ROOTDIR)/blocked 
include $(ROOTDIR)/Makefile.options

OBJFILES= $(UTILITYDIR)/*.o $(MATRIXDIR)/*.o $(FILETYPESDIR)/*.o ./*.o  
CFLAGS= -I $(ROOTDIR) -I $(UTILITYDIR) -I $(FILETYPESDIR) -I $(MATRIXDIR) -I $(EIGENDIR) -std=c++11 
LIBS= -lstdc++ -lm


BlockedRemap.o: BlockedRemap.hpp BlockedRemap.cpp   
	$(CC) -c BlockedRemap.cpp $(CFLAGS)

BlockedVector.o: BlockedVector.hpp BlockedVector.cpp 
	$(CC) -c BlockedVector.cpp $(CFLAGS)

Tower.o: Tower.hpp Tower.cpp 
	$(CC) -c Tower.cpp $(CFLAGS)

Street.o: Street.hpp Street.cpp 
	$(CC) -c Street.cpp $(CFLAGS)

BlockedMatrix.o: BlockedMatrix.hpp BlockedMatrix.cpp Tower.o Street.o  
	$(CC) -c BlockedMatrix.cpp $(CFLAGS)


testTower: testTower.cpp Tower.o 
	$(CC) -o testTower testTower.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testBlockedMatrix: testBlockedMatrix.cpp BlockedMatrix.o  
	$(CC) -o testBlockedMatrix testBlockedMatrix.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testReblocking: testReblocking.cpp BlockedMatrix.o  
	$(CC) -o testReblocking testReblocking.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testBlockedRemap: testBlockedRemap.cpp BlockedRemap.o  
	$(CC) -o testBlockedRemap testBlockedRemap.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


testBofstream: testBofstream.cpp BlockedMatrix.o 
	$(CC) -o testBofstream testBofstream.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


objects: BlockedRemap.o BlockedVector.o Tower.o Street.o BlockedMatrix.o  

tests: testTower testBlockedMatrix testReblocking testBlockedRemap testBofstream 


all: objects tests

clean: 
	rm -f *.o testTower testBlockedMatrix testReblocking testBlockedRemap testRstream testBofstream 


anew: clean all 


#GenericBlockedVector.o: GenericBlockedVector.hpp GenericBlockedVector.cpp 
#	$(CC) -c GenericBlockedVector.cpp $(CFLAGS)
#GenericTower.o: GenericTower.hpp GenericTower.cpp  
#	$(CC) -c GenericTower.cpp $(CFLAGS)



