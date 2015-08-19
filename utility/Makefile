ROOTDIR= ..
UTILITYDIR= $(ROOTDIR)/utility
MATRIXDIR= $(ROOTDIR)/matrices
include $(ROOTDIR)/Makefile.options

OBJFILES=  ./*.o 
CFLAGS= -I $(ROOTDIR) -I $(UTILITYDIR) -I $(MATRIXDIR) -std=c++11
LIBS= -lstdc++ -lm


#MMFglobal.o: MMFglobal.cpp
#	$(CC) -c MMFglobal.cpp $(CFLAGS) 

Log.o: Log.hpp Log.cpp
	$(CC) -c Log.cpp $(CFLAGS)

testLog: testLog.cpp Log.o 
	$(CC) -o testLog testLog.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


Rstream.o: Rstream.hpp Rstream.cpp
	$(CC) -c Rstream.cpp $(CFLAGS)

testRstream: testRstream.cpp Rstream.o 
	$(CC) -o testRstream testRstream.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


Bofstream.o: Bofstream.hpp Bofstream.cpp
	$(CC) -c Bofstream.cpp $(CFLAGS)

Bifstream.o: Bifstream.hpp Bifstream.cpp
	$(CC) -c Bifstream.cpp $(CFLAGS)

Serializable.o: Serializable.hpp Serializable.cpp
	$(CC) -c Serializable.cpp $(CFLAGS)

testBofstream: testBofstream.cpp Bofstream.o 
	$(CC) -o testBofstream testBofstream.cpp $(OBJFILES) $(CFLAGS) $(LIBS)



TopkList.o: TopkList.hpp TopkList.cpp 
	$(CC) -c TopkList.cpp $(CFLAGS)

GivensRotation.o: GivensRotation.hpp GivensRotation.cpp ElementaryRotation.hpp Rotation.hpp 
	$(CC) -c GivensRotation.cpp $(CFLAGS)

KpointRotation.o: KpointRotation.hpp KpointRotation.cpp ElementaryRotation.hpp Rotation.hpp 
	$(CC) -c KpointRotation.cpp $(CFLAGS)

testRotation: testRotation.cpp GivensRotation.o KpointRotation.o 
	$(CC) -o testRotation testRotation.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


#WorkQueue.o: WorkQueue.hpp WorkQueue.cpp 
#	$(CC) -c WorkQueue.cpp $(CFLAGS)

#testWorkQueue: testWorkQueue.cpp WorkQueue.o 
#	$(CC) -o testWorkQueue testWorkQueue.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


ThreadBank.o: ThreadBank.hpp ThreadBank.cpp 
	$(CC) -c ThreadBank.cpp $(CFLAGS)

ThreadManager.o: ThreadManager.hpp ThreadManager.cpp ThreadBank.hpp 
	$(CC) -c ThreadManager.cpp $(CFLAGS)

testThreadBank: testThreadBank.cpp ThreadBank.o 
	$(CC) -o testThreadBank testThreadBank.cpp $(OBJFILES) $(CFLAGS) $(LIBS)



MultiLoop.o: MultiLoop.hpp MultiLoop.cpp 
	$(CC) -c MultiLoop.cpp $(CFLAGS)

Clusters.o: Clusters.hpp Clusters.cpp 
	$(CC) -c Clusters.cpp $(CFLAGS)





objects: Log.o Rstream.o Bofstream.o Bifstream.o Serializable.o TopkList.o \
GivensRotation.o KpointRotation.o ThreadBank.o ThreadManager.o MultiLoop.o Clusters.o  

tests: testLog testRstream testBofstream testRotation testThreadBank 

all: objects tests

clean: 
	rm -f *.o testLog testRstream testBofstream testRotation testThreadBank 

anew: clean all 
