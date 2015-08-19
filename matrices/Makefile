ROOTDIR=..
UTILITYDIR= $(ROOTDIR)/utility
FILETYPESDIR= $(ROOTDIR)/filetypes
MATRIXDIR= $(ROOTDIR)/matrices
include $(ROOTDIR)/Makefile.options

OBJFILES= $(UTILITYDIR)/*.o $(FILETYPESDIR)/*.o ./*.o  
CFLAGS= -I $(ROOTDIR) -I $(UTILITYDIR) -I $(FILETYPESDIR) -I $(MATRIXDIR) -I $(EIGENDIR) -std=c++11 
LIBS= -lstdc++ -lm


Remap.o: Remap.hpp Remap.cpp
	$(CC) -c Remap.cpp $(CFLAGS)


Vector.o: Vector.hpp Vector.cpp
	$(CC) -c Vector.cpp $(CFLAGS)

Cvector.o: Cvector.hpp Cvector.cpp Vector.hpp 
	$(CC) -c Cvector.cpp $(CFLAGS)

Vectorv.o: Vectorv.hpp Vectorv.cpp SparseVector.hpp  
	$(CC) -c Vectorv.cpp $(CFLAGS)

Vectorl.o: Vectorl.hpp Vectorl.cpp SparseVector.hpp 
	$(CC) -c Vectorl.cpp $(CFLAGS)

Vectorh.o: Vectorh.hpp Vectorh.cpp SparseVector.hpp 
	$(CC) -c Vectorh.cpp $(CFLAGS)



Matrix.o: Matrix.hpp Matrix.cpp Matrix.hpp 
	$(CC) -c Matrix.cpp $(CFLAGS)

Cmatrix.o: Cmatrix.hpp Cmatrix.cpp Matrix.o MatrixX.hpp 
	$(CC) -c Cmatrix.cpp $(CFLAGS)



MatrixX.o: MatrixX.hpp MatrixX.cpp Matrix.o Vectorv.o Vectorl.o Vectorh.o  
	$(CC) -c MatrixX.cpp $(CFLAGS)


testRemap: testRemap.cpp Remap.o 
	$(CC) -o testRemap testRemap.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testVector: testVector.cpp Vector.o Cvector.o 
	$(CC) -o testVector testVector.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testMatrix: testMatrix.cpp Cmatrix.o SparseMatrix.o 
	$(CC) -o testMatrix testMatrix.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testCmatrix: testCmatrix.cpp Cmatrix.o    
	$(CC) -o testCmatrix testCmatrix.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testMatrixX: testMatrixX.cpp MatrixX.o 
	$(CC) -o testMatrixX testMatrixX.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

timeMatrixX: timeMatrixX.cpp MatrixX.o 
	$(CC) -o timeMatrixX timeMatrixX.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testRotations: testRotations.cpp Vectorv.o Vectorl.o Vectorh.o 
	$(CC) -o testRotations testRotations.cpp $(OBJFILES) $(CFLAGS) $(LIBS)



helpers: Remap.o 

vectors: Vector.o Cvector.o Vectorv.o Vectorl.o Vectorh.o 

matrices: Matrix.o Cmatrix.o MatrixX.o 

objects: helpers vectors matrices 

tests: testRemap testVector testMatrixX timeMatrixX testRotations testCmatrix 

all: objects tests

clean: 
	rm -f *.o testRemap testVector testMatrixX timeMatrixX testIterator testRotations \
testRstream testBofstream testCmatrix 


anew: clean all 
