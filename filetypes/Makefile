ROOTDIR= ..
UTILITYDIR= $(ROOTDIR)/utility
MATRIXDIR= $(ROOTDIR)/matrices
FILETYPESDIR= $(ROOTDIR)/filetypes 
include $(ROOTDIR)/Makefile.options

OBJFILES=  $(UTILITYDIR)/*.o $(MATRIXDIR)/*.o ./*.o 
CFLAGS= -I $(ROOTDIR) -I $(UTILITYDIR) -I $(FILETYPESDIR) -I $(MATRIXDIR) -I$(EIGENDIR) -I $(MATIO)/include -std=c++11
LIBS= -lstdc++ -lm -L$(MATIO)/lib -lmatio 




SparseMatrixFile.o: SparseMatrixFile.hpp SparseMatrixFile.cpp
	$(CC) -c SparseMatrixFile.cpp $(CFLAGS)

SparseMatrixFileASCII.o: SparseMatrixFileASCII.hpp SparseMatrixFileASCII.cpp SparseMatrixFile.o 
	$(CC) -c SparseMatrixFileASCII.cpp $(CFLAGS)

SparseMatrixFileMatlab.o: SparseMatrixFileMatlab.hpp SparseMatrixFileMatlab.cpp SparseMatrixFile.o 
	$(CC) -c SparseMatrixFileMatlab.cpp $(CFLAGS)

testSparseMatrixFile: testSparseMatrixFile.cpp SparseMatrixFileASCII.o 
	$(CC) -o testSparseMatrixFile testSparseMatrixFile.cpp $(OBJFILES) $(CFLAGS) $(LIBS)



DenseMatrixFile.o: DenseMatrixFile.hpp DenseMatrixFile.cpp
	$(CC) -c DenseMatrixFile.cpp $(CFLAGS)

DenseMatrixFileASCII.o: DenseMatrixFileASCII.hpp DenseMatrixFileASCII.cpp DenseMatrixFile.o 
	$(CC) -c DenseMatrixFileASCII.cpp $(CFLAGS)

DenseMatrixFileMatlab.o: DenseMatrixFileMatlab.hpp DenseMatrixFileMatlab.cpp DenseMatrixFile.o 
	$(CC) -c DenseMatrixFileMatlab.cpp $(CFLAGS)

testDenseMatrixFile: testDenseMatrixFile.cpp DenseMatrixFileASCII.o 
	$(CC) -o testDenseMatrixFile testDenseMatrixFile.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


testMatio: testMatio.cpp
	$(CC) -o testMatio testMatio.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testMatlabFile: testMatlabFile.cpp DenseMatrixFileMatlab.o 
	$(CC) -o testMatlabFile testMatlabFile.cpp $(OBJFILES) $(CFLAGS) $(LIBS)

testSparseMatlabFile: testSparseMatlabFile.cpp SparseMatrixFileMatlab.o 
	$(CC) -o testSparseMatlabFile testSparseMatlabFile.cpp $(OBJFILES) $(CFLAGS) $(LIBS)


sparsefiles: SparseMatrixFile.o SparseMatrixFileASCII.o SparseMatrixFileMatlab.o 

densefiles: DenseMatrixFile.o DenseMatrixFileASCII.o DenseMatrixFileMatlab.o 

objects: sparsefiles densefiles 

tests: testSparseMatrixFile testDenseMatrixFile testMatio testMatlabFile testSparseMatlabFile 

all: objects tests

clean: 
	rm -f *.o testSparseMatrixFile testDenseMatrixFile testMatio testMatlabFile testSparseMatlabFile 

anew: clean all 

