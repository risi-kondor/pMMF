/* 
    pMMF is an open source software library for efficiently computing
    the MMF factorization of large matrices.

    Copyright (C) 2015 Risi Kondor, Nedelina Teneva, Pramod Mudrakarta
    
    This file is part of pMMF.

    pMMF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _MMFbase
#define _MMFbase

#define _UTILCOPYWARNING
#define _MATRIXCOPYWARNING
#define _BLOCKEDCOPYWARNING 
#define _MMFCOPYWARNING
#define _MSSCOPYWARNING

//#define _THREADBANKVERBOSE

#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <initializer_list>
#include <assert.h>
#include <string>
#include <algorithm>
#include <mutex>
#include <atomic>

using namespace std;

typedef int INDEX;
typedef double FIELD;

class Zero{};
class Identity{};
class Uninitialized{};
class Dangling{};
class Nullcols{};
class Dense{};
class Sparse{};
class Symmetric{};
class Transpose{};

class Remap;
class Indexset;
class Rotation;
class ElementaryRotation;
class GivensRotation;
class KpointRotation; 
template<class TYPE> class Clusters;

class Vector;
class Cvector;
class Vectorv;
class Vectorl;
class Vectorh;

class Matrix;
class Cmatrix; 
template<class COLUMNTYPE> class MatrixX;

template<class BLOCK> class BlockedVector;
template<class BLOCK> class Tower;
template<class BLOCK> class Street;
template<class BLOCK> class BlockedMatrix;

class EigenVectorXdAdaptor;
class EigenMatrixXdAdaptor;

class SparseMatrixFile; 
class DenseMatrixFile;

class Bofstream;
class Bifstream;

template<class BLOCK> class MMFprocess;
template<class BLOCK> class MMFmatrix;
class MMFchannel;
class MMFstage;
class MMF;



// Helper classes ---------------------------------------------------------------------------------


class Random{
public: 
  Random(){} 
  Random(const double _p):p(_p){} 
  Random(const int _n):n(_n){}
public:
  double p=0.5;
  int n=0;
};



struct IndexValuePair{
  IndexValuePair():i(0),value(0){}
  string str(){ostringstream result; result<<"("<<i<<","<<value<<")"; return result.str();}
  INDEX i; FIELD value; };

struct IndexValueTriple{
  IndexValueTriple():i(0),j(0),value(0){}
  template<class ITERATOR> 
  IndexValueTriple(const ITERATOR& it):i(it.i()),j(it.j()),value(*it){}
  string str(){ostringstream result; result<<"("<<i<<","<<j<<","<<value<<")"; return result.str();}
  INDEX i; INDEX j; FIELD value; 
};

struct BlockIndexPair{
  BlockIndexPair(): block(-1), index(-1){}
  BlockIndexPair(const INDEX& _block, const INDEX& _index): block(_block), index(_index){}
  string str(){ostringstream result; result<<"("<<block<<","<<index<<")"; return result.str();}
  INDEX block; INDEX index; };

class IndexSet{
public:
  IndexSet(const int _k):k(_k){ix=new INDEX[k];}
  INDEX& operator[](const int i){return ix[i];}
  INDEX operator[](const int i) const {return ix[i];}
  void sort(){std::sort(ix,ix+k);}
  string str() const{ostringstream result; result<<"("; for(int i=0; i<k-1; i++) result<<ix[i]<<","; result<<ix[k-1]<<")"; return result.str();}
  int k;
  int* ix;
};

typedef vector<int> BlockStructure; 

class CoutLock{
public:
  CoutLock(): lock(mx){}
  lock_guard<mutex> lock;
  static mutex mx;
};


// Helper functions -------------------------------------------------------------------------------
template<typename TYPE>
inline void swapp(TYPE& x, TYPE& y) {TYPE t=x; x=y; y=t;}

template<typename TYPE>
inline void replace(TYPE*& ptr, TYPE* result) {delete ptr; ptr=result;}

template<typename TYPE>
inline void move_over(TYPE*& ptr, TYPE*& result) {delete ptr; ptr=result; result=nullptr;}

template<typename TYPE> // this is what std::move is for 
inline TYPE&& pilfer(TYPE& x){return reinterpret_cast<TYPE&&>(x);}

template<typename TYPE>
inline TYPE zombie(TYPE* xptr){
  TYPE x(reinterpret_cast<TYPE&&>(*xptr)); delete xptr; return x;}


// ------------------------------------------------------------------------------------------------
typedef MatrixX<Vectorv> MatrixXv;
typedef MatrixX<Vectorl> MatrixXl;
typedef MatrixX<Vectorh> MatrixXh;

typedef BlockedVector<Cvector> BlockedCvector;
typedef BlockedVector<Vectorv> BlockedVectorv;
typedef BlockedVector<Vectorl> BlockedVectorl;
typedef BlockedVector<Vectorh> BlockedVectorh;

typedef BlockedMatrix<Cmatrix> BlockedCmatrix;
typedef BlockedMatrix<MatrixX<Vectorv> > BlockedMatrixXv;
typedef BlockedMatrix<MatrixX<Vectorl> > BlockedMatrixXl;
typedef BlockedMatrix<MatrixX<Vectorh> > BlockedMatrixXh;

typedef BlockStructure Bstructure;


#include "Rstream.hpp"
#include "Serializable.hpp"
#include "Log.hpp"

#endif
