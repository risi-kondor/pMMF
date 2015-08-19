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



#ifndef _Cvector
#define _Cvector

#include "DenseVector.hpp"
#include "Remap.hpp"
#include "GivensRotation.hpp"
#include "KpointRotation.hpp"


class Cvector: public DenseVector, public Serializable{
public: 

  class Virtual;

  Cvector(const Cvector& x): Cvector(x.n){
    #ifdef _MATRIXCOPYWARNING
      cout<<"WARNING: Cvector copied."<<endl;
    #endif
    for(int i=0; i<n; i++) array[i]=x.array[i]; // {std::copy(x.array,x.array+n,array);} 
  }

  Cvector(Cvector&& x): DenseVector(x.n){
    #ifdef _MATRIXCOPYWARNING
    //  cout<<"WARNING: Cvector moved."<<endl;
    #endif
    array=x.array; x.array=nullptr;
  }

  Cvector& operator=(const Cvector& x){
    #ifdef _MATRIXCOPYWARNING
      cout<<"WARNING: Cvector assigned."<<endl;
    #endif
    n=x.n; delete[] array; array=new FIELD[n]; //std::copy(x.array,x.array+n,array); 
    for(int i=0; i<n; i++) array[i]=x.array[i]; return *this;
  }

  Cvector& operator=(Cvector&& x){
    #ifdef _MATRIXCOPYWARNING
    //  cout<<"WARNING: Cvector move-assigned."<<endl;
    #endif
    n=x.n; delete[] array; array=x.array; x.array=nullptr; return *this;
  }

  ~Cvector(){delete[] array;}


public: // constructors

  Cvector(): DenseVector(0){array=NULL;}
  Cvector(const int _n): DenseVector(_n) {array=new FIELD[n];}
  // Cvector(const int _n, const Uninitialized& dummy): DenseVector(_n) {array=new FIELD[n];}
  // Cvector(const int _n, const Zero& dummy): Cvector(_n) {for(int i=0; i<n; i++) array[i]=0;}
  Cvector(const int _n, const Random& dummy); 
  Cvector(const int _n, const FIELD* _array);
  Cvector(const initializer_list<FIELD> list);


public: // named constructors

  static Cvector Random(const int n);
  static Cvector Zero(const int n){Cvector v(n); for(int i=0; i<n; i++) v.array[i]=0; return v;}


public: // converters 

  Cvector(const Cvector& x, const Remap& remap, const bool inverse=false): Cvector(x.n){
    if(!inverse) for(int i=0; i<n; i++) array[i]=x.array[remap.forward[i]];
    else for(int i=0; i<n; i++) array[i]=x.array[remap.backward[i]];}
  Cvector(const Vectorv& x);
  Cvector(const Vectorl& x);
  Cvector(const Vectorh& x);
  Cvector(const EigenVectorXdAdaptor& x);

  template<class TYPE> TYPE convert() const;

  static Cvector merge(const Cvector& x, const Cvector& y);


public: // element access 

  FIELD& operator()(const int i){return array[i];}
  FIELD operator()(const int i) const {return array[i];}

  void foreach(std::function<void(INDEX,FIELD&)> lambda){for(int i=0; i<n; i++) lambda(i,array[i]);}

public: // scalar-valued operations 

  bool operator==(const Cvector& x) {assert(x.n==n);
    for(int i=0; i<n; i++) if (array[i]!=x.array[i]) return false; return true;}

  int nnz() const{int t=0; for(int i=0; i<n; i++) if (array[i]!=0) t++; return t;}

  INDEX argmax() const{ 
    if(n==0) return 0; int best=0; FIELD max=array[best];  
    for(int i=0; i<n; i++) 
      if(array[i]>max) {best=i; max=array[i];}
    return best;
  }

  INDEX argmax_abs() const{
    if(n==0) return 0; int best=0; FIELD max=fabs(array[best]);  
    for(int i=0; i<n; i++) 
      if(fabs(array[i])>max) {best=i; max=fabs(array[i]);}
    return best;
  }

  FIELD norm2() const {FIELD t=0; for(int i=0; i<n; i++) t+=array[i]*array[i]; return t;}
  
  FIELD norm() const { return sqrt(norm2()); }

  FIELD diff2(const Cvector& x) const{
    assert(x.n==n); FIELD t=0; for(int i=0; i<n; i++) t+=(array[i]-x.array[i])*(array[i]-x.array[i]); return t;}

  FIELD dot(const Cvector& x) const{assert(x.n==n);
    FIELD t=0; for(int i=0; i<n; i++) t+=array[i]*x.array[i]; return t;}


public: // in-place operations 

  Cvector& operator*=(const FIELD& x){
    for(int i=0; i<n; i++) array[i]*=x; return *this;}

  Cvector& operator/=(const FIELD& x){
    for(int i=0; i<n; i++) array[i]/=x; return *this;}

  Cvector& operator*=(const Cvector& x){
    assert(x.n==n); for(int i=0; i<n; i++) array[i]*=x.array[i]; return *this;}

  Cvector& operator/=(const Cvector& x){
    assert(x.n==n); for(int i=0; i<n; i++) array[i]/=x.array[i]; return *this;}

  Cvector& add(const Cvector& x){
    assert(n==x.n); for(int i=0; i<n; i++) array[i]+=x.array[i]; return *this;}

  Cvector& add(const Cvector& x, const FIELD c){
    assert(n==x.n); for(int i=0; i<n; i++) array[i]+=c*x.array[i]; return *this;}


public: // vector-valued operations 

  Cvector operator-(const Cvector& x){
    assert(x.n==n); Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]-x.array[i]; return v;}


public: // rotations 

  void apply(const GivensRotation& Q){
    assert(Q.i1<n); assert(Q.i2<n);
    FIELD x1=array[Q.i1]; 
    FIELD x2=array[Q.i2]; 
    array[Q.i1]=Q.cos*x1-Q.sin*x2;
    array[Q.i2]=Q.sin*x1+Q.cos*x2; 
  }

  void apply(const KpointRotation& Q){
    for(int i=0; i<Q.k; i++) assert(Q.ix[i]<n);
    for(int i=0; i<Q.k; i++) {Q.tempp[i]=&array[Q.ix[i]]; Q.temp[i]=*Q.tempp[i];}
    for(int i=0; i<Q.k; i++){
      double s=0; for(int j=0; j<Q.k; j++) s+=Q.q[j*Q.k+i]*Q.temp[j]; *Q.tempp[i]=s;}
  }
  
  void applyInverse(const GivensRotation& Q){
    assert(Q.i1<n); assert(Q.i2<n);
    FIELD x1=array[Q.i1]; 
    FIELD x2=array[Q.i2]; 
    array[Q.i1]=Q.cos*x1+Q.sin*x2;
    array[Q.i2]=-Q.sin*x1+Q.cos*x2; 
  }

  void applyInverse(const KpointRotation& Q){
    for(int i=0; i<Q.k; i++) assert(Q.ix[i]<n);
    for(int i=0; i<Q.k; i++) {Q.tempp[i]=&array[Q.ix[i]]; Q.temp[i]=*Q.tempp[i];}
    for(int i=0; i<Q.k; i++){
      double s=0; for(int j=0; j<Q.k; j++) s+=Q.q[i*Q.k+j]*Q.temp[j]; *Q.tempp[i]=s;}
  }
  

public:

  Cvector(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;


public:  

  FIELD* array;

};



class Cvector::Virtual : public Cvector{
public:
  Virtual(const int _n, FIELD* _array): Cvector(){n=_n; array=_array;}
  Virtual(const Cvector::Virtual& x): Virtual(x.n,x.array){} 
  ~Virtual(){array=NULL;}
};



ostream& operator<<(ostream& stream, const Cvector& v);




#endif
