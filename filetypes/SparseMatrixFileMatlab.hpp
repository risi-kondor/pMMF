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



#ifndef _SparseMatrixFileMatlab
#define _SparseMatrixFileMatlab

#include <matio.h>
#include "SparseMatrixFile.hpp"

class SparseMatrixFile::Matlab: public SparseMatrixFile {
public:

	Matlab(const char* filename) {
		matfile = Mat_Open(filename, MAT_ACC_RDONLY);
		if (matfile == NULL) {
			cout << "Error: file cannot be opened" << endl;
			return;
		}
		matvar = Mat_VarReadNext(matfile);
		if (matvar == NULL) {
			cout << "Error: no variables" << endl;
			return;
		}
		mat_sparse_t *sparse;
		if (matvar->class_type == MAT_C_SPARSE) {
			sparse = (mat_sparse_t*) matvar->data;
		}
		nrows = matvar->dims[0];
		ncols = matvar->dims[1];
		next = reinterpret_cast<double*>(sparse->data);
		Mat_VarPrint(matvar, 1);

		Jc = sparse->jc;
		Ir = sparse->ir;
		njc = sparse->njc;
		ndata = sparse->ndata;
	}

	template<class MATRIX>
	Matlab(const char* filename, const MATRIX& M) {
		cout << M.str() << endl;
		ofstream f;
		f.open(filename);
		for (auto it = M.begin(); it != M.end(); it++)
			f << it.i() + 1 << " " << it.j() + 1 << " " << *it << endl;
		f.close();
		cout
				<< "sparse matrix format created, readable by spconvert() in Matlab"
				<< endl;
	}

	~Matlab() {
		Mat_VarFree(matvar);
		Mat_Close(matfile);
	}

	SparseMatrixFile::iterator begin() {
		return SparseMatrixFile::iterator(*this);
	}

	void get(SparseMatrixFile::iterator& it) {
	}

	SparseMatrixFile& operator>>(IndexValueTriple& dest) {
		int i = indIr;
		int j = indJc;
		int c = 0;
		for (; i < njc - 1; i++) {
			c = 0;
			for (; j < Jc[i + 1] && j < ndata; j++) {
				c++;
				dest.i = Ir[j];
				dest.j = i;
				indJc = j + 1;
				indIr = i;
				break;
			}
			if (c > 0) {
				break;
			}
			indIr = i;
		}
		dest.value = *next++;
		if (!dest.value) {
			dest.i = -1;
			return *this;
		}
		return *this;
	}

public:
	mat_t* matfile;
	matvar_t* matvar;
	double *next;
	int indIr = 0;
	int indJc = 0;
	int* Ir;
	int* Jc;
	int njc;
	int ndata;
};
#endif
