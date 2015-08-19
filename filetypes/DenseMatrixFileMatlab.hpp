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



#ifndef _DenseMatrixFileMatlab
#define _DenseMatrixFileMatlab

#include <matio.h>
#include "DenseMatrixFile.hpp"
#include <algorithm>

class DenseMatrixFile::Matlab: public DenseMatrixFile {
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
		nrows = matvar->dims[0];
		ncols = matvar->dims[1];
		next = reinterpret_cast<double*>(matvar->data);

		int count = nrows * ncols;
		for (int x = 0; x < nrows; ++x) {
			int count_adjustment = nrows - x - 1;
			for (int y = 0, step = 1; y < ncols;
					++y, step += count_adjustment) {
				int last = count - (y + x * ncols);
				int first = last - step;
				std::rotate(next + first, next + first + 1, next + last);
			}
		}
	}

	template<class MATRIX>
	Matlab(const char* filename, const MATRIX& M) {
		size_t dims[2];
		dims[0] = M.nrows;
		dims[1] = M.ncols;
		matfile = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT); //ACC_RDWR);
		if (matfile == NULL) {
			cout << "Cannot open file for write" << endl;
			return;
		}
		matvar = Mat_VarCreate("M", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
				M.array, 0);
		Mat_VarWrite(matfile, matvar, MAT_COMPRESSION_NONE);
		cout << "done" << endl;
	}

	~Matlab() {
		Mat_VarFree(matvar);
		Mat_Close(matfile);
	}

	DenseMatrixFile::iterator begin() {
		return DenseMatrixFile::iterator(*this);
	}

	void get(DenseMatrixFile::iterator& it) {
	}

	DenseMatrixFile& operator>>(FIELD& dest) {
		dest = *next++;
		return *this;
	}

public:

	mat_t* matfile;
	matvar_t* matvar;
	double *next;

};
#endif

