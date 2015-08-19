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



#ifndef _SparseMatrixFileASCII
#define _SparseMatrixFileASCII

#include "SparseMatrixFile.hpp"

class SparseMatrixFile::ASCII: public SparseMatrixFile {
public:

	ASCII(const char* filename) {
		ifs.open(filename);
		char buffer[255];
		ifs.get(buffer, 255);
		ifs.close();
		ifs.open(filename);
		int nextracted = 0;
		while (ifs.good() && ifs.tellg() < strlen(buffer)) {
			float b;
			ifs >> b;
			if (!ifs.fail())
				nextracted++;
		}
		if (nextracted == 2) {
			ifs.close();
			ifs.open(filename);
			ifs >> nrows >> ncols;
			return;
		}
		if (nextracted == 3) {
			ifs.close();
			ifs.open(filename);
			nrows = 0;
			ncols = 0;
			int a;
			int b;
			float f;
			while (ifs.good()) {
				ifs >> a >> b >> f;
				if (a > nrows - 1)
					nrows = a + 1;
				if (b > ncols - 1)
					ncols = b + 1;
			}
			ifs.close();
			ifs.open(filename);
			return;
		}
		cout << "Error: could not parse first line" << endl;
	}

	template<class MATRIX>
	ASCII(const char* filename, const MATRIX& M) {
		cout << M.str() << endl;

		ofstream f;
		f.open(filename);
		f << M.ncols << " " << M.nrows << endl;
		for (auto it = M.begin(); it != M.end(); it++)
			f << it.i() << " " << it.j() << " " << *it << endl;
		f.close();
	}

	SparseMatrixFile::iterator begin() {
		ifs.seekg(0);
		int dummy1;
		int dummy2;
		cout << "---" << dummy1 << "-" << dummy2 << endl;
		ifs >> dummy1 >> dummy2;
		return SparseMatrixFile::iterator(*this);
	}

	void get(SparseMatrixFile::iterator& it) {
		if (!ifs.good()) {
			it.endflag = true;
			return;
		}
		it.endflag = false;
		ifs >> it.obj.i >> it.obj.j >> it.obj.value;
	}

	SparseMatrixFile& operator>>(IndexValueTriple& dest) {
		if (!ifs.good()) {
			dest.i = -1;
			return *this;
		}
		ifs >> dest.i >> dest.j >> dest.value;
		return *this;
	}

};

#endif
