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



#ifndef _DenseMatrixFileASCII
#define _DenseMatrixFileASCII

#include "DenseMatrixFile.hpp"

class DenseMatrixFile::ASCII: public DenseMatrixFile {
public:

	ASCII(const char* filename) {
		ifs.open(filename);
		if (ifs.fail()) {
			cout << "Failed to open " << filename << "." << endl;
			return;
		}
		char buffer[1024];
		int linelength = 0;
		do {
			ifs.get(buffer, 1024);
			linelength += strlen(buffer);
		} while (strlen(buffer) > 0);
		cout << "Line length=" << linelength << endl;
		ifs.close();
		ifs.open(filename); // why does ifs.seekg(0) not work?
		float b;
		ncols = 0;
		while (ifs.good() && ifs.tellg() < linelength) {
			ifs >> b;
			ncols++;
		}
		nrows = 0;
		while (ifs.good()) {
			for (int i = 0; i < ncols; i++)
				ifs >> b;
			nrows++;
		}
		cout << nrows << " " << ncols << endl;
		ifs.close();
		ifs.open(filename);
	}

	DenseMatrixFile::iterator begin() {
		ifs.seekg(0);
		return DenseMatrixFile::iterator(*this);
	}

	void get(DenseMatrixFile::iterator& it) {
		if (!ifs.good()) {
			it.endflag = true;
			return;
		}
		it.endflag = false;
		ifs >> it.value;
	}

	DenseMatrixFile& operator>>(FIELD& dest) {
		if (ifs.good())
			ifs >> dest;
		return *this;
	}

};
// note: there is a potential problem when the file marker moves to eof, because seekg cannot reset it to beginning
#endif
