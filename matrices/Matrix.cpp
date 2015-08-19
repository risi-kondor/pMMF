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



#include "Matrix.hpp"
#include "Cmatrix.hpp"

Cmatrix Matrix::dot(const Matrix& x) const {
	return Cmatrix(0, 0);
}

string Matrix::str(const Dense dummy) const {
	ostringstream stream;
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	for (int i = 0; i < nrows; i++) {
		stream << "[ ";
		for (int j = 0; j < ncols; j++) {
			stream.width(6);
			stream << this->read(i, j) << " ";
		}
		stream << " ]\n";
	}
	return stream.str();
}

string Matrix::str(const Sparse dummy) const {
	ostringstream stream;
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++)
			if ((*this)(i, j) != 0)
				stream << "(" << i << "," << j << ") : " << (*this)(i, j)
						<< endl;
	return stream.str();
}

ostream& operator<<(ostream& stream, const Matrix& x) {
	stream << x.str();
	return stream;
}
