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



#include "KpointRotation.hpp"
#include "Cmatrix.hpp"

KpointRotation::KpointRotation(const int _k, const Identity dummy) :
		KpointRotation(_k) {
	for (int i = 0; i < k; i++)
		ix[i] = i;
	for (int i = 0; i < k; i++)
		for (int j = 0; j < k; j++)
			if (i == j)
				q[j * k + i] = 1;
			else
				q[j * k + i] = 0;
}

KpointRotation::KpointRotation(const int _k, const Random dummy, const int n) :
		KpointRotation(_k) {
} // unimplemented

KpointRotation::KpointRotation(const IndexSet& I, const Cmatrix& M) :
		k(I.k) {
	ix = new int[k];
	for (int i = 0; i < k; i++)
		ix[i] = I.ix[i];
	q = new double[k * k];
	for (int i = 0; i < k; i++)
		for (int j = 0; j < k; j++)
			q[j * k + i] = M.array[i * k + j];
	temp = new double[k];
	tempp = new FIELD*[k];
}

string KpointRotation::str() const {
	ostringstream stream;
	stream << "Rotation on (";
	for (int i = 0; i < k - 1; i++)
		stream << ix[i] << ",";
	stream << ix[k - 1] << "):\n";
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	for (int i = 0; i < k; i++) {
		stream << "[ ";
		for (int j = 0; j < k; j++) {
			stream.width(6);
			stream << q[j * k + i] << " ";
		}
		stream << " ]\n";
	}
	return stream.str();
}

void KpointRotation::serialize(Rstream& s) const {
	s << "KpointRotation{" << Rstream::endl;
	s.var("k", k);
	s << "}" << Rstream::endl;
}

void KpointRotation::serialize(Bofstream& ofs) const {
	ofs.write(k);
	ofs.write_array(ix, k);
	ofs.write_array(q, k * k);
}

KpointRotation::KpointRotation(Bifstream& ifs) {
	ifs.read(k);
	ifs.read_array(ix);
	ifs.read_array(q);
	temp = new double[k];
	tempp = new FIELD*[k];
}
