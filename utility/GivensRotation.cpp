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



#include "GivensRotation.hpp"
#include <random>
#include <math.h>

extern default_random_engine randomNumberGenerator;

GivensRotation::GivensRotation(int _i1, int _i2, double _alpha) {
	if (_i1 < _i2) {
		i1 = _i1;
		i2 = _i2;
		cos = ::cos(_alpha);
		sin = ::sin(_alpha);
	} else {
		i2 = _i1;
		i1 = _i2;
		cos = ::cos(_alpha);
		sin = -::sin(_alpha);
	}
}

GivensRotation::GivensRotation(const Random dummy, const int n) {
	uniform_int_distribution<int> distri(0, n - 1);
	uniform_real_distribution<FIELD> distrr(0, 1);
	i1 = distri(randomNumberGenerator);
	i2 = i1;
	while (i2 == i1) {
		i2 = distri(randomNumberGenerator);
	}
	if (i2 > i1)
		swapp(i1, i2);
	cos = distrr(randomNumberGenerator);
	sin = sqrt(1 - cos * cos);
}

string GivensRotation::str() const {
	ostringstream stream;
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	stream << "Rotation on (" << i1 << "," << i2 << "):" << endl;
	stream << "[ ";
	stream.width(6);
	stream << cos << " ";
	stream.width(6);
	stream << -sin << "  ]" << endl;
	stream << "[ ";
	stream.width(6);
	stream << sin << " ";
	stream.width(6);
	stream << cos << "  ]" << endl;
	return stream.str();
}

void GivensRotation::serialize(Rstream& s) const {
	s << "GivensRotation{" << Rstream::endl;
	s.var("i1", i1);
	s.var("i2", i2);
	s.var("cos", cos);
	s.var("sin", sin);
	s << "}" << Rstream::endl;
}

void GivensRotation::serialize(Bofstream& ofs) const {
	ofs.write(i1);
	ofs.write(i2);
	ofs.write(cos);
	ofs.write(sin);
}

GivensRotation::GivensRotation(Bifstream& ifs) {
	ifs.read(i1);
	ifs.read(i2);
	ifs.read(cos);
	ifs.read(sin);
}
