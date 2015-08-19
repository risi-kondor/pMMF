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



#include "Cvector.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "EigenAdaptors.hpp"

extern default_random_engine randomNumberGenerator;

Cvector::Cvector(const int _n, const class Random& dummy) :
		DenseVector(_n) {
	array = new FIELD[n];
	uniform_real_distribution<FIELD> distr;
	for (int i = 0; i < n; i++)
		array[i] = distr(randomNumberGenerator);
}

Cvector Cvector::Random(const int n) {
	Cvector v(n);
	uniform_real_distribution<FIELD> distr;
	for (int i = 0; i < n; i++)
		v.array[i] = distr(randomNumberGenerator);
	return v;
}

Cvector::Cvector(const initializer_list<FIELD> list) :
		DenseVector(list.size()) {
	array = new FIELD[n];
	int i = 0;
	for (FIELD v : list)
		array[i++] = v;
}

Cvector::Cvector(const int _n, const FIELD* _array) :
		DenseVector(_n) {
	array = new FIELD[n];
	for (int i = 0; i < n; i++)
		array[i] = _array[i];
}

// ---- Conversions ---------------------------------------------------------------------------------------------
Cvector::Cvector(const Vectorv& x) :
		DenseVector(x.n) {
	array = new FIELD[n];
	for (int i = 0; i < n; i++)
		array[i] = 0;
	for (auto& it : x)
		array[it.first] = it.second;
}

Cvector::Cvector(const Vectorl& x) :
		DenseVector(x.n) {
	array = new FIELD[n];
	for (int i = 0; i < n; i++)
		array[i] = 0;
	for (auto& it : x)
		array[it.first] = it.second;
}

Cvector::Cvector(const Vectorh& x) :
		DenseVector(x.n) {
	array = new FIELD[n];
	for (int i = 0; i < n; i++)
		array[i] = 0;
	for (auto& it : x)
		array[it.first] = it.second;
}

Cvector::Cvector(const EigenVectorXdAdaptor& x) :
		Cvector(x.size()) {
	for (int i = 0; i < n; i++)
		array[i] = x(i);
}

template<>
Eigen::VectorXd Cvector::convert() const {
	Eigen::VectorXd v(n);
	for (int i = 0; i < n; i++)
		v(i) = array[i];
	return v;
}

Cvector Cvector::merge(const Cvector& x, const Cvector& y) {
	Cvector r(x.n + y.n);
	for (int i = 0; i < x.n; i++)
		r.array[i] = x.array[i];
	for (int i = 0; i < y.n; i++)
		r.array[i + x.n] = y.array[i];
	return r;
}

// ---- I/O ------------------------------------------------------------------------------------------------------
ostream& operator<<(ostream& stream, const Cvector& x) {
	stream << "(";
	for (int i = 0; i < x.n - 1; i++)
		stream << x.array[i] << ",";
	stream << x.array[x.n - 1] << ")";
	return stream;
}

void Cvector::serialize(Rstream& rstream) const {
	rstream << "Cvector{" << Rstream::endl;
	rstream.var("n", n);
	rstream.var("array", " *OMITTED*");
	rstream << "}" << Rstream::endl;
}

void Cvector::serialize(Bofstream& ofs) const {
	ofs.tag("Cvector", 0);
	ofs.write(n);
	ofs.write_array(array, n);
}

Cvector::Cvector(Bifstream& ifs) :
		DenseVector(0) {
	ifs.check("Cvector", 0);
	ifs.read(n);
	ifs.read_array(array);
}
