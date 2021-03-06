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



#include "Remap.hpp"
#include <random>

extern default_random_engine randomNumberGenerator;

Remap::Remap(const int _n, const Random random) :
		Remap(_n) {
	for (int i = 0; i < n - 1; i++) {
		uniform_int_distribution<int> distr(i + 1, n - 1);
		int j = distr(randomNumberGenerator);
		swap(i, j);
	}
}

string Remap::str() const {
	ostringstream stream;
	for (int i = 0; i < n; i++)
		stream << i << " -> " << forward[i] << endl;
	return stream.str();
}
