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
#include "MMFglobal.hpp"

int main(int argc, char** argv) {
	Cvector v(6, Random());
	Vectorv v1(v);
	Vectorl v2(v);
	Vectorh v3(v);
	cout << v.str(Dense()) << endl;
	cout << v1.str(Dense()) << endl;
	cout << v2.str(Dense()) << endl;
	cout << v3.str(Dense()) << endl;
	cout << Vectorv(v1) << endl;

	for (auto& p : v1)
		cout << p.first << " " << p.second << endl;
	cout << endl;
	for (auto& p : v2)
		cout << p.first << " " << p.second << endl;
	cout << endl;
	for (auto& p : v3)
		cout << p.first << " " << p.second << endl;
	cout << endl;

	v1.insert(3, 74);
	v1(2) = 0;
	cout << v1 << endl;
	v1.tidy();
	cout << v1 << endl;

	v1.foreach([](int i, FIELD& v) {cout<<"ff:"<<i<<" "<<v<<endl;});
}
