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



#include "Log.hpp"
#include <iostream>
#include <mutex>

//extern mutex cout_mutex;

Log& Log::operator<<(const string& s) {
	CoutLock lock;
	cout << s << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::operator<<(const char* s) {
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const string& s) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const char* s, const int i) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s << i << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const char* s1, const int i, const char* s2) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s1 << i << s2 << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const char* s1, const int i1, const char* s2,
		const int i2, const char* s3) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s1 << i1 << s2 << i2 << s3 << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const char* s1, const int i1, const char* s2,
		const int i2, const char* s3, const int i3, const char* s4) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s1 << i1 << s2 << i2 << s3 << i3 << s4 << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const char* s1, const int i1, const char* s2,
		const int i2, const char* s3, const int i3, const char* s4,
		const int i4, const char* s5) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s1 << i1 << s2 << i2 << s3 << i3 << s4 << i4 << s5
			<< endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const char* s, const double f) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << "s: " << s << f << endl;
	skippedlines = 0;
	return *this;
}

Log& Log::log(const int v, const int i, const double f) {
	if (verbosity < v)
		return *this;
	CoutLock lock;
	double now =
			chrono::duration<double>(chrono::system_clock::now() - t).count();
	cout.precision(5);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << now << " " << i << " " << f << " " << endl;
	skippedlines = 0;
	return *this;
}

void Log::startClock(const int i) {
	t = chrono::system_clock::now();
}

double Log::clock(const int i) {
	return chrono::duration<double>(chrono::system_clock::now() - t).count();
}
