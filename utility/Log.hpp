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



#ifndef _Log
#define _Log

#include <chrono>
#include "MMFbase.hpp"

class Log {
public:

	Log() {
		startClock();
	}

public:

	Log& operator<<(const string& s);
	Log& operator<<(const char* s);

	Log& skip(const int v = 0, const int n = 1) {
		if (verbosity < v)
			return *this;
		if (skippedlines >= n)
			return *this;
		for (int i = 0; i < n - skippedlines; i++)
			cout << endl;
		skippedlines = n;
		return *this;
	}

	Log& log(const int v, const string& s);
	Log& log(const int v, const char* s, const int i);
	Log& log(const int v, const char* s1, const int i, const char* s2);
	Log& log(const int v, const char* s1, const int i1, const char* s2,
			const int i2, const char* s3);
	Log& log(const int v, const char* s1, const int i1, const char* s2,
			const int i2, const char* s3, const int i3, const char* s4);
	Log& log(const int v, const char* s1, const int i1, const char* s2,
			const int i2, const char* s3, const int i3, const char* s4,
			const int i4, const char* s5);
	Log& log(const int v, const char* s, const double f);
	Log& log(const int v, const int i, const double f);

	void startClock(const int i = 0);
	double clock(const int i = 0);

public:
	chrono::time_point<chrono::system_clock> t;

	int verbosity = 0;
	int skippedlines = 0;
};

#endif
