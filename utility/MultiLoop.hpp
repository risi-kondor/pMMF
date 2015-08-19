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



#ifndef _MultiLoop
#define _MultiLoop

#include "MMFbase.hpp"
#include "ThreadBank.hpp"

extern bool multithreading;

class MultiLoop {
public:

	MultiLoop(const int n, std::function<void(int)> lambda) {
		if (multithreading) {
			ThreadBank threads(1000);
			for (int i = 0; i < n; i++)
				threads.add(lambda, i);
		} else
			for (int i = 0; i < n; i++)
				lambda(i);
	}

	MultiLoop(const int n1, const int n2,
			std::function<void(int, int)> lambda) {
		if (multithreading) {
			for (int i1 = 0; i1 < n1; i1++) {
				ThreadBank threads(1000);
				for (int i2 = 0; i2 < n2; i2++)
					threads.add(lambda, i1, i2);
			}
		} else
			for (int i1 = 0; i1 < n1; i1++)
				for (int i2 = 0; i2 < n2; i2++)
					lambda(i1, i2);
	}

};

template<class OBJ>
class MultiForeach {
public:

	template<class CONTAINER>
	MultiForeach(CONTAINER& V, std::function<void(OBJ)> lambda) {
		if (multithreading) {
			vector < thread > workers;
			for (auto& x : V)
				workers.push_back(thread(lambda, x));
			for (auto& w : workers)
				w.join();
		} else
			for (auto& x : V)
				lambda(x);
	}

};

template<class ARG1>
class MultiLoop1 { // should be able to do it with context
public:
	MultiLoop1(const int n, std::function<void(int, ARG1)> lambda, ARG1 arg1) {
		if (multithreading) {
			ThreadBank threads(1000);
			for (int i = 0; i < n; i++)
				threads.add(lambda, i, arg1);
		} else
			for (int i = 0; i < n; i++)
				lambda(i, arg1);
	}
};

template<class ARG1, class ARG2>
class MultiLoop2 { // should be able to do it with context
public:
	MultiLoop2(const int n, std::function<void(int, ARG1, ARG2)> lambda,
			ARG1 arg1, ARG2 arg2) {
		if (multithreading) {
			ThreadBank threads(1000);
			for (int i = 0; i < n; i++)
				threads.add(lambda, i, arg1, arg2);
		} else
			for (int i = 0; i < n; i++)
				lambda(i, arg1, arg2);
	}
};

#endif

