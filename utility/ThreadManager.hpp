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



#ifndef _ThreadManager
#define _ThreadManager

#include <list>
#include <queue>
#include <mutex>

class ThreadBank;

using namespace std;

class ThreadManager {
public:

	ThreadManager(const int _maxthreads) :
			maxthreads(_maxthreads), nthreads(0) {
	}
	~ThreadManager() {
	}

public:

	void enqueue(ThreadBank* bank);
	void release(ThreadBank* bank);

	int get_nthreads() {
		lock_guard < mutex > lock(mx);
		return nthreads;
	}

private:

	bool is_runnable(ThreadBank* bank);
	void launch(ThreadBank* bank);

public:

	int maxthreads;

private:

	mutex mx;
	int nthreads;
	list<ThreadBank*> queue;

};

#endif
