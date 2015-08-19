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



#include <iostream>
#include <chrono>

#include "ThreadBank.hpp"
#include "MMFglobal.hpp"

void parent(const int i) {
	{
		CoutLock lock;
		cout << "Parent " << i << " starting." << endl;
	}
	{
		ThreadBank threads(7);
		uniform_int_distribution<int> distr(0, 500);
		this_thread::sleep_for(
				chrono::milliseconds(distr(randomNumberGenerator)));
		for (int j = 0; j < 10; j++) {
			threads.add(
					[i,&threads](int _j) {
						{	CoutLock lock;
							cout<<"  Child "<<i<<":"<<_j<<" starting."<<endl;
							cout<<"     "<<threads.nthreads-threads.nprivileged;
							cout<<"+"<<threads.nprivileged;
							cout<<"("<<threadManager.get_nthreads()<<") threads running"<<endl;}
						uniform_int_distribution<int> distr(0,500);
						this_thread::sleep_for(chrono::milliseconds(distr(randomNumberGenerator)));
						{	CoutLock lock; cout<<"  Child "<<i<<":"<<_j<<" done."<<endl;}
					}, j);
		}
	}
	{
		CoutLock lock;
		cout << "Parent " << i << " done." << endl;
	}
}

int main(int argc, char** argv) {

	ThreadBank threads(7);
	for (int i = 0; i < 7; i++)
		threads.add(parent, i);

}
