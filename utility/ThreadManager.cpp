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



#include "ThreadManager.hpp"
#include "ThreadBank.hpp"

void ThreadManager::enqueue(ThreadBank* bank) {
	lock_guard < mutex > lock(mx);
	if (is_runnable(bank))
		launch(bank);
	else
		queue.push_back(bank);
}

void ThreadManager::release(ThreadBank* bank) {
	lock_guard < mutex > lock(mx);
	if (bank->nprivileged > 0)
		bank->nprivileged--;
	else
		nthreads--;
	for (auto it = queue.begin(); it != queue.end(); it++)
		if (is_runnable(*it)) {
			launch(*it);
			it = queue.erase(it);
		}
}

bool ThreadManager::is_runnable(ThreadBank* bank) {
	return bank->is_ready()
			&& (bank->nprivileged < bank->maxprivileged || nthreads < maxthreads);
}

void ThreadManager::launch(ThreadBank* bank) {
	if (bank->nprivileged < bank->maxprivileged)
		bank->nprivileged++;
	else
		nthreads++;
	bank->gate.unlock();
}
