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



#ifndef _Rstream
#define _Rstream

#include "MMFbase.hpp"

class Rstream {
public:

	Rstream(ostream& _out, const int _depth = 16) :
			out(_out), indent(0), depth(_depth), bol(true) {
	}

	~Rstream() {
		out << std::endl;
	}
	template<class T>
	Rstream& operator<<(const T& x) {
		if (bol) {
			for (int i = 0; i < indent; i++)
				out << "  ";
			bol = false;
		}
		out << x;
		return *this;
	}

	typedef Rstream& (*RstreamManipulator)(Rstream&);
	Rstream& operator<<(const RstreamManipulator& manip) {
		return manip(*this);
	}

	static Rstream& endl(Rstream& x) {
		x.out << std::endl;
		x.bol = true;
		return x;
	}

	template<class T>
	Rstream& write(const T& x) {
		if (depth < 0) {
			out << std::endl;
			bol = true;
			return *this;
		}
		indent++;
		depth--;
		x.serialize(*this);
		indent--;
		depth++;
		return *this;
	}

	template<class T>
	const Rstream& var(const char* name, const T& x) {
		if (bol) {
			for (int i = 0; i < indent; i++)
				out << "  ";
			bol = false;
		}
		out << "  " << name << "=" << x << std::endl;
		bol = true;
		return *this;
	}

public:

	int indent;
	int depth;
	mutable bool bol;

	ostream& out;

};

#endif
