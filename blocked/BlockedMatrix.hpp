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



#ifndef _BlockedMatrix
#define _BlockedMatrix

#include "Tower.hpp"
#include "Street.hpp"
#include "BlockedRemap.hpp"
#include "BlockedVector.hpp"
#include "MatrixX.hpp"
#include "MultiLoop.hpp"

extern Log mlog;

template<class BLOCK>
class BlockedMatrix {
public:

	BlockedMatrix(const BlockedMatrix& X);
	BlockedMatrix(BlockedMatrix&& X);
	BlockedMatrix<BLOCK>& operator=(const BlockedMatrix& X);
	BlockedMatrix<BLOCK>& operator=(BlockedMatrix&& X);
	BlockedMatrix<BLOCK> clone() const;
	~BlockedMatrix() {
		if (!nodestruct)
			for (int i = 0; i < nstreets * ntowers; i++)
				delete block[i];
		delete[] block;
	}
	bool nodestruct = false;

public:

	BlockedMatrix() :
			BlockedMatrix(0, 0) {
	}
	;
	BlockedMatrix(const int _nstreets, const int _ntowers,
			const bool _nodestruct = false);
	BlockedMatrix(const Bstructure& rstruct, const Bstructure& cstruct);
	BlockedMatrix(const Bstructure& rstruct, const Bstructure& cstruct,
			const Zero& dummy);

	BlockedMatrix(const BLOCK& M);
	BlockedMatrix(BLOCK&& M);
	BlockedMatrix(BLOCK& M, bool nocopy = false);

	BlockedMatrix(const Cmatrix& M, const Bstructure& rstruct,
			const Bstructure& cstruct);
	BlockedMatrix(const MatrixXv& M, const Bstructure& rstruct,
			const Bstructure& cstruct);
	BlockedMatrix(const MatrixXl& M, const Bstructure& rstruct,
			const Bstructure& cstruct);
	BlockedMatrix(const MatrixXh& M, const Bstructure& rstruct,
			const Bstructure& cstruct);

	template<class TYPE>
	TYPE convert() const;

public:
	// named constructors
	static BlockedMatrix Identity(const Bstructure& structure);
	static BlockedMatrix RandomSymmetric(const Bstructure& structure);
	static BlockedMatrix Random(const Bstructure& rstruct,
			const Bstructure& cstruct);

public:
	// indexing
	BLOCK& operator()(const int I, const int J) {
		assert(I < nstreets);
		assert(J < ntowers);
		return *block[J * nstreets + I];
	}
	const BLOCK& operator()(const int I, const int J) const {
		assert(I < nstreets);
		assert(J < ntowers);
		return *block[J * nstreets + I];
	}

	BLOCK*& ptr(const int I, const int J) {
		assert(I < nstreets);
		assert(J < ntowers);
		return block[J * nstreets + I];
	}
	const BLOCK* ptr(const int I, const int J) const {
		assert(I < nstreets);
		assert(J < ntowers);
		return block[J * nstreets + I];
	}

	FIELD& operator()(const int I, const int i, const int J, const int j) {
		return (*block[J * nstreets + I])(i, j);
	}
	FIELD operator()(const int I, const int i, const int J, const int j) const {
		return (*block[J * nstreets + I])(i, j);
	}
	FIELD read(const int I, const int i, const int J, const int j) const {
		return block[J * nstreets + I]->read(i, j);
	}
	FIELD isFilled(const int I, const int i, const int J, const int j) const {
		assert(I < nstreets);
		assert(J < ntowers);
		return block[J * nstreets + I]->isFilled(i, j);
	}

	int nnz() const {
		int t = 0;
		for (int i = 0; i < nstreets * ntowers; i++)
			t += block[i]->nnz();
		return t;
	}

public:
	// virtual views and info
	Tower<BLOCK> virtual_tower(const int J);
	const Tower<BLOCK> virtual_tower(const int J) const;
	Street<BLOCK> virtual_street(const int I);
	const Street<BLOCK> virtual_street(const int I) const;

	Bstructure rowstructure() const;
	Bstructure colstructure() const;
	int nrows() const;
	int ncols() const;

	int findI(int& i, const Bstructure& bstruct) const {
		int I = 0;
		while (i >= bstruct[I]) {
			i -= bstruct[I];
			I++;
		}
		return I;
	}

public:
	// operations
	template<class VECTOR> BlockedVector<VECTOR> operator*(
			const BlockedVector<VECTOR>& x) const;
	template<class BLOCK2> BlockedMatrix<BLOCK> operator*(
			BlockedMatrix<BLOCK2>& X) const;

	template<class BLOCK2> FIELD diff2(const BlockedMatrix<BLOCK2>& X) const;

public:
	// solvers
	BlockedMatrix<BLOCK> inverse() const;
	BlockedMatrix<BLOCK> inverseSqrt() const;

public:
	// reblocking
	BlockedMatrix(const BlockedMatrix& X, const class Identity& dummy,
			const BlockedRemap& cmap, const bool inverse = false);
	BlockedMatrix(const BlockedMatrix& X, const BlockedRemap& rmap,
			const class Identity& dummy, const bool inverse = false);
	template<class BLOCK2>
	BlockedMatrix(const BlockedMatrix<BLOCK2>& X, const BlockedRemap& rmap,
			const BlockedRemap& cmap, const bool inverse = false,
			const bool push = false);

	BlockedMatrix<BLOCK> pullColumns(const BlockedRemap& cmap,
			const bool inverse = 0) const;

public:
	// I/O
	string str(const Dense dummy) const;
	string str(const Sparse dummy) const;
	string str() const {
		return str(Dense());
	}
	;

	BlockedMatrix(Bifstream& ifs);
	void serialize(Bofstream& ofs) const;
	void serialize(Rstream& rstream) const;

	BlockedMatrix(DenseMatrixFile& file, const BlockStructure& rstruct,
			const BlockStructure& cstruct);
	BlockedMatrix(SparseMatrixFile& file, const BlockStructure& rstruct,
			const BlockStructure& cstruct);
	BlockedMatrix(DenseMatrixFile& file, const int _nstreets = 1,
			const int _ntowers = 1);
	BlockedMatrix(SparseMatrixFile& file, const int _nstreets = 1,
			const int _ntowers = 1);

public:

	int nstreets;
	int ntowers;
	BLOCK** block;

};

template<class BLOCK>
ostream& operator<<(ostream& stream, const BlockedMatrix<BLOCK>& x) {
	stream << x.str();
	return stream;
}

// ---- Constructors ---------------------------------------------------------------------------------------------

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const int _nstreets, const int _ntowers,
		const bool _nodestruct) :
		nstreets(_nstreets), ntowers(_ntowers), nodestruct(_nodestruct) {
	block = new BLOCK*[nstreets * ntowers];
	for (int i = 0; i < nstreets * ntowers; i++)
		block[i] = nullptr;
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const Bstructure& rstruct,
		const Bstructure& cstruct) :
		BlockedMatrix(rstruct.size(), cstruct.size()) {
	for (int i = 0; i < nstreets; i++)
		for (int j = 0; j < ntowers; j++)
			block[j * nstreets + i] = new BLOCK(rstruct[i], cstruct[j]);
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const Bstructure& rstruct,
		const Bstructure& cstruct, const Zero& dummy) :
		BlockedMatrix(rstruct.size(), cstruct.size()) {
	for (int i = 0; i < nstreets; i++)
		for (int j = 0; j < ntowers; j++)
			block[j * nstreets + i] = new BLOCK(rstruct[i], cstruct[j], dummy);
}

template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::Random(const Bstructure& rstruct,
		const Bstructure& cstruct) {
	int nstr = rstruct.size();
	int ntow = cstruct.size();
	BlockedMatrix M(nstr, ntow);
	for (int i = 0; i < nstr; i++)
		for (int j = 0; j < ntow; j++)
			M.ptr(i, j) = new BLOCK(BLOCK::Random(rstruct[i], cstruct[j]));
	return M;
}


template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::Identity(
		const Bstructure& structure) {
	int nblocks = structure.size();
	BlockedMatrix M(nblocks, nblocks);
	for (int i = 0; i < nblocks; i++)
		for (int j = 0; j < nblocks; j++)
			if (i == j)
				M.ptr(i, j) = new BLOCK(BLOCK::Identity(structure[i]));
			else
				M.ptr(i, j) = new BLOCK(structure[i], structure[j], Zero());
	return M;
}

template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::RandomSymmetric(
		const Bstructure& structure) {
	int nblocks = structure.size();
	BlockedMatrix M(nblocks, nblocks);
	for (int i = 0; i < nblocks; i++) {
		for (int j = 0; j < i; j++) {
			M.ptr(i, j) = new BLOCK(BLOCK::Random(structure[i], structure[j]));
			M.ptr(j, i) = new BLOCK(*M.ptr(i, j), Transpose());
		}
		M.ptr(i, i) = new BLOCK(BLOCK::RandomSymmetric(structure[i]));
	}
	return M;
}

// ---- Copying --------------------------------------------------------------------------------------------------

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const BlockedMatrix<BLOCK>& X) :
		BlockedMatrix(X.nstreets, X.ntowers) {
#ifdef _BLOCKEDCOPYWARNING 
	cout << "WARNING: BlockedMatrix copied." << endl;
#endif
	MultiLoop(nstreets * ntowers,
			[this,&X](const int i) {block[i]=new BLOCK(*X.block[i]);});
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(BlockedMatrix<BLOCK> && X) :
		BlockedMatrix(X.nstreets, X.ntowers, X.nodestruct) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedMatrix moved." << endl;
#endif
	for (int i = 0; i < nstreets * ntowers; i++) {
		block[i] = X.block[i];
		X.block[i] = nullptr;
	};
}

template<class BLOCK>
BlockedMatrix<BLOCK>& BlockedMatrix<BLOCK>::operator=(
		const BlockedMatrix<BLOCK>& X) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedMatrix assigned" << endl;
#endif
	if (!nodestruct)
		for (int i = 0; i < nstreets * ntowers; i++)
			delete block[i];
	delete[] block;
	ntowers = X.ntowers;
	nstreets = X.nstreets;
	nodestruct = false;
	block = new BLOCK*[nstreets * ntowers];
	MultiLoop(nstreets * ntowers,
			[this,X](const int i) {block[i]=new BLOCK(*X.block[i]);});
	return *this;
}

template<class BLOCK>
BlockedMatrix<BLOCK>& BlockedMatrix<BLOCK>::operator=(
		BlockedMatrix<BLOCK> && X) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedMatrix move-assigned" << endl;
#endif
	if (!nodestruct)
		for (int i = 0; i < nstreets * ntowers; i++)
			delete block[i];
	delete[] block;
	ntowers = X.ntowers;
	nstreets = X.nstreets;
	nodestruct = X.nodestruct;
	block = new BLOCK*[nstreets * ntowers];
	for (int i = 0; i < nstreets * ntowers; i++) {
		block[i] = X.block[i];
		X.block[i] = nullptr;
	};
	return *this;
}

template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::clone() const {
#ifdef _BLOCKEDCOPYWARNING
	//  cout<<"WARNING: BlockedMatrix cloned."<<endl;
#endif
	BlockedMatrix M(nstreets, ntowers, true);
	for (int i = 0; i < nstreets * ntowers; i++)
		M.block[i] = block[i];
	return M;
}

// ---- Converters -----------------------------------------------------------------------------------------------

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const BLOCK& M) :
		BlockedMatrix(1, 1) {
	block[0] = new BLOCK(M);
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(BLOCK&& M) :
		BlockedMatrix(1, 1) {
	block[0] = new BLOCK(M);
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(BLOCK& M, bool nocopy) :
		BlockedMatrix(1, 1) {
	if (nocopy) {
		block[0] = &M;
		nodestruct = true;
	} else
		block[0] = new BLOCK(M);
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const MatrixXv& X,
		const Bstructure& rstruct, const Bstructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	int joffset = 0;
	for (int tower = 0; tower < ntowers; tower++) {
		for (int j = 0; j < cstruct[tower]; j++) {
			assert(joffset + j < X.ncols);
			auto& col = *X.column[joffset + j];
			col.sort();
			auto it = col.begin();
			int ioffset = 0;
			for (int street = 0; street < nstreets; street++) {
				while (it != col.end() && it->first < ioffset + rstruct[street]) {
					block[tower * nstreets + street]->column[j]->append(
							it->first - ioffset, it->second);
					it++;
				}
				ioffset += rstruct[street];
			}
		}
		joffset += cstruct[tower];
	}
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const MatrixXl& X,
		const Bstructure& rstruct, const Bstructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	int joffset = 0;
	for (int tower = 0; tower < ntowers; tower++) {
		for (int j = 0; j < cstruct[tower]; j++) {
			assert(joffset + j < X.ncols);
			auto& col = *X.column[joffset + j];
			auto it = col.begin();
			int ioffset = 0;
			for (int street = 0; street < nstreets; street++) {
				while (it != col.end() && it->first < ioffset + rstruct[street]) {
					block[tower * nstreets + street]->column[j]->append(
							it->first - ioffset, it->second);
					it++;
				}
				ioffset += rstruct[street];
			}
		}
		joffset += cstruct[tower];
	}
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const MatrixXh& X,
		const Bstructure& rstruct, const Bstructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	int joffset = 0;
	for (int tower = 0; tower < ntowers; tower++) {
		for (int j = 0; j < cstruct[tower]; j++) {
			assert(joffset + j < X.ncols);
			auto& col = *X.column[joffset + j];
			auto it = col.begin();
			int ioffset = 0;
			for (int street = 0; street < nstreets; street++) {
				while (it != col.end() && it->first < ioffset + rstruct[street]) {
					block[tower * nstreets + street]->column[j]->insert(
							it->first - ioffset, it->second);
					it++;
				}
				ioffset += rstruct[street];
			}
		}
		joffset += cstruct[tower];
	}
}

// ---- Virtual views and info -----------------------------------------------------------------------------------


template<class BLOCK>
Tower<BLOCK> BlockedMatrix<BLOCK>::virtual_tower(const int J) {
	assert(J < ntowers);
	Tower<BLOCK> T(nstreets, block[J * nstreets]->ncols, true);
	for (int I = 0; I < nstreets; I++)
		T.block[I] = block[J * nstreets + I];
	return T;
}

template<class BLOCK>
const Tower<BLOCK> BlockedMatrix<BLOCK>::virtual_tower(const int J) const {
	assert(J < ntowers);
	Tower<BLOCK> T(nstreets, block[J * nstreets]->ncols, true);
	for (int I = 0; I < nstreets; I++)
		T.block[I] = block[J * nstreets + I];
	return T;
}

template<class BLOCK>
Street<BLOCK> BlockedMatrix<BLOCK>::virtual_street(const int I) {
	assert(I < nstreets);
	Street<BLOCK> S(ntowers, block[I]->nrows, true);
	for (int J = 0; J < ntowers; J++)
		S.block[J] = block[J * nstreets + I];
	return S;
}

template<class BLOCK>
const Street<BLOCK> BlockedMatrix<BLOCK>::virtual_street(const int I) const {
	assert(I < nstreets);
	Street<BLOCK> S(ntowers, block[I]->nrows, true);
	for (int J = 0; J < ntowers; J++)
		S.block[J] = block[J * nstreets + I];
	return S;
}

template<class BLOCK>
Bstructure BlockedMatrix<BLOCK>::rowstructure() const {
	Bstructure result;
	for (int i = 0; i < nstreets; i++) {
		assert(block[i] != NULL);
		result.push_back(block[i]->nrows);
	}
	return result;
}

template<class BLOCK>
Bstructure BlockedMatrix<BLOCK>::colstructure() const {
	Bstructure result;
	for (int i = 0; i < ntowers; i++) {
		assert(ptr(0, i) != NULL);
		result.push_back(ptr(0, i)->ncols);
	}
	return result;
}

template<class BLOCK>
int BlockedMatrix<BLOCK>::nrows() const {
	int result = 0;
	for (int i = 0; i < nstreets; i++) {
		assert(block[i] != NULL);
		result += block[i]->nrows;
	}
	return result;
}

template<class BLOCK>
int BlockedMatrix<BLOCK>::ncols() const {
	int result = 0;
	for (int i = 0; i < ntowers; i++) {
		assert(ptr(0, i) != NULL);
		result += ptr(0, i)->ncols;
	}
	return result;
}

// ---- Vector Operations ---------------------------------------------------------------------------------------

template<class BLOCK>
template<class VECTOR>
BlockedVector<VECTOR> BlockedMatrix<BLOCK>::operator*(
		const BlockedVector<VECTOR>& x) const {
	assert(ntowers == x.nblocks);
	BlockedVector<VECTOR> r(nstreets);
	MultiLoop(nstreets, [this,&r,&x](const int i) {
		r.block[i]=new VECTOR((*block[i])*(*x.block[0]));
		for(int J=1; J<ntowers; J++)
		r.block[i]->add((*block[J*nstreets+i])*(*x.block[J]));});
	return r;
}

// ---- Matrix Operations ----------------------------------------------------------------------------------------

template<class BLOCK>
template<class BLOCK2>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::operator*(
		BlockedMatrix<BLOCK2>& X) const {
	BlockedMatrix<BLOCK> M(nstreets, X.ntowers);
	assert(ntowers == X.nstreets);
	mlog.log(4, "Blocked matrix multiplication:[", nstreets, "x", ntowers,
			"] by [", ntowers, "x", X.ntowers, "] ...");
	MultiLoop(nstreets, X.ntowers, [this,&M,&X](const int i, const int j) {
		mlog.log(5,"  Computing block [",i,",",j,"] ...");
		M.ptr(i,j)=new BLOCK((*this)(i,0).nrows,X(0,j).ncols,Zero());
		for(int k=0; k<ntowers; k++)
		M(i,j)+=(*this)(i,k)*X(k,j);
		mlog.log(5,"  Block [",i,",",j,"] done.");});
	mlog.skip(5);
	return M;
}

// ---- Scalar Valued Operations ---------------------------------------------------------------------------------

template<class BLOCK>
template<class BLOCK2>
FIELD BlockedMatrix<BLOCK>::diff2(const BlockedMatrix<BLOCK2>& X) const {
	assert(nstreets == X.nstreets);
	assert(ntowers == X.ntowers);
	vector<FIELD> sub(nstreets * ntowers);
	MultiLoop(nstreets * ntowers,
			[this,&sub,&X](const int i) {sub[i]=block[i]->diff2(*X.block[i]);});
	FIELD T = 0;
	for (FIELD t : sub)
		T += t;
	return T;
}

// ---- Solvers --------------------------------------------------------------------------------------------------

template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::inverse() const {
	if (nstreets * ntowers == 1)
		return static_cast<BLOCK>(static_cast<Cmatrix>(*block[0]).inverse());
	return BlockedMatrix<BLOCK>(
			static_cast<BLOCK>(convert<Cmatrix>().inverse()), rowstructure(),
			colstructure());
}

template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::inverseSqrt() const {
	if (nstreets * ntowers == 1)
		return static_cast<BLOCK>(static_cast<Cmatrix>(*block[0]).inverseSqrt());
	return BlockedMatrix<BLOCK>(
			static_cast<BLOCK>(convert<Cmatrix>().inverseSqrt()),
			rowstructure(), colstructure());
}

// ---- Reblocking -----------------------------------------------------------------------------------------------

template<class BLOCK>
template<class BLOCK2>
BlockedMatrix<BLOCK>::BlockedMatrix(const BlockedMatrix<BLOCK2>& X,
		const BlockedRemap& rmap, const BlockedRemap& cmap, const bool inverse,
		const bool push) :
		BlockedMatrix((1 - inverse) * rmap.ndest + inverse * rmap.nsource,
				(1 - inverse) * cmap.ndest + inverse * cmap.nsource) {

	vector<Tower<BLOCK>*> newtower(X.ntowers);
	MultiLoop(X.ntowers, [&X,&rmap,&inverse,&push,&newtower](int i) {
		mlog.log(4,"  Remapping tower ",i," ...");
		const Tower<BLOCK2> oldtower=X.virtual_tower(i);
		(newtower)[i] =new Tower<BLOCK>(oldtower,rmap,inverse,push);
		mlog.log(4,"  Tower ",i," done.");});

	Street<BLOCK>* oldstreet[nstreets];
	for (int u = 0; u < nstreets; u++)
		oldstreet[u] = new Street<BLOCK>(X.ntowers,
				newtower[0]->block[u]->nrows);

	for (int u = 0; u < nstreets; u++)
		for (int i = 0; i < X.ntowers; i++) {
			oldstreet[u]->block[i] = newtower[i]->block[u];
			newtower[i]->block[u] = NULL;
		}
	for (int u = 0; u < X.ntowers; u++)
		delete newtower[u];

	Street<BLOCK>* newstreet[nstreets];
	MultiLoop2<Street<BLOCK>**, Street<BLOCK>**>(nstreets,
			[&cmap,&inverse,&push](const int i, Street<BLOCK>** oldstreetp, Street<BLOCK>** newstreetp) {
				mlog.log(4,"  Remapping street ",i," ...");
				(newstreetp)[i] =new Street<BLOCK>(*oldstreetp[i],cmap,inverse,push);
				mlog.log(4,"  Street ",i," done.");}, oldstreet, newstreet);

	for (int u = 0; u < nstreets; u++)
		delete oldstreet[u];

	for (int i = 0; i < nstreets; i++)
		for (int j = 0; j < ntowers; j++) {
			block[j * nstreets + i] = newstreet[i]->block[j];
			newstreet[i]->block[j] = NULL;
		}
	for (int u = 0; u < nstreets; u++)
		delete newstreet[u];

}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const BlockedMatrix& X,
		const class Identity& dummy, const BlockedRemap& cmap,
		const bool inverse) :
		BlockedMatrix(X.nstreets,
				(1 - inverse) * cmap.ndest + inverse * cmap.nsource) {

	Street<BLOCK>* newstreet[nstreets];
	MultiLoop1<Street<BLOCK>**>(nstreets,
			[&X,&cmap,&inverse](const int i, Street<BLOCK>** newstreetp) {
				//const VirtualStreet<BLOCK> oldstreet(X.street(i));
				const Street<BLOCK> oldstreet=X.virtual_street(i);
				(newstreetp)[i] =new Street<BLOCK>(oldstreet,cmap,inverse);},
			newstreet);

	for (int i = 0; i < nstreets; i++)
		for (int j = 0; j < ntowers; j++) {
			block[j * nstreets + i] = newstreet[i]->block[j];
			newstreet[i]->block[j] = NULL;
		}

	for (int u = 0; u < nstreets; u++)
		delete newstreet[u];
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(const BlockedMatrix& X,
		const BlockedRemap& rmap, const class Identity& dummy,
		const bool inverse) :
		BlockedMatrix((1 - inverse) * rmap.ndest + inverse * rmap.nsource,
				X.ntowers) {

	Tower<BLOCK>* newtower[ntowers];
	MultiLoop1<Tower<BLOCK>**>(ntowers,
			[&X,&rmap,&inverse](const int i, Tower<BLOCK>** newtowerp) {
				const Tower<BLOCK> oldtower=X.virtual_tower(i);
				(newtowerp)[i] =new Tower<BLOCK>(oldtower,rmap,inverse);},
			newtower);

	for (int i = 0; i < nstreets; i++)
		for (int j = 0; j < ntowers; j++) {
			block[j * nstreets + i] = newtower[j]->block[i];
			newtower[j]->block[i] = NULL;
		}

	for (int u = 0; u < ntowers; u++)
		delete newtower[u];
}

template<class BLOCK>
BlockedMatrix<BLOCK> BlockedMatrix<BLOCK>::pullColumns(const BlockedRemap& cmap,
		const bool inverse) const {
	int nnewtowers = (!inverse) * cmap.ndest + (inverse) * cmap.nsource;
	assert((!inverse) * cmap.nsource + (inverse) * cmap.ndest == ntowers);
	BlockedMatrix result(nstreets, nnewtowers);

	Street<BLOCK>* newstreet[nstreets];
	MultiLoop1<Street<BLOCK>**>(nstreets,
			[this,&cmap,&inverse](const int i, Street<BLOCK>** newstreetp) {
				//const VirtualStreet<BLOCK> oldstreet(street(i));
				const Street<BLOCK> oldstreet=this->virtual_street(i);
				newstreetp[i]=new Street<BLOCK>(oldstreet,cmap,inverse);},
			newstreet);

	for (int i = 0; i < nstreets; i++)
		for (int j = 0; j < nnewtowers; j++) {
			result.block[j * nstreets + i] = newstreet[i]->block[j];
			newstreet[i]->block[j] = NULL;
		}

	for (int u = 0; u < nstreets; u++)
		delete newstreet[u];

	return result;
}

// ---- I/O ------------------------------------------------------------------------------------------------------

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(DenseMatrixFile& file, const int _nstreets,
		const int _ntowers) :
		BlockedMatrix(_nstreets, _ntowers) {

	assert(_nstreets > 0);
	int rows = file.nrows / _nstreets;
	BlockStructure rstruct(nstreets);
	for (int i = 0; i < nstreets - 1; i++)
		rstruct[i] = rows;
	rstruct[nstreets - 1] = file.nrows - (nstreets - 1) * rows;

	assert(_ntowers > 0);
	int cols = file.ncols / _ntowers;
	BlockStructure cstruct(ntowers);
	for (int i = 0; i < ntowers - 1; i++)
		cstruct[i] = cols;
	cstruct[ntowers - 1] = file.ncols - (ntowers - 1) * cols;

	BlockedMatrix<BLOCK> B(file, rstruct, cstruct);

	B.nodestruct = true;
	for (int I = 0; I < nstreets; I++)
		for (int J = 0; J < ntowers; J++)
			ptr(I, J) = B.ptr(I, J);
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(SparseMatrixFile& file, const int _nstreets,
		const int _ntowers) :
		BlockedMatrix(_nstreets, _ntowers) {

	assert(_nstreets > 0);
	int rows = file.nrows / _nstreets;
	BlockStructure rstruct(nstreets);
	for (int i = 0; i < nstreets - 1; i++)
		rstruct[i] = rows;
	rstruct[nstreets - 1] = file.nrows - (nstreets - 1) * rows;

	assert(_ntowers > 0);
	int cols = file.ncols / _ntowers;
	BlockStructure cstruct(ntowers);
	for (int i = 0; i < ntowers - 1; i++)
		cstruct[i] = cols;
	cstruct[ntowers - 1] = file.ncols - (ntowers - 1) * cols;

	BlockedMatrix<BLOCK> B(file, rstruct, cstruct);

	B.nodestruct = true;
	for (int I = 0; I < nstreets; I++)
		for (int J = 0; J < ntowers; J++)
			ptr(I, J) = B.ptr(I, J);
}

template<class BLOCK>
string BlockedMatrix<BLOCK>::str(const Dense dummy) const {
	ostringstream stream;
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	for (int I = 0; I < nstreets; I++) {
		for (int i = 0; i < block[I]->nrows; i++) {
			stream << "[ ";
			for (int J = 0; J < ntowers; J++) {
				for (int j = 0; j < block[J * nstreets]->ncols; j++)
					if (isFilled(I, i, J, j)) {
						stream.width(6);
						stream << this->read(I, i, J, j) << " ";
					} else
						stream << "      ";
				if (J < ntowers - 1)
					stream << " | ";
			}
			stream << " ]" << endl;
		}
		if (I < nstreets - 1) {
			stream << "[ ";
			for (int J = 0; J < ntowers; J++) {
				for (int j = 0; j < block[J * nstreets]->ncols; j++)
					stream << "-------";
				if (J < ntowers - 1)
					stream << "-+-";
			}
			stream << "-]" << endl;
		}
	}
	return stream.str();
}

template<class BLOCK>
string BlockedMatrix<BLOCK>::str(const Sparse dummy) const {
	ostringstream stream;
	for (int I = 0; I < ntowers; I++)
		for (int J = 0; J < nstreets; J++) {
			stream << "<" << I << "," << J << ">" << endl;
			stream << block[J * nstreets + I]->str(Sparse());
		}
	return stream.str();
}

template<class BLOCK>
void BlockedMatrix<BLOCK>::serialize(Rstream& rstream) const {
	rstream << "BlockedMatrix{" << Rstream::endl;
	rstream.var("nstreets", nstreets);
	rstream.var("ntowers", ntowers);
	for (int J = 0; J < ntowers; J++)
		for (int I = 0; I < nstreets; I++) {
			rstream << "  block[" << I << "," << J << "]=";
			rstream.write(*block[J * nstreets + I]);
		}
	rstream << "}" << Rstream::endl;
}

template<class BLOCK>
void BlockedMatrix<BLOCK>::serialize(Bofstream& ofs) const {
	ofs.tag("BlockedMatrix", 0);
	ofs.write(nstreets);
	ofs.write(ntowers);
	ofs.seriate_array(block, nstreets * ntowers);
}

template<class BLOCK>
BlockedMatrix<BLOCK>::BlockedMatrix(Bifstream& ifs) {
	ifs.check("BlockedMatrix", 0);
	ifs.read(nstreets);
	ifs.read(ntowers);
	ifs.unseriate_array(block);
}

#endif
