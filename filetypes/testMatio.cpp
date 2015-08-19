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



#include <stdlib.h>
#include <stdio.h>
#include <matio.h>
#include <iostream>
#include "MMFglobal.hpp"

using namespace std;

int main(int argc, char **argv) {
	mat_t *matfp;
	matvar_t *matvar;

	matfp = Mat_Open(argv[1], MAT_ACC_RDONLY);

	if (NULL == matfp) {
		fprintf(stderr, "Error opening MAT file \"%s\"!\n", argv[1]);
		return EXIT_FAILURE;
	}

	while ((matvar = Mat_VarReadNext(matfp)) != NULL) {
		printf("%s\n", matvar->name);
		cout << matvar->rank << endl;
		cout << matvar->dims[0] << "," << matvar->dims[1] << endl;
		for (int i = 0; i < 10; i++)
			cout << reinterpret_cast<double*>(matvar->data)[i] << endl;
		Mat_VarFree(matvar);
		matvar = NULL;
	}
	Mat_Close(matfp);
	return EXIT_SUCCESS;
}
