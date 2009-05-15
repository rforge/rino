/*################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <math.h>
#include <stdbool.h>

#include <string>
#include <vector>
using namespace std;
// start extern "C"
extern "C" {
	SEXP listElt(SEXP, char*);
	SEXP makeList(SEXP *, char**);
	void setdims(SEXP, int, int*);
	int *getdims(SEXP obj);
	void setclass(SEXP, char*);
	SEXP C2Rint(int*);
	SEXP C2Rdouble(double*);
	void rm2cm_double(double*, int, int);
	void cm2rm_double(double*, int, int);
}
// end extern "C"
