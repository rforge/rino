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
/*Collection of helper functions for use in interfacing with R*/

#include "rsolnpfun.h"

extern "C" {
	
  SEXP listElt(SEXP list, char* str)
  {
    R_len_t i;
    SEXP elmt=R_NilValue, names=getAttrib(list, R_NamesSymbol);
    for (i=0; i<length(list); i++){
      if(!(strcmp(CHAR(STRING_ELT(names,i)), str))){
	elmt=VECTOR_ELT(list,i);
	break;
      }
    }
    return elmt;
  }

   SEXP makeList(SEXP *objs, char **names)
  {
    int len = sizeof(objs)/sizeof(SEXP), i;
    SEXP Rlist, Rnames;
    PROTECT(Rlist = allocVector(VECSXP, len));
    PROTECT(Rnames = allocVector(STRSXP, len));
    for(i=0; i<len; i++){
      SET_VECTOR_ELT(Rlist, i, objs[i]);
      SET_STRING_ELT(Rnames, i, mkChar(names[i]));
    }
    setAttrib(Rlist, R_NamesSymbol, Rnames);
    UNPROTECT(1);
    return Rlist;
  }

   void setdims(SEXP obj, int num_dims, int* dims){
    SEXP dim;
    PROTECT(dim = allocVector(INTSXP, num_dims));
    int *idim=INTEGER(dim);
    for(int i=0; i<num_dims; i++) idim[i] = dims[i];
    setAttrib(obj, R_DimSymbol, dim);
    UNPROTECT(1);
  }

  int *getdims(SEXP obj)
  { return INTEGER(getAttrib(obj,R_DimSymbol)); }

  void setclass(SEXP obj, char *name){
     SEXP cls;
     PROTECT(cls = allocVector(STRSXP, 1));
     SET_STRING_ELT(cls, 0, mkChar(name));
     classgets(obj, cls);
     UNPROTECT(1);
  }

    SEXP C2Rint(int *x)
  {
    SEXP R; int i, n=sizeof(x)/sizeof(int);
    PROTECT(R = allocVector(INTSXP, n));
    int *r=INTEGER(R); for(i=0;i<n;i++) r[i]=x[i];
    UNPROTECT(1); return R;
  }

  SEXP C2Rdouble(double *x)
  {
    SEXP R; int i, n=sizeof(x)/sizeof(double);
    PROTECT(R=allocVector(REALSXP, n));
    double *r=REAL(R); for(i=0;i<n;i++) r[i]=x[i];
    UNPROTECT(1); return R;
  }

  void rm2cm_double(double *m, int nr, int nc)
  { for(int i=0;i<nc;i++) for(int j=i+1;j<nr;j++) swap(m[i*nr+j],m[j*nc+i]); }

  void cm2rm_double(double *m, int nr, int nc)
  { for(int i=0;i<nr;i++) for(int j=i+1;j<nc;j++) swap(m[i*nr+j],m[j*nc+i]); }
}