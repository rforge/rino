#ifndef SOLFUN_H
#define SOLFUN_H

#include <math.h>
#include <stdbool.h>
#include <limits.h>

#include <string>
#include <vector>
using namespace std;

#include "newmatap.h"
#include "newmat.h"
using namespace NEWMAT;

extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
 // R object creation/manipulation functions
  void setdims(SEXP, int, int*);
  int *getdims(SEXP);
  void setclass(SEXP, char*);
  SEXP makeList(SEXP *, char**);
  SEXP listElt(SEXP, char*);

  // R to C object translation
  ReturnMatrix R2CRV(SEXP);
  ReturnMatrix R2CCV(SEXP);
  ReturnMatrix R2Cmat(SEXP, int, int);
  ReturnMatrix R2Cmat2(SEXP);
  ReturnMatrix R2Carr(SEXP);

  // C to R object translation
  SEXP C2Rmat(const Matrix&);
  SEXP C2R3D(const Matrix&, int *);
  SEXP C2R3D2(const Matrix&, int *);
  SEXP C2Rint(int*);
  SEXP C2Rdouble(double *);

  // C <==> FORTRAN object translation
  ReturnMatrix F2C(double *, int, int);
  double *C2F(const Matrix &);

  // Fast, in-place row-major/col-major translation
  void rm2cm_double(double *, int, int);
  void cm2rm_double(double *, int, int);
  void rm2cm_Matrix(Matrix &mat);

  // Simple math util functions for class Matrix
  ReturnMatrix absmat(const Matrix&);
  ReturnMatrix sqrtVec(const Matrix&);
  ReturnMatrix cumprod(const Matrix&);
  ReturnMatrix colsums(const Matrix&);

  // Extra Matrix Subset Functions
  ReturnMatrix diag(const Matrix&);
  ReturnMatrix sqrtdiagmat(const Matrix&);
  ReturnMatrix diagmat(const Matrix&);
  ReturnMatrix divmat(const Matrix&);
  ReturnMatrix dotmatmult(const Matrix&, const Matrix&);
  ReturnMatrix dotcvecmult(const ColumnVector&, const ColumnVector&);
  ReturnMatrix div2mat(const Matrix&, const Matrix&);
  ReturnMatrix div2cvec(const ColumnVector& , const ColumnVector&);
  // Matrix Print Function Overloads
  void printMatrix(const Matrix&);
  void printRVector(const RowVector&);
  void printCVector(const ColumnVector&);

  // Simple math util functions for numerical stability
  Real sqrtHyp(const Real&, const Real&);

  // Functions to draw from different distributions
  ReturnMatrix rnorms(int);
  ReturnMatrix rnorms_mat(int, int);
  SEXP CV2RV(const ColumnVector&);
  SEXP C2RD(const double);
  ReturnMatrix makediag(const int , const double );
  ReturnMatrix vec2diag(const ColumnVector& );
  ReturnMatrix mcbind(const Matrix& , const Matrix& );
  ReturnMatrix mrbind(const Matrix& , const Matrix& );
  ReturnMatrix row2msort(const Matrix& mat);
  ReturnMatrix solveqr(const Matrix& , const Matrix& );
  ReturnMatrix add1vec(const ColumnVector& , const double );
  double solvecond(const Matrix& );
  SEXP solnp_chol(SEXP);
} // end extern "C"
#endif
