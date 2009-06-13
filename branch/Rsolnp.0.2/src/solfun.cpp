#include "solfun.h"

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

  ReturnMatrix R2Carr(SEXP Robj)
  {
    int i, ndims=length(getAttrib(Robj, R_DimSymbol));
    int *dims=INTEGER(getAttrib(Robj, R_DimSymbol));
    int len=1; for(i=0;i<ndims;i++) len*=dims[i];
    ColumnVector vals(len); Real *Rptr=REAL(Robj); vals << Rptr;
    //    for(i=1;i<=len;i++) vals(i)= Rptr[i-1];

    if(ndims==3){
      int x=dims[0], y=dims[1], z=dims[2];
      Matrix out(z,x*y);
      for(int i=0; i<z; i++)
	out.Row(i) = vals.Rows(1+x*y*i,x*y*(i+1));
      out.Release(); return out.ForReturn();
    }
    else if(ndims==2){
      vals = vals.AsMatrix(dims[0], dims[1]);
      vals.Release(); return vals.ForReturn();
    }
    else if(ndims==1){
      vals.Release(); return vals.ForReturn();
    }
    else{
      Rprintf("Invalid number of dims (%d)... \nvalid range: 1<=ndim<=3)", ndims);
      Matrix out;
      out.Release(); return out.ForReturn();
    }
  }

  ReturnMatrix R2Cmat(SEXP mat, int nr, int nc)
  {
    // Catch...
    int i, j, *dims=INTEGER(getAttrib(mat, R_DimSymbol));
    int rw=dims[0], cl=dims[1];
    Matrix X(rw,cl);
    Real *rmat = REAL(mat), *x=X.Store();
    for (i=0; i<cl; i++)
      for (j=0; j<rw; j++) x[j*cl+i] = rmat[i*rw+j];
    if (rw*cl!=nr*nc) X = X.SubMatrix(1,nr,1,nc);

    // ...and release.
    X.Release(); return X.ForReturn();
  }

  ReturnMatrix R2Cmat2(SEXP matR)
  {
    int i, j, *d=getdims(matR), rw=d[0], cl=d[1];
    Matrix X(rw,cl); double *mat=REAL(matR), *x=X.Store();
    for(i=0;i<cl;i++) for(j=0;j<rw;j++) x[j*cl+i]=mat[i*rw+j];
    X.Release(); return X.ForReturn();
  }

  ReturnMatrix R2CCV(SEXP vec)
  {
    PROTECT(vec = coerceVector(vec,REALSXP));
    int len = length(vec);
    ColumnVector output(len); output=0.0;
    Real *Rptr = REAL(vec), *Cptr = output.Store();
    for(int i=0; i<len; i++) Cptr[i] = Rptr[i];
    UNPROTECT(1);
    output.Release(); return output.ForReturn();
  }

  ReturnMatrix R2CRV(SEXP vec)
  {
    PROTECT(vec = coerceVector(vec,REALSXP));
    int len = length(vec);
    RowVector output(len); output=0.0;
    Real *Rptr = REAL(vec), *Cptr = output.Store();
    for(int i=0; i<len; i++) Cptr[i] = Rptr[i];
    UNPROTECT(1);
    output.Release(); return output.ForReturn();
  }

  SEXP C2Rmat(const Matrix& Cmat)
  {
    SEXP Rmat;
    int i, j, nr=Cmat.Nrows(), nc=Cmat.Ncols();
    PROTECT(Rmat = allocVector(REALSXP, nr*nc));
    Real *Rptr = REAL(Rmat);
    for(i=0; i<nc; i++)
      for(j=0; j<nr; j++)
	Rptr[i*nr+j] = Cmat.Store()[j*nc+i];
    int dims[]={nr,nc};  setdims(Rmat,2,dims);
    UNPROTECT(1);
    return Rmat;
  }

  SEXP C2R3D(const Matrix& obj, int *dims)
  {
    SEXP Rarr;
    int i, j, k, x=dims[0], y=dims[1], z=dims[2];
    PROTECT(Rarr = allocVector(REALSXP, x*y*z));
    Real *Rptr = REAL(Rarr);
    Matrix tmp(x,y);
    for(k=0; k<z; k++){
      tmp = obj.SubMatrix(1,x,1+y*k,y*(k+1));
      for(i=0; i<y; i++)
	for(j=0; j<x; j++)
	  Rptr[i*x+j+(k*x*y)] = tmp.Store()[j*y+i];
    }
    setdims(Rarr, 3, dims);
    UNPROTECT(1);
    return Rarr;
  }

  SEXP C2R3D2(const Matrix& obj, int *dims)
  {
    Matrix tmp = obj.t();
    int i, x=dims[0], y=dims[1], z=dims[2];
    SEXP output;

    PROTECT(output=allocVector(REALSXP,x*y*z));
    Real *Rptr=REAL(output), *Cptr=tmp.Store();
    for(i=1; i<=tmp.Storage(); i++) Rptr[i]=Cptr[i];
    setdims(output, 3, dims);
    UNPROTECT(1);

    return output;
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

  // Convert NEWMAT Matrix to a column-major double* array
  // representation for use with Fortran subroutines
  //
  // NOTE: Variables declared using this function must be deleted
  //       from memory after final use using 'Free(obj_name);'

  double* C2F(const Matrix &C)
  {
    int i, j, nr=C.Nrows(), nc=C.Ncols();
    double *c=C.Store(), *f=(double*)Calloc(nr*nc,double);
    for(i=0;i<nc;i++) for(j=0;j<nr;j++) f[i*nr+j]=c[j*nc+i];
    return f;
  }

  ReturnMatrix F2C(double *f, int nr, int nc)
  {
    Matrix X(nr,nc); double *x=X.Store(); int i, j;
    for (i=0; i<nc; i++)
      for (j=0; j<nr; j++)
	x[j*nc+i] = f[i*nr+j];
    X.Release(); return X.ForReturn();
  }

  void printMatrix(const Matrix& mat)
  {
    int i,j;
    Rprintf("\n");
    for(i=1; i<=mat.Nrows(); i++){
      for(j=1; j<=mat.Ncols(); j++){
	if(!((j-1)%6)) Rprintf("\n[%d,%d]\t",i,j);
	Rprintf("%f\t", mat(i,j));
      }
    }
    Rprintf("\n\n");
  }

  void printRVector(const RowVector& vec)
  {
    int i;
    for(i=1;i<=vec.Storage();i++){
      if(!(i%7) || i==1) Rprintf("\n[%d] ",i);
      Rprintf(" %f ",vec(i));
    }
    Rprintf("\n\n");
  }

  void printCVector(const ColumnVector& vec)
  {
    Rprintf("\n");
    for(int i=1; i<=vec.Storage(); i++)
      Rprintf("[%d,] %f\n",i,vec(i));
  }

  ReturnMatrix absmat(const Matrix& mat)
  {
    Matrix out=mat; Real *ptr=out.Store();
    for(int i=0; i<out.Storage(); i++) if(ptr[i]<0) ptr[i]=-ptr[i];
    out.Release(); return out.ForReturn();
  }

  ReturnMatrix diag(const Matrix& mat)
  {
    int i, m=mat.Nrows(); RowVector D(m); double *d=D.Store();
    for(i=0; i<m; i++) d[i]=mat(i+1,i+1); D.Release(); return D.ForReturn();
  }

  ReturnMatrix sqrtdiagmat(const Matrix& mat)
  {
    int i, m=mat.Nrows(); Matrix D(m,1);
    for(i=0; i<m; i++){
    	if(mat(i+1,i+1)>=0){
    		D(i+1,1)=sqrt(mat(i+1,i+1));
    	} else{
    		D(i+1,1)=0;
    	}
    }
    	D.Release(); return D.ForReturn();
  }

  ReturnMatrix diagmat(const Matrix& mat)
  {
	  int i, m=mat.Nrows();
	  Matrix out(m,m);
	  out=0.0;
	  for(i=0; i<m; i++) out(i+1,i+1)=mat(i+1,i+1);
	  out.Release(); return out.ForReturn();
  }

  ReturnMatrix dotmatmult(const Matrix& mat1, const Matrix& mat2)
  {
	  int i, j, m=mat1.Nrows(), n=mat1.Ncols();
	  Matrix out(m,n);
	  out=0.0;
	  for(i=0; i<m; i++) for(j=0;j<n;j++) out(i+1,j+1)=mat1.SubMatrix(i+1,i+1,j+1,j+1).AsScalar()*mat2.SubMatrix(i+1,i+1,j+1,j+1).AsScalar();
	  out.Release(); return out.ForReturn();
  }

  ReturnMatrix dotcvecmult(const ColumnVector& vec1, const ColumnVector& vec2)
  {
	  int i, n1=vec1.Nrows(), n2=vec2.Nrows();
	  if(n1!=n2){
		  Rprintf("Invalid dims of vectors in dotcvecmult routine");
		  ColumnVector out;
		  out.Release(); return out.ForReturn();
	  }
	  ColumnVector out(n1);
	  out=0.0;
	  for(i=0; i<n1; i++) out.Row(i+1)=vec1.Row(i+1).AsScalar()*vec2.Row(i+1).AsScalar();
	  out.Release(); return out.ForReturn();
  }

  ReturnMatrix divmat(const Matrix& mat)
  {
	  int i,j, m=mat.Nrows(), n=mat.Ncols();
	  Matrix out(m,n);
	  double tmp;
	  out=0.0;
	  for(i=0; i<m; i++){
		  for(j=0;j<n;j++){
			  tmp=mat.SubMatrix(i+1,i+1,j+1,j+1).AsScalar();
			  if(tmp!=0.0){
				  out(i+1,j+1)=1.0/tmp;
			  } else{
				  out(i+1,j+1)=1.0/2.220446e-16;
			  }
		  }
	  }
	  out.Release();
	  return out.ForReturn();
  }

  ReturnMatrix div2mat(const Matrix& mat1, const Matrix& mat2)
  {
	  /* assumes equal dimensions*/

	  int i,j, m1=mat1.Nrows(), n1=mat1.Ncols();
	  int m2=mat2.Nrows(), n2=mat2.Ncols();
	  if(m1!=m2 || n1!=n2){
	      Rprintf("Invalid dims of matrices in div2mat routine");
	      Matrix out;
	      out.Release(); return out.ForReturn();
	  }
	  Matrix out(m1,n1);
	  out=0.0;
	  double tmp;
	  for(i=0; i<m1; i++){
		  for(j=0;j<n1;j++){
			  tmp = mat2.SubMatrix(i+1,i+1,j+1,j+1).AsScalar();
			  if(tmp!=0){
				  out.SubMatrix(i+1,i+1,j+1,j+1)=mat1.SubMatrix(i+1,i+1,j+1,j+1).AsScalar()/mat2.SubMatrix(i+1,i+1,j+1,j+1).AsScalar();
			  } else{
				  out.SubMatrix(i+1,i+1,j+1,j+1)=1/2.220446e-16;
			  }
		  }
	  }
	  out.Release();
	  return out.ForReturn();
  }

  ReturnMatrix div2cvec(const ColumnVector& vec1, const ColumnVector& vec2)
  {
	  /* assumes equal dimensions*/

	  int i, n1=vec1.Nrows(), n2=vec2.Nrows();
	  if(n1!=n2){
		  Rprintf("Invalid dims of vectors in div2cvec routine");
		  ColumnVector out;
		  out.Release(); return out.ForReturn();
		  }
	  ColumnVector out(n1);
	  out=0.0;
	  double tmp;
	  for(i=0; i<n1; i++){
			  tmp = vec2.Row(i+1).AsScalar();
			  if(tmp!=0.0){
				  out.Row(i+1)=vec1.Row(i+1).AsScalar()/vec2.Row(i+1).AsScalar();
			  } else{
				  out.Row(i+1)=1/2.220446e-16;
			  }
		  }
	  out.Release();
	  return out.ForReturn();
  }

  ReturnMatrix sqrtVec(const Matrix& mat)
  {
    RowVector vec = mat.AsRow(); Real *ptr = vec.Store();
    for(int i=0;i<mat.Storage();i++) ptr[i] = sqrt(ptr[i]);
    vec.Release(); return vec.ForReturn();
  }

  Real sqrtHyp(const Real& x, const Real& y)
  {
    Real z;
    if(fabs(x) > fabs(y)){z = fabs(x)*sqrt(1+(y/x)*(y/x));}
    else if (y!=0){z = fabs(y)*sqrt(1+(x/y)*(x/y));} else {z = 0.0;}
    return z;
  }

  void rm2cm_double(double *m, int nr, int nc)
  { for(int i=0;i<nc;i++) for(int j=i+1;j<nr;j++) swap(m[i*nr+j],m[j*nc+i]); }

  void cm2rm_double(double *m, int nr, int nc)
  { for(int i=0;i<nr;i++) for(int j=i+1;j<nc;j++) swap(m[i*nr+j],m[j*nc+i]); }

  void rm2cm_Matrix(Matrix &mat)
  {
    double *m=mat.Store(); int i, j, nr=mat.Nrows(), nc=mat.Ncols();
    for(i=1;i<=nc;i++) for(j=i+1;j<=nr;j++) swap(m[i*nr+j],m[j*nc+i]);
  }

  ReturnMatrix cumprod(const Matrix& mat)
  {
    int i, n=mat.Storage(); RowVector tmp(n); tmp << mat.Store();
    for(i=2;i<=n;i++) tmp(i)=i*tmp(i-1);
    tmp.Release(); return tmp.ForReturn();
  }

  ReturnMatrix colsums(const Matrix& mat)
  {
    int i, j, m=mat.Nrows(), n=mat.Ncols();
    RowVector out(n); double *csums=out.Store();
    for(i=0;i<n;i++){  csums[i]=0.0;
      for(j=0;j<m;j++) csums[i]+=mat(j,i);
    }
    out.Release(); return out.ForReturn();
  }

  SEXP CV2RV(const ColumnVector& vec)
  {
	SEXP RV;
	int i, nr=vec.Nrows();
	PROTECT(RV = allocVector(REALSXP, nr));
	double *r=REAL(RV);
	for(i=0; i<nr; i++) r[i]=vec.Row(i+1).AsScalar();
	UNPROTECT(1);
	return RV;
  }

  SEXP C2RD(const double x){
	  SEXP RV;
	  PROTECT(RV = allocVector(REALSXP, 1));
	  double *r=REAL(RV);
	  r[0]=x;
	  UNPROTECT(1);
	  return RV;
  }
  ReturnMatrix makediag(const int m, const double d)
  {
    Matrix out(m,m);
    out=0.0;
    for(int i=0; i<m; i++) out(i+1,i+1)=d;
    out.Release(); return out.ForReturn();
  }

  ReturnMatrix vec2diag(const ColumnVector& vec)
  {
	int m=vec.Nrows();
    Matrix out(m,m);
    out=0.0;
    for(int i=0; i<m; i++) out.SubMatrix(i+1,i+1,i+1,i+1)=vec.Row(i+1).AsScalar();
    out.Release(); return out.ForReturn();
  }

  ReturnMatrix mcbind(const Matrix& mat1, const Matrix& mat2)
  {
	  int m1 = mat1.Nrows();
	  int m2 = mat2.Nrows();
	  int n1 = mat1.Ncols();
	  int n2 = mat2.Ncols();
	  if(m1!=m2){
		  Rprintf("Invalid dims of rows in mcbind routine");
		  Matrix out;
		  out.Release(); return out.ForReturn();
	  	  }
	  Matrix out(m1, n1+n2);
	  out.SubMatrix(1,m1,1,n1)=mat1;
	  out.SubMatrix(1,m1,n1+1,n1+n2)=mat2;
	  out.Release(); return out.ForReturn();
  }

  ReturnMatrix mrbind(const Matrix& mat1, const Matrix& mat2)
  {
	  int m1 = mat1.Nrows();
	  int m2 = mat2.Nrows();
	  int n1 = mat1.Ncols();
	  int n2 = mat2.Ncols();
	  if(n1!=n2){
		  Rprintf("Invalid dims of cols in rcbind routine");
		  Matrix out;
		  out.Release(); return out.ForReturn();
		  }
	  Matrix out(m1+m2, n1);
	  out.SubMatrix(1,m1,1,n1)=mat1;
	  out.SubMatrix(m1+1,m1+m2,1,n1)=mat2;
	  out.Release(); return out.ForReturn();
  }

  ReturnMatrix row2msort(const Matrix& mat)
  {
	  int m1 = mat.Nrows();
	  int n1 = mat.Ncols();
	  double r1, r2;
	  if(n1!=2){
		  Rprintf("Invalid dims of cols in row2msort routine");
		  Matrix out;
		  out.Release(); return out.ForReturn();
	  }
	  Matrix out(m1,n1);
	  out=mat;
	  for(int i=1;i<=m1;i++){
		  r1=mat.SubMatrix(i,i,1,1).AsScalar();
		  r2=mat.SubMatrix(i,i,2,2).AsScalar();
		  if(r1>r2){
			  out.SubMatrix(i,i,1,1)=r2;
			  out.SubMatrix(i,i,2,2)=r1;

		  }
	  }
	  out.Release(); return out.ForReturn();
  }

  ReturnMatrix solveqr(const Matrix& mat1, const Matrix& mat2)
  {
	  //need some additional try/catch for singular matrix
	  // in qr procedure.
	  Matrix mm1=mat1;
	  Matrix mm2=mat2;
	  UpperTriangularMatrix U;
	  Matrix M;
	  try{
	  QRZ(mm1, U);
	  QRZ(mm1, mm2, M);
	  }
	  catch(Exception){
		  Rprintf(Exception::what());
	  }
	  Matrix out = U.i() * M;
	  U.Release();
	  M.Release();
	  mm1.Release();
	  mm2.Release();
	  out.Release();
	  return out.ForReturn();
  }

  ReturnMatrix add1vec(const ColumnVector& vec, const double d)
  {
	  int m=vec.Nrows();
	  Matrix out(m+1,1);
	  out.SubMatrix(1,m,1,1)=vec;
	  out.SubMatrix(m+1,m+1,1,1)=d;
	  out.Release();
	  return out.ForReturn();
  }

  double solvecond(const Matrix& mat)
  {
	  int j;
	  int m=mat.Nrows();
	  int n=mat.Ncols();
	  int l=fmin(m,n);
	  double x=0.0;
	  Matrix mat2;
	  if(m<n){
		  mat2=mat.t();
	  } else{
		  mat2=mat;
	  }
	  Matrix U, V; DiagonalMatrix D;
	  SVD(mat2, D, U, V);
	  double tmp;
	  for(j=0;j<l;j++){
		  tmp=D.SubMatrix(j+1,j+1,j+1,j+1).AsScalar();
		  if(tmp == 0.0){
			  x=1/2.220446e-16;
			  j=l;
			  }
	  }
	  if(x==0.0){
		  x = D.Maximum()/D.Minimum();
	  }
	  D.Release();
	  U.Release();
	  V.Release();
	  mat2.Release();
	  return x;
  }

  SEXP solnp_chol(SEXP A)
  {
  	SEXP ans = PROTECT((TYPEOF(A) == REALSXP)?duplicate(A):
  			   coerceVector(A, REALSXP));
  	SEXP adims = getAttrib(A, R_DimSymbol);
  	int m = INTEGER(adims)[0];
  	int n = INTEGER(adims)[1];
  	int i, j;
  	for (j = 0; j < n; j++) {	/* zero the lower triangle */
  	    for (i = j+1; i < n; i++) {
  		REAL(ans)[i + j * n] = 0.;
  	    }
  	}
  	F77_CALL(dpotrf)("Upper", &m, REAL(ans), &m, &i);
  	if (i != 0) {
  	    if (i > 0)
  		Rprintf("the leading minor of order %d is not positive definite",i);
  	    Rprintf("argument %d of Lapack routine %s had invalid value",-i, "dpotrf");
  	}
  	unprotect(1);
  	return ans;
  }
}
