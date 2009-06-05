#include "solnp.h"

/* inputs:
 *  * parameters :  the parameter vector
 *  * control : the solver control vector
 *  (rho, minit, delta, tol, neq, nineq, np, lpb, trace)
 *  * info : problem information vector
 *  (tc = total constraints, nineq = ineq constraints,
 *  neq = equality constraints, ...)
 *  * fn : the names of the function to be called
 *  (F=obj, J=ineq, H=eq)
 *  * Constraint = Inequality  & Equality Constraint Vector
 *  * Bnds = Parameter Bounds
 *  * valf = vector of function and constraint values
 *  * hessian = vector for hessian (nxn)
 */
extern "C" {
SEXP solnp(SEXP pars, SEXP Ofn, SEXP Efn, SEXP Ifn, SEXP xeq, SEXP xpb, SEXP xscale,
		SEXP control, SEXP rho, SEXP xlgr, SEXP xmu, SEXP xhess, SEXP xn)
{
	double *n=REAL(xn);
	/* length of scale vector*/
	int scln = Rf_length(xscale);
	int i;

	/* no. of parameters, equality and inequality constraints*/
	int np=(int) n[0];
	int neq=(int) n[1];
	int nineq=(int) n[2];

	/* outpar is the resulting list
	 * hess is the resulting hessian
	 * lagrange is the resulting lagrange multipliers
	 * xmu is the mean*/
	SEXP outpar, outhess, outlagrange, outmu;
	PROTECT(outpar = allocVector(REALSXP, np));
	PROTECT(outhess = allocMatrix(REALSXP, np, np));
	PROTECT(outlagrange = allocVector(REALSXP, np));
	PROTECT(outmu = allocVector(REALSXP, np));
	/* kp is the no. of SEXP we are protecting*/
	int kp = 4;
	/* control vector holding optimization parameters*/
	double *ctrl = REAL(control);
	int xrho   = (int) ctrl[ 0 ];
	int maxit = (int) ctrl[ 1 ];
	double delta = ctrl[ 2 ];
	double tol   = ctrl[ 3 ];
	int *lpb;
	double *alp;
	lpb = (int *) malloc(2);
	lpb[0] = (int) ctrl[4];
	lpb[1] = (int) ctrl[5];

	alp = (double *) malloc(3);
	alp[0] = 0.0; alp[1] = 0.0; alp[2] = 0.0;

	/*change parameter*/
	double ch = 1.0;
	/*nc = no of constraints*/
	int nc = neq + nineq;
	/*npic parameters + inequalities*/
	int npic = np + nineq;

	Matrix hess(npic, npic);
	hess=R2Cmat(xhess,npic, npic);

	ColumnVector scale(scln);
	scale = R2CCV(xscale);

	ColumnVector ob(1 + neq + nineq);
	ob=R2CCV(solnpeval(pars, Ofn, Efn, Ifn, rho, xn));
	// need special operator for columvector division
	ob=div2cvec(ob,scale.Rows(1, nc + 1));

	ColumnVector p0( nineq + np );
	p0=R2CCV(pars);
	// need special operator for columvector division
	p0=div2cvec(p0,scale.Rows(neq + 2, nc + np + 1));

	ColumnVector tempv(np);
	tempv=0.0;
	/*parameter bounds (inequality and parameter)*/

	int mm = npic;
	if( lpb[1] == 1 ) {
		if( lpb[0] == 0 ) {
			mm = nineq;
		}
		else {
			mm = npic;
		}
		Matrix pb(mm,2);
		pb=R2Cmat(xpb,mm,2);
		Matrix cscale(mm,2);
		cscale.SubMatrix(1,mm,1,1) = scale.Rows(neq+2,neq+mm+1);
		cscale.SubMatrix(1,mm,2,2) = scale.Rows(neq+2,neq+mm+1);
		pb = div2mat(pb,cscale);
	}

	/* scale the lagrange multipliers and the Hessian*/
	if( nc > 0 ){
		ColumnVector yy(nc);
		yy = dotcvecmult((1/scale.Row(1).AsScalar())*yy,scale.Rows(2,nc+1));
	}
	Matrix tmp1(npic, npic);
	tmp1 = scale.Rows((neq + 2),(nc + np + 1))*scale.Rows((neq + 2),(nc + np + 1)).t();
	hess = dotmatmult(hess, tmp1);
	hess = hess*(1/scale.Row(1).AsScalar());
	double j = ob.Row(1).AsScalar();

	if( nineq > 0){
		if( neq == 0){
			Matrix a(nineq, nineq+np);
			a=0.0;
			a.SubMatrix(1,nineq,1,nineq) = makediag(nineq,-1.0);
		}
		else{
			Matrix a(neq+nineq, nineq+np);
			a=0.0;
			a.SubMatrix(neq+1,neq+nineq,1,nineq) = makediag(nineq,-1.0);
		}
	}

	if( neq > 0 && nineq == 0 ) {
		Matrix a(neq, np);
		a = 0.0;
	}

	if( neq == 0 && nineq == 0 ) {
		RowVector  a(np);
		a = 0.0;
	}

	ColumnVector g(npic);
	g=0.0;

	if( nc > 0 ) {
		ColumnVector constraint (np);
		constraint = 0.0;
		constraint = ob.Rows( 2,(nc + 1));
		for(i=1;i<=np;i++){
			p0.Row(nineq + i)+= delta;
			tempv=dotcvecmult(p0.Rows(nineq+1,npic),scale.Rows(nc + 2, nc + np + 1));
			ob=R2CRV(solnpeval(CV2RV(tempv), Ofn, Efn, Ifn, rho, xn));
			//continue here ob = ob/scale....;
			//create temp vector to pass to
			// evalfei function
		}
	}


}

SEXP solnpeval(SEXP pars, SEXP Ofn, SEXP Efn, SEXP Ifn, SEXP rho, SEXP n)
{
	if(!isEnvironment(rho)) error("'rho' should be an environment");
	SEXP fx, xpt, outpar;
	int k;
	PROTECT(xpt = findVarInFrame(rho, install(".par")));
	PROTECT(fx = allocVector(REALSXP, 1));
	double *nn=REAL(n);
	int np=(int) nn[0];
	int neq=(int) nn[1];
	int nineq=(int) nn[2];
	PROTECT(outpar = allocVector(REALSXP, 1+neq+nineq));
	k=3;
	ColumnVector ob(1+nineq+neq);
	//ColumnVector newpars(np);
	//newpars=R2CCV(xpt);
	ob=0.0;
	fx=eval(Ofn, rho);
	ob.Row(1)=R2CCV(fx);
	if(neq>0){
		SEXP ex;
		PROTECT(ex = allocVector(REALSXP, neq));
		k+=1;
		ex = eval(Efn, rho);
		ob.Rows(2,neq+1)=R2CCV(ex);
	}
	if(nineq>0){
		SEXP ix;
		PROTECT(ix = allocVector(REALSXP, nineq));
		k+=1;
		ix = eval(Ifn, rho);
		ob.Rows(2+neq,nineq+neq+1)=R2CCV(ix);
	}
	//defineVar(install(".par"), CV2RV(newpars), rho);
	//UNPROTECT(2);
	outpar=CV2RV(ob);
	ob.Release();
	UNPROTECT(k);
	return outpar;
}
}
