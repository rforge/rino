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
#define EPSX 2.220446e-16
SEXP solnp(SEXP pars, SEXP Ofn, SEXP Efn, SEXP Ifn, SEXP xeq, SEXP xpb, SEXP xscale,
		SEXP control, SEXP rho, SEXP xlagrange1, SEXP xlagrange2, SEXP xhess, SEXP xn,
		SEXP obin)
{
	//try{
	int i, minit, mm;
	//int lgn1 = Rf_nrows(xlagrange1);
	SEXP MU, SM;
	PROTECT(MU = allocVector(REALSXP, 1));
	PROTECT(SM = allocVector(REALSXP, 1));
	double lagrange2=REAL(xlagrange2)[0];
	/* no. of parameters, equality and inequality constraints*/
	int np = int(REAL(xn)[0]);
	//Rprintf("np",C2Rint(&np));
	int neq = int(REAL(xn)[1]);
	int nineq = int(REAL(xn)[2]);
	//Rprintf("nineq",C2Rint(&nineq));
	/*nc = no of constraints*/
	int nc = neq + nineq;
	/*npic parameters + inequalities*/
	int npic = np + nineq;
	/* outvec is the resulting list*/
	SEXP outvec;
	PROTECT(outvec = allocVector(VECSXP, 4));
	/* control vector holding optimization parameters*/
	double xrho   = double ( REAL(control)[0] );
	int maxit = int ( REAL(control)[1] );
	double delta = double ( REAL(control)[2] );
	double tol   = double ( REAL(control)[3] );
	ColumnVector lagrange1;
	lagrange1 = R2CCV(xlagrange1);
	/*parameter and inequality bounds signal*/
	int lpb1 = int ( REAL(control)[4] );
	int lpb2 = int ( REAL(control)[5] );
	double alp1= 0.0, alp2 = 0.0, alp3 = 0.0;
	double sc1, sc2;

	/*holder for function return value*/
	ColumnVector sob(3);
	/*change parameter*/
	double ch = 1.0;

	/* hessian matrix*/
	Matrix hess;
	hess=R2Cmat(xhess,npic, npic);

	/* scaling vector*/
	ColumnVector scale;
	scale = R2CCV(xscale);

	/*parameter bounds (inequality and parameter)*/
	Matrix pb;
	/*gap matrix*/
	Matrix gap, gap2, ptt;
	/* vectors to hold function values*/
	ColumnVector ob(1 + neq + nineq);
	ColumnVector obm(1 + neq + nineq);
	ob=R2CCV(obin);
	ob=div2cvec(ob,scale.Rows(1, nc + 1));
	RowVector cx;
	ColumnVector dx, yg, sx, u, ob1, ob2, ob3;
	Matrix cz;
	/*parameter vector after p0*/
	ColumnVector p;

	/* vectors to hold parameter values*/
	ColumnVector p0(nineq + np);
	p0=R2CCV(pars);
	p0=div2cvec(p0,scale.Rows(neq + 2, nc + np + 1));

	/* vectors to hold temporary par values*/
	ColumnVector tempv(np);
	tempv=0.0;

	/*constraint vector*/
	ColumnVector constraint;

	/*temporary vectors and matrices*/
	Matrix mtmp1, mtmp2, mtmp3, mtmp4, mtmp5;
	Matrix a;
	ColumnVector b, y, v, vtmp1, vtmp2, vtmp3, vtmp4;
	SymmetricMatrix S; LowerTriangularMatrix L;

	double j=0.0, obmm=0.0, obn=0.0, z, reduce=0.0, go=0.0;
	double dtmp1=0.0, dtmp2=0.0, dtmp3=0.0, dtmp4=0.0;
	mm = npic;
	if( lpb2 == 1) {
		if( lpb1 == 0 ) {
			mm = nineq;
		}
		else {
			mm = npic;
		}
		pb=R2Cmat(xpb,mm,2);
		Matrix cscale(mm,2);
		cscale.SubMatrix(1,mm,1,1) = scale.Rows(neq+2,neq+mm+1);
		cscale.SubMatrix(1,mm,2,2) = scale.Rows(neq+2,neq+mm+1);
		pb = div2mat(pb,cscale);
		//will not use cscale again so release.
		cscale.Release();
	}
	/* scale the lagrange multipliers and the Hessian*/
	if( nc > 0 ){
		vtmp1 = lagrange1 * (1/scale.Row(1).AsScalar());
		lagrange1 = dotcvecmult(vtmp1,scale.Rows(2,nc+1));
		vtmp1.CleanUp();
	}
	mtmp1 = scale.Rows(neq + 2,nc + np + 1)*scale.Rows(neq + 2,nc + np + 1).t();
	mtmp1 = mtmp1 * (1/scale.Row(1).AsScalar());
	mtmp2 = dotmatmult(hess, mtmp1);
	hess = mtmp2;
	j = ob.Row(1).AsScalar();
	int arow=1;
	mtmp1.CleanUp();
	mtmp2.CleanUp();
	R_CheckUserInterrupt();
	if( nineq > 0){
		if( neq == 0){
			a.ReSize(nineq, npic);
			a=0.0;
			mtmp1.ReSize(nineq,nineq);
			mtmp1=0.0;
			mtmp1=makediag(nineq,-1.0);
			mtmp2.ReSize(nineq,np);
			mtmp2=0.0;
			a=mtmp1|mtmp2;
			mtmp1.CleanUp();
			mtmp2.CleanUp();
			arow=nineq;
		} else{
			a.ReSize(neq+nineq, nineq+np);
			a=0.0;
			arow=neq + nineq;
			a.SubMatrix(neq+1,neq+nineq,1,nineq) = makediag(nineq,-1.0);
		}
	}

	if( neq > 0 && nineq == 0 ) {
		a.ReSize(neq, np);
		a = 0.0;
		arow=neq;
	}

	if( neq == 0 && nineq == 0 ) {
		a.CleanUp();
		RowVector  a(np);
		a = 0.0;
		arow=1;
	}
	/*evaluate gradients using forward differences*/
	/*gradient*/
	ColumnVector g(npic);
	g=0.0;

	if( nc > 0 ) {
		constraint.ReSize(nc);
		constraint = ob.Rows( 2, nc + 1);
		for(i=1;i<=np;i++){
			ob=0.0;
			dtmp2=p0.Row(nineq + i).AsScalar();
			p0.Row(nineq + i)=dtmp2+delta;
			tempv=dotcvecmult(p0.Rows(nineq+1,npic),scale.Rows(nc + 2, nc + np + 1));
			ob = solnpeval(CV2RV(tempv), Ofn, Efn, Ifn, rho, xn);
			ob = div2cvec(ob,scale.Rows(1,nc+1));
			g.Row( nineq + i)=(ob.Row(1).AsScalar() - j)/delta;
			vtmp1 = ob.Rows(2,nc+1)-constraint;
			a.SubMatrix(1,arow,nineq+i,nineq+i)=vtmp1/delta;
			p0.Row(nineq+i)= dtmp2;
			tempv.CleanUp();
			vtmp1.CleanUp();
		}
		// check
		R_CheckUserInterrupt();
		if( nineq > 0 ) {
			vtmp1=constraint.Rows( neq + 1, neq + nineq) - p0.Rows( 1, nineq );
			constraint.Rows( neq + 1, neq + nineq ) = vtmp1;
			vtmp1.CleanUp();
		}
		dtmp1 = solvecond(a);
		if( dtmp1 > (1.0/EPSX) ) {
			Rprintf("SOLNP-->Redundant constraints were found. Poor\n intermediate "
					"results may result.Suggest that you\n remove redundant constraints "
					"and re-OPTIMIZE\n");
		}
		dtmp1=0.0;
		b.ReSize(nc);
		b=0.0;
		b = (a * p0) - constraint;
	}
	if( nc > 0){
		ch = -1.0;
		alp1 = tol - constraint.MaximumAbsoluteValue();

		if( alp1 <= 0){
			ch = 1.0;

			if( lpb2 == 0 ) {
				mtmp2 = a*a.t();
				vtmp2.ReSize(npic);
				vtmp2=0.0;
				vtmp2 = p0 - ( a.t() * (mtmp2.i()*constraint) );
				p0 = vtmp2;
				mtmp2.CleanUp();
				vtmp2.CleanUp();
				alp1=1;
			}
		}

		if(alp1 <= 0 ){
			vtmp1.ReSize(npic+1);
			vtmp1.Rows(1,npic)=p0.Rows(1,npic);
			vtmp1.Row(npic+1)=1;
			p0.CleanUp();
			p0=vtmp1;
			vtmp1.CleanUp();
			vtmp2 = -1.0*constraint;
			a = a|vtmp2;
			vtmp2.CleanUp();
			cx.ReSize(npic+1);
			cx=0.0;
			cx.Column(npic+1)=1;
			dx.ReSize(npic+1);
			dx=1.0;
			go=1.0;
			minit = 0;

			while( go >= tol)
			{
				minit+=1;
				gap.ReSize(mm,2);
				mtmp1.ReSize(mm,2);
				mtmp1.SubMatrix(1,mm,1,1)= ( p0.Rows(1,mm)-pb.Column(1) );
				mtmp1.SubMatrix(1,mm,2,2)= ( pb.Column(2)-p0.Rows(1,mm) );
				gap=row2msort(mtmp1);
				dx.Rows(1,mm) = gap.Column(1);
				dx.Row(npic+1) = p0.Row(npic+1).AsScalar();
				mtmp1.CleanUp();
				if( lpb1 == 0 ) {
					vtmp1.ReSize(mm+1);
					vtmp2.ReSize(npic-mm);
					vtmp1.Rows(1,mm)=dx.Rows(1,mm);
					vtmp1.Row(mm+1)=100;
					dtmp1 = vtmp1.Maximum();
					vtmp2 = dtmp1;
					dx.Rows(mm + 1,npic)=vtmp2;
					vtmp2.CleanUp();
					vtmp1.CleanUp();
				}
				mtmp1 = vec2diag(dx);
				mtmp3 = (a*mtmp1).t();
				vtmp4 = dotcvecmult(dx, cx.t());
				y << solveqr(mtmp3,vtmp4);
				vtmp2 = cx.t()-(a.t()*y);
				vtmp3 = dotcvecmult(dx, vtmp2);
				v = dotcvecmult(dx, vtmp3);
				mtmp1.CleanUp();
				mtmp3.CleanUp();
				vtmp2.CleanUp();
				vtmp3.CleanUp();
				vtmp4.CleanUp();

				dtmp3=v.Row(npic + 1).AsScalar();
				if( dtmp3 > 0.0 ){
					z = p0.Row(npic + 1).AsScalar()/v.Row( npic + 1).AsScalar();
				for(i=1;i<=mm;i++){
					dtmp4=v.Row(i).AsScalar();
					if(dtmp4 < 0.0){
						dtmp1 = (pb.SubMatrix(i,i,2,2).AsScalar()-p0.Row(i).AsScalar())/v.Row(i).AsScalar();
						dtmp2=z;
						z = min(dtmp2 , -1.0*dtmp1);
					}
					else if(dtmp4 > 0.0){
						dtmp1 = (p0.Row(i).AsScalar()-pb.SubMatrix(i,i,1,1).AsScalar())/v.Row(i).AsScalar();
						dtmp2=z;
						z = min(dtmp2 , dtmp1);
					}
				}
				//defineVar(install(".z"), C2RD(z), rho);
				R_CheckUserInterrupt();
				dtmp1 = p0.Row(npic+1).AsScalar()/v.Row(npic+1).AsScalar();
				if( z >=  dtmp1){
					vtmp1 = p0 - (z * v);
					p0 = vtmp1;
					//defineVar(install(".p0"), CV2RV(p0), rho);
					vtmp1.CleanUp();
				} else{
					vtmp1 = p0 - (0.9 * z * v);
					p0 = vtmp1;
					//defineVar(install(".p0"), CV2RV(p0), rho);
					vtmp1.CleanUp();
				}
				go = p0.Row(npic+1).AsScalar();

				if(minit >=10){
					go  = 0.0;
				}

			} else{
				go = 0.0;
				minit = 10;
			}
		}

		if(minit>=10){
			Rprintf("SOLNP-->The linearized problem has no feasible "
					"solution.\n The problem may not be feasible\n");
		}

		mtmp1 = a.SubMatrix(1,arow,1,npic);
		a.CleanUp();
		a << mtmp1;
		b = a * p0.Rows(1,npic);
		mtmp1.CleanUp();
		}
	}
	p.CleanUp();
	p.ReSize(npic);
	p = p0.Rows(1,npic);
	z=0.0;
	y=0.0;

	if(ch > 0){
		tempv=dotcvecmult(p.Rows(nineq+1,npic),scale.Rows(nc + 2, nc + np + 1));
		ob=solnpeval(CV2RV(tempv), Ofn, Efn, Ifn, rho, xn);
		ob = div2cvec(ob,scale.Rows(1,nc+1));
		tempv.CleanUp();
	}

	j = ob.Row(1).AsScalar();

	if( nineq > 0){
		vtmp1 = ob.Rows(neq + 2,nc + 1) - p.Rows( 1, nineq );
		ob.Rows(neq + 2, nc + 1) = vtmp1;
		vtmp1.CleanUp();
	}

	if( nc > 0 ) {
		vtmp1 = (ob.Rows( 2, nc + 1) - (a * p)) +  b;
		ob.Rows( 2,nc + 1) = vtmp1;
		dtmp1 = ob.Rows( 2, nc + 1).SumSquare();
		vtmp3 = ( lagrange1.t() * ob.Rows( 2 , nc + 1));
		j = ob.Row(1).AsScalar() - vtmp3.Row(1).AsScalar()  + xrho * dtmp1;
		vtmp1.CleanUp();
		vtmp3.CleanUp();
	}

	minit = 0;

	while( minit < maxit ) {

		minit = minit + 1;

		if( ch > 0 ) {

			for(i=1; i<=np;i++) {
				dtmp2=p.Row( nineq + i ).AsScalar();
				p.Row( nineq + i ) = dtmp2 + delta;
				tempv=dotcvecmult(p.Rows(nineq+1,npic),scale.Rows(nc + 2, nc + np + 1));
				obm=solnpeval(CV2RV(tempv), Ofn, Efn, Ifn, rho, xn);
				obm = div2cvec(obm,scale.Rows(1,nc+1));
				tempv.CleanUp();

			if( nineq > 0){
				vtmp1 = obm.Rows(neq + 2,nc + 1) - p.Rows( 1, nineq );
				obm.Rows(neq + 2, nc + 1)= vtmp1;
				vtmp1.CleanUp();
			}

			if( nc > 0){
				vtmp1 = (obm.Rows(2, nc + 1) - (a * p)) + b;
				obm.Rows(2,nc + 1) = vtmp1;
				vtmp2 = lagrange1.t() * obm.Rows( 2 , nc + 1);
				dtmp1 = obm.Rows( 2, nc + 1).SumSquare();
				obmm = obm.Row(1).AsScalar() - vtmp2.Row(1).AsScalar()  + xrho * dtmp1;
				vtmp1.CleanUp();
			}
			g.Row(nineq + i) = (obmm - j)/delta;
			p.Row(nineq + i) = dtmp2;
		}

		if( nineq > 0 ) {
			g.Rows( 1, nineq ) = 0.0* lagrange1.Rows(neq + 1,nc );
			}
		}

		R_CheckUserInterrupt();
		if( minit > 1 ) {

			vtmp3 = g - yg;
			yg << vtmp3;
			vtmp4 = p - sx;
			sx << vtmp4;
			vtmp3.CleanUp();
			vtmp4.CleanUp();
			vtmp1 = sx.t() * (hess * sx);
			sc1 = vtmp1.Row(1).AsScalar();
			vtmp2 = sx.t() * yg;
			sc2 = vtmp2.Row(1).AsScalar();
			vtmp1.CleanUp();
			vtmp2.CleanUp();
			dtmp4=sc1*sc2;
			if( dtmp4 > 0.0 ) {
				vtmp1 = hess * sx;
				sx << vtmp1;
				mtmp1 = (sx * sx.t()) * (1.0/sc1);
				mtmp2 = (yg * yg.t()) * (1.0/sc2);
				mtmp3 = (hess - mtmp1);
				mtmp3 = mtmp3 + mtmp2;
				hess<<mtmp3;
				vtmp1.CleanUp();
				mtmp1.CleanUp();
				mtmp2.CleanUp();
				mtmp3.CleanUp();
			}
		}

		dx.CleanUp();
		dx.ReSize(npic);
		dx=0.01;

		if( lpb2 == 1 ) {
			gap.CleanUp();
			ColumnVector gap(mm);
			gap = 0.0;
			mtmp2.ReSize(mm,2);
			mtmp3.ReSize(mm,2);
			mtmp2 = 0.0;
			mtmp2.Column(1)=p.Rows(1,mm)-pb.Column(1);
			mtmp2.Column(2)=pb.Column(2)-p.Rows(1,mm);
			mtmp3=row2msort(mtmp2);
			gap = mtmp3.Column(1) + sqrt(EPSX);
			vtmp1.ReSize(mm);
			vtmp1=1.0;
			vtmp2=div2cvec(vtmp1,gap);
			dx.Rows(1,mm) = vtmp2;
			vtmp1.CleanUp();
			vtmp2.CleanUp();
			mtmp2.CleanUp();
			mtmp3.CleanUp();

			if( lpb1 == 0 ){
				vtmp1 = add1vec(dx.Rows(1,mm),0.01);
				vtmp2.ReSize(npic-mm);
				vtmp2=1.0;
				dx.Rows(mm+1,npic) = vtmp1.Minimum()*vtmp2;
				vtmp1.CleanUp();
				vtmp2.CleanUp();
			}
		}

		R_CheckUserInterrupt();
		go = -1.0;
		lagrange2 = lagrange2/10.0;

		while( go <= 0){
			vtmp1 = dotcvecmult(dx, dx);
			vtmp1 = lagrange2 * vtmp1;
			mtmp1 << vec2diag(vtmp1);
			mtmp2 = hess + mtmp1;
			/*S<< mtmp2;
			try{
				L = Cholesky(S);
			}
			catch(Exception){
				Rprintf(Exception::what());
			}*/
			//make it upper triangular
			SM = solnp_chol(C2Rmat(mtmp2));
			//mtmp4 << L.t();
			mtmp4 = R2Cmat(SM,Rf_nrows(SM),Rf_ncols(SM));
			cz = mtmp4.i();
			yg = cz.t() * g;
			/// check code here
			mtmp1.CleanUp();
			mtmp2.CleanUp();
			S.CleanUp();
			vtmp1.CleanUp();
			mtmp4.CleanUp();

			if( nc == 0 ) {
				u = (-cz) * yg;
			} else{
				mtmp3 = cz.t() * a.t();
				y << solveqr(mtmp3, yg);
				u = (-cz) * (yg - ( (cz.t() * a.t())  * y));
				mtmp3.CleanUp();
			}
			p0= (p + u.Rows(1,npic));

			if( lpb2 == 0 ) {
				go = 1.0;
			} else {
				vtmp1.ReSize(2*mm);
				vtmp1.Rows(1,mm)=p0.Rows(1,mm)- pb.Column(1);
				vtmp1.Rows(mm+1,2*mm)=pb.Column(2)- p0.Rows(1,mm);
				go = vtmp1.Minimum();
				lagrange2 = 3.0 * lagrange2;
				vtmp1.CleanUp();
			}
		}
		R_CheckUserInterrupt();
		alp1 = 0.0;
		ob1 = ob;
		ob2 = ob1;
		sob.Row(1) = j;
		sob.Row(2) = j;
		ptt = p|p;
		alp3 = 1.0;
		ptt = ptt|p0;
		tempv=dotcvecmult(ptt.SubMatrix(nineq+1,npic,3,3),scale.Rows(nc + 2, nc + np + 1));
		ob3 = solnpeval(CV2RV(tempv), Ofn, Efn, Ifn, rho, xn);
		ob3 = div2cvec(ob3,scale.Rows(1,nc+1));
		sob.Row(3) = ob3.Row(1).AsScalar();

		if( nineq > 0 ) {
			vtmp1 = ob3.Rows(neq + 2,nc + 1)- ptt.SubMatrix(1,nineq,3,3);
			ob3.Rows(neq + 2,nc + 1) = vtmp1;
			vtmp1.CleanUp();
		}
		if( nc > 0 ) {
			vtmp1 = (ob3.Rows( 2,nc + 1) - (a * ptt.Columns(3 ,3))) + b;
			ob3.Rows(2,nc + 1) = vtmp1;
			vtmp2 = lagrange1.t() * ob3.Rows( 2 , nc + 1);
			dtmp1 = ob3.Rows( 2, nc + 1).SumSquare();
			sob.Row(3) = ob3.Row(1).AsScalar() - vtmp2.Row(1).AsScalar()  + xrho * dtmp1;
			vtmp1.CleanUp();
			vtmp2.CleanUp();
			vtmp3.CleanUp();
		}
		go = 1;

		while( go > tol ) {
			alp2 = (alp1 + alp3) / 2.0;
			ptt.Columns(2,2) = ( ((1.0 - alp2) * p) + (alp2 * p0 ));
			tempv=dotcvecmult(ptt.SubMatrix(nineq+1,npic,2,2),scale.Rows(nc + 2, nc + np + 1));
			ob2 = solnpeval(CV2RV(tempv), Ofn, Efn, Ifn, rho, xn);
			ob2 = div2cvec(ob2,scale.Rows(1,nc+1));
			sob.Row(2) = ob2.Row(1).AsScalar();

			if( nineq > 0 ) {
				vtmp1 = ob2.Rows(neq + 2,nc + 1)- ptt.SubMatrix(1,nineq,2,2);
				ob2.Rows(neq + 2,nc + 1) = vtmp1;
				vtmp1.CleanUp();
			}
			if( nc > 0 ) {
				vtmp1 = (ob2.Rows( 2,nc + 1) - (a * ptt.Columns(2 ,2))) + b;
				ob2.Rows(2,nc + 1) =  vtmp1;
				vtmp2 = lagrange1.t() * ob2.Rows( 2 , nc + 1);
				dtmp1 = ob2.Rows( 2, nc + 1).SumSquare();
				sob.Row(2) = ob2.Row(1).AsScalar() - vtmp2.Row(1).AsScalar()  + xrho * dtmp1;
				vtmp1.CleanUp();
				vtmp2.CleanUp();
			}

			obmm = sob.Maximum();

			if( obmm < j ) {
				obn = sob.Minimum();
				go = tol * (obmm - obn) / (j - obmm);
			}

			R_CheckUserInterrupt();
			dtmp1=sob.Row(1).AsScalar();
			dtmp2=sob.Row(2).AsScalar();
			dtmp3=sob.Row(3).AsScalar();
			if( dtmp2 >= dtmp1){
				sob.Row(3) = sob.Row(2).AsScalar();
				ob3 = ob2;
				alp3 = alp2;
				ptt.Column(3) = ptt.Column(2);
			} else if(dtmp1 <= dtmp3) {
				sob.Row(3) = sob.Row(2).AsScalar();
				ob3 = ob2;
				alp3 = alp2;
				ptt.Column(3) = ptt.Column(2);
			} else{
				sob.Row(1) = sob.Row(2).AsScalar();
				ob1 = ob2;
				alp1 = alp2;
				ptt.Column(1) = ptt.Column(2);
			}

			if( go >= tol ) {
				go = alp3 - alp1;
			}

		}

		sx << p;
		yg << g;
		ch = 1.0;
		obn = sob.Minimum();

		if( j <=obn ) {
			maxit = minit;
		}

		reduce = (j - obn) / ( 1 + fabs(j) );
		if( reduce < tol ) {
			maxit = minit;
		}

		dtmp1=sob.Row(1).AsScalar();
		dtmp2=sob.Row(2).AsScalar();
		dtmp3=sob.Row(3).AsScalar();

		if(dtmp1 <  dtmp2) {
			j = sob.Row(1).AsScalar();
			p << ptt.Column(1);
			ob = ob1;
		} else if(dtmp3 <  dtmp2) {
			j = sob.Row(3).AsScalar();
			p << ptt.Column(3);
			ob = ob3;
		} else{
			j = sob.Row(2).AsScalar();
			p << ptt.Column(2);
			ob = ob2;
		}
		ob1.CleanUp();
		ob2.CleanUp();
		ob3.CleanUp();
		ptt.CleanUp();
	}
	R_CheckUserInterrupt();
	vtmp1=dotcvecmult(p, scale.Rows(neq + 2,nc + np + 1));
	p << vtmp1;
	vtmp1.CleanUp();

	//unscale the lagrange multipliers
	if( nc > 0 ){
		vtmp1 = div2cvec(y , scale.Rows( 2,nc + 1));
		y << scale.Row(1).AsScalar()*vtmp1;
		vtmp1.CleanUp();
	}

	mtmp1 = scale.Rows(neq + 2,nc + np + 1)*scale.Rows(neq + 2,nc + np + 1).t();
	mtmp2 = div2mat(hess,mtmp1);
	hess = scale.Row(1).AsScalar()*mtmp2;
	mtmp1.CleanUp();
	mtmp2.CleanUp();
	if( reduce > tol ) {
		Rprintf("SOLNP-->Minor optimization routine did not converge in the \n"
				"specified number of minor iterations.  You may need\n"
				"to increase the number of minor iterations.        \n");
	}
	/*Cleanup and return output*/
	/*pos 1 : pars*/
	REAL(MU)[0]=lagrange2;
	SET_VECTOR_ELT(outvec,0, CV2RV(p));
	/*pos 2 : lagrange1-out*/
	SET_VECTOR_ELT(outvec,1, CV2RV(y));
	/*pos 3 : hessian-out*/
	SET_VECTOR_ELT(outvec,2, C2Rmat(hess));
	/*pos 4 : mu-out*/
	SET_VECTOR_ELT(outvec,3,MU);
	UNPROTECT(3);
	/*free newmat*/
	lagrange1.Release();
	sob.Release(); hess.Release(); scale.Release();
	pb.Release(); gap.Release(); gap2.Release(); ptt.Release();
	ob.Release(); obm.Release(); cx.Release();
	dx.Release(); yg.Release(); sx.Release(); cz.Release();
	u.Release(); ob1.Release(); ob2.Release(); ob3.Release();
	p.Release(); p0.Release(); constraint.Release(); g.Release();
	mtmp1.Release(); mtmp2.Release(); mtmp3.Release(); mtmp4.Release(); mtmp5.Release();
	a.Release(); b.Release(); y.Release(); v.Release();
	vtmp1.Release(); vtmp2.Release(); vtmp3.Release(); vtmp4.Release();
	S.Release(); L.Release();
	return(outvec);
	/*}
	catch(Exception){
			Rprintf(Exception::what());
	}*/
}

ReturnMatrix solnpeval(SEXP pars, SEXP Ofn, SEXP Efn, SEXP Ifn, SEXP rho, SEXP n)
{
	if(!isEnvironment(rho)) error("'rho' should be an environment");
	SEXP fx;
	int k;
	PROTECT(fx = allocVector(REALSXP, 1));
	defineVar(install(".par"), pars, rho);
	int neq=int(REAL(n)[1]);
	int nineq=int(REAL(n)[2]);
	k=1;
	ColumnVector ob(1+nineq+neq);
	ob=0.0;
	fx=eval(Ofn, rho);
	ob.Row(1)=(R2CCV(fx)).AsScalar();
	if(neq>=1){
		SEXP ex;
		PROTECT(ex = allocVector(REALSXP, neq));
		k+=1;
		ex = eval(Efn, rho);
		ob.Rows(2,neq+1)=R2CCV(ex);
	}
	if(nineq>=1){
		SEXP ix;
		PROTECT(ix = allocVector(REALSXP, nineq));
		k+=1;
		ix = eval(Ifn, rho);
		ob.Rows(2+neq,nineq+neq+1)=R2CCV(ix);
	}
	//defineVar(install(".par"), CV2RV(newpars), rho);
	UNPROTECT(k);
	ob.Release();
	return ob.ForReturn();
}
// solve condition for svd function needed.
}
