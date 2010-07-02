#################################################################################
##
##   R package rsqp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2008
##   This file is part of the R package rsqp.
##
##   The R package rsqp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rsqp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

sqp = function(pars, fun, grad = NULL, hess = NULL, ineqfun = NULL, gradineq = NULL, 
		ineqLB = NULL,  ineqUB = NULL, eqfun = NULL, gradeq = NULL, eqB = NULL, 
		LB = NULL, UB = NULL,  control = list(maxiter = 100, tol = 1e-6, trace = 0, 
				initiate.with.solnp = TRUE, solnp.iter = 2), ...)
{
	
	# make one pass from solnp to get improved parameters:
	if(!is.null(control$tol)) tmptol = control$tol else tmptol = 1e-6
	if(!is.null(control$solnp.iter)) solnp.iter = control$solnp.iter else solnp.iter = 2
	if(is.null(control$initiate.with.solnp)) initiate.with.solnp = FALSE else initiate.with.solnp = control$initiate.with.solnp
	if(initiate.with.solnp){ 
		tmps = try(solnp( pars = pars, fun = fun, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
			eqfun = eqfun, eqB = eqB, LB = LB, UB = UB, control = list(outer.iter = solnp.iter, 
					tol = tmptol, trace = 0), ...), silent = TRUE )
	if( !inherits(tmps, "try-error") ) pars = tmps$pars
	}
	tic = Sys.time()
	.sqpenv <- environment()
	
	.sqp_nfun  <<- 0
	# lower/upper bounds and pars check
	check1 = .checkpars(pars, LB, UB, .sqpenv)
	# .sqp_LB and .sqp_UB vectors assigned
	
	# initial function check
	funv = .checkfun(pars, fun, .sqpenv, ...)
	# .sqp_fun(pars, ...) assigned

	# information indicator
	ind = vector(mode = "numeric", length = 10)
	# [1] length of pars
	# [2] has function gradient?
	# [3] has hessian?
	# [4] has ineq?
	# [5] ineq length
	# [6] has jacobian (inequality)
	# [7] has eq?
	# [8] eq length
	# [9] has jacobian (equality)
	ind[1] = length(pars)
	
	# gradient and hessian checks
	if(!is.null(grad)){
		gradv = .checkgrad(pars, grad, .sqpenv, ...)
		ind[2] = 1
	} else{
		.sqp_gradfun = function(pars, ...) .fdgrad(pars, fun = .sqp_fun, ...)
		ind[2] = 0
		gradv = .sqp_gradfun(pars, ...)
	}
	# .sqp_gradfun(pars, ...) assigned
	
	if(!is.null(hess)){
		hessv = .checkhess(pars, grad, .sqpenv, ...)
		.sqp_hessfun = hess
		ind[3] = 1
	} else {
		.sqp_hessfun = NULL
		ind[3] = 0
		hessv = diag(ind[1])
	}
	# .sqp_hessfun(pars, ...) assigned
	
	# inequality checks
	if(!is.null(ineqfun)){
		ineqv = .checkineq(pars, ineqfun, ineqLB, ineqUB, .sqpenv, ...)
		ind[4] = 1
		nineq = 2 * length(ineqLB)
		ind[5] = nineq
		if(!is.null(gradineq)){
			ineqjacv = .cheqjacineq(pars, gradineq, ineqUB, ineqLB, .sqpenv, ...)
			ind[6] = 1
		} else{
			.sqp_ineqjac = function(pars, ...) .fdjac(pars, fun = .sqp_ineqfun, ...)
			ind[6] = 0
			ineqjacv = .sqp_ineqjac(pars, ...)
		}
	} else {
		ineqv = NULL
		.sqp_ineqfun = function(x, ...) .emptyfun(x, ...)
		.sqp_ineqjac = function(x, ...) .emptyjac(x, ...)
		nineq = 0
		ind[4] = 0
		ind[5] = 0
		ind[6] = 0
	}
	# .sqp_ineqfun(pars, ...) assigned
	# .sqp_ineqjac(pars, ...) assigned

	# equality checks
	if(!is.null(eqfun)){
		eqv = .checkeq(pars, eqfun, eqB, .sqpenv, ...)
		ind[7] = 1
		neq = length(eqB)
		ind[8] = neq
		if(!is.null(gradeq)){
			eqjacv = .cheqjacineq(pars, gradineq, ineqUB, ineqLB, .sqpenv, ...)
			ind[9] = 1
		} else{
			.sqp_eqjac = function(pars, ...) .fdjac(pars, fun = .sqp_eqfun, ...)
			eqjacv = .sqp_eqjac(pars, ...)
			ind[9] = 0
		}
	} else {
		eqv = NULL
		eqjacv = NULL
		.sqp_eqfun = function(x, ...) .emptyfun(x, ...)
		.sqp_eqjac = function(x, ...) .emptyjac(x, ...)
		ind[7] = 0
		neq = 0
		ind[8] = 0
		ind[9] = 0
	}
	# .sqp_eqfun(pars, ...) assigned
	# .sqp_eqjac(pars, ...) assigned
	
	# Total Constraints (ncon) : neq + nineq + LB + UB
	nineqlb = nineq + 2 * length(.sqp_LB)
	ncon = neq + nineqlb
	# finally assign and merge ineqfun and ineqfun gradient to include UB/LB bounds
	# if not null
	if(!is.null(.sqp_LB)){
		.sqp_ineqlbfun = function(pars, ...) .ineqlbfun(pars, .env = .sqpenv, ...)
		.sqp_ineqlbjac = function(pars, ...) .ineqlbjac(pars, .env = .sqpenv, ...)
	} else{
		.sqp_ineqlbfun = .sqp_ineqfun
		.sqp_ineqlbjac = .sqp_ineqjac 
	}
	# assign index
	.sqp_index = ind

	# start optimization
	iter = 0
	if(is.null(control$maxiter)) maxiter = 100 else maxiter = control$maxiter
	if(is.null(control$tol)) tol = 1e-8 else tol = control$tol
	if(is.null(control$trace)) trace = 0 else trace = control$trace
	# funv
	.sqp_info = 0
	.qp_iter =  0
	# funv = obj
	# gradv  = c
	# hessv = B
	# eqv = ce
	# eqjacv = F
	# ineqv = ci
	# ineqjac = C
	ineqv = .sqp_ineqlbfun(pars, ...)
	ineqjacv = .sqp_ineqlbjac(pars, ...)
	
	A = rbind(eqjacv, ineqjacv)
	
	## Choose an initial lambda (pars are provided by the caller)
	
	lambda = 2 * .ones(ncon, 1)
	qp_iter = 1
	alpha = 1	
	info = 0
	while (iter < maxiter)
	{
		## Check convergence.
		if(neq>0) lambda_eq = lambda[1:neq] else lambda_e = NULL
		
		if(nineqlb>0) lambda_ineq = lambda[(neq+1):ncon] else lambda_ineq = NULL
		conv = c(eqv, ineqv)
		
		t0 = .vnorm (gradv - t(A)%*%lambda)
		if(neq>0) t1 = .vnorm (eqv) else t1 = NULL
		if(nineqlb>0) t2 = all (ineqv >= 0) else t2 = 0
		if(nineqlb>0) t3 = all (lambda_ineq >= 0) else t3 = 0
		t4 = .vnorm(lambda*conv)
		if (t2 && t3 && max (c(t0,t1,t4) < tol)) break()
		## Compute search direction p by solving QP.
		if(neq > 0) g = -eqv else g = NULL
		if(nineqlb > 0) d = -ineqv else d = NULL
		Amat = rbind(eqjacv, ineqjacv)
		bvec = c(g, d)
		solqp = try(solve.QP(Dmat = hessv, dvec = -gradv, Amat = t(Amat), bvec = bvec, meq = ind[8]),
				silent = TRUE)
		if(inherits(solqp,"try-error"))
			stop("\nsqp-->error: qp part infeasible\n", call. = FALSE)

		p = solqp$solution
		objqp = solqp$value
		.qp_iter = .qp_iter + solqp$iterations[1]
		# we need to calculate the lagrange multipliers which are not explicitly
		# returned to R!
		tmpx = hessv %*% p - (-gradv)
		Activ = t(Amat)[,solqp$iact]
		tmplambda= solve(crossprod(Activ,Activ), crossprod(Activ, tmpx))
		lambda = rep(0, ncol(t(Amat)))
		lambda[solqp$iact] = tmplambda
		
		tmp = .constrlinesearch(pars, p, lambda = lambda, funv = funv, .sqpenv,  ...)
		parsnew = tmp$xnew
		alpha = tmp$alpha
		funvnew = tmp$funv
		## Evaluate objective function, constraints, and gradients at
		## x_new.
		
		gradvnew = .sqp_gradfun(parsnew, ...)
		eqvnew = .sqp_eqfun(parsnew, ...)
		eqjacvnew = .sqp_eqjac(parsnew, ...)
		ineqvnew = .sqp_ineqlbfun(parsnew, ...)
		ineqjacvnew = .sqp_ineqlbjac(parsnew, ...)
		Anew = rbind(eqjacvnew, ineqjacvnew)
		
		## Set
		##
		## s = alpha * p
		## y = grad_x L (x_new, lambda) - grad_x L (x, lambda})
		y = gradvnew - gradv;
		
		if (!is.null(A)){
			xt = (t(Anew - A)%*%lambda)
			y = y - xt
		}
		delx = parsnew - pars
		if (.vnorm (delx) < tol * .vnorm (pars)){
			.sqp_info <<- 101
			break()
		}
		if (ind[3]){
			hessv = .sqp_hessfun(pars, ...)
		} else{
			delxt = t(delx)
			d1 = as.numeric(delxt%*%hessv%*%delx)
			t1 = 0.2 * d1
			t2 = as.numeric(delxt%*%y)
			if (t2 < t1){
				theta = 0.8*d1/(d1 - t2)
			} else{
				theta = 1
			}
			r = theta*y + (1-theta)*hessv%*%delx
			d2 = as.numeric(delxt%*%r)
			if (d1 == 0 || d2 == 0){
				.sqp_info <<- 102
				break()
			}
			hessv = hessv - (hessv%*%delx%*%delxt%*%hessv)/d1 + (r%*%t(r))/d2
		}
		pars = parsnew
		funv = funvnew
		gradv = gradvnew
		
		eqv = eqvnew
		eqjacv = eqjacvnew
		ineqv = ineqvnew
		ineqjacv = ineqjacvnew
		A = Anew
		iter = iter + 1
		if(trace) .report(iter, .qp_iter, funv)
	}
	cat("\n\n")
	if (iter >= maxiter) .sqp_info = 103
	nf = .sqp_nfun
	elapsed = Sys.time() - tic

	solution = list(pars = pars, obj = funv, info = .sqp_info, iter = iter,
			nfun.calls = nf, lambda = lambda, elapsed = elapsed)
	return(solution)
}



.checkpars = function(pars, LB, UB, .env)
{
	if(is.null(pars))
		stop("\nsqp-->error: must supply starting parameters\n", call. = FALSE)
	if(!is.null(LB)){
		if(length(pars)!=length(LB))
			stop("\nsqp-->error: LB length not equal to parameter length\n", call. = FALSE)
		if(is.null(UB)) UB = rep(.Machine$double.xmax/2, length(LB))
	} else{
		LB = NULL
	}
	if(!is.null(UB)){
		if(length(pars)!=length(UB))
			stop("\nsqp-->error: UB length not equal to parameter length\n", call. = FALSE)
		if(is.null(LB)) LB = rep(-.Machine$double.xmax/2, length(UB))
	} else{
		UB = NULL
	}
	if(!is.null(UB) && any(LB>UB))
		stop("\nsqp-->error: UB must be greater than LB\n", call. = FALSE)
	
	if(!is.null(UB) && any(LB==UB))
		warning("\nsqp-->warning: Equal Lower/Upper Bounds Found. Consider\n
						excluding fixed parameters.\n", call. = FALSE)
	# deal with infinite values as these are not accepted by solve.QP
	if(!is.null(LB) && !any(is.finite(LB))){
		idx = which(!is.finite(LB))
		LB[idx] = sign(LB[idx])*.Machine$double.xmax/2
	}
	if(!is.null(UB) && !any(is.finite(UB))){
		idx = which(!is.finite(UB))
		UB[idx] = sign(UB[idx])*.Machine$double.xmax/2
	}	
	assign(".sqp_LB", LB, envir = .env)
	assign(".sqp_UB", UB, envir = .env)
	return(list(LB = LB, UB = UB))
}

.checkfun = function(pars, fun, .env, ...)
{
	val = fun(pars, ...)
	if(length(val)!=1)
		stop("\nsqp-->error: objective function returns value of length greater than 1!\n", call. = FALSE)
	assign(".sqp_fun", fun, envir = .env)
	.sqp_nfun <<- .sqp_nfun + 1
	return(val)
}

.checkgrad = function(pars, fun, .env, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(val)!=n)
		stop("\nsqp-->error: gradient vector length must be equal to length(pars)\n", call. = FALSE)
	assign(".sqp_gradfun", fun, envir = .env)
	return(val)
}

.checkhess = function(pars, fun, .env, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(as.vector(val))!=(n*n))
		stop("\nsqp-->error: hessian must be of length length(pars) x length(pars)\n", call. = FALSE)
	assign(".sqp_hessfun", fun, envir = .env)
	return(val)
}

.checkeq = function(pars, fun, eqB, .env, ...)
{
	n = length(eqB)
	val = fun(pars, ...)
	if(length(val)!=n)
		stop("\nsqp-->error: equality function returns vector of different length
to equality value\n", call. = FALSE)
	.sqp_eqfun = function(x, ...) fun(x, ...) - eqB
	assign(".sqp_eqB", eqB, envir = .env)
	assign(".sqp_eqfun", .sqp_eqfun, envir = .env)
	return(val)
}

.checkineq = function(pars, fun, ineqLB, ineqUB, .env, ...)
{
	val = fun(pars, ...)
	n = length(val)
	if(!is.null(ineqLB)){
		if(length(ineqLB)!=n)
			stop("\nsqp-->error: inequality function returns vector of different length to
inequality lower bounds\n", call. = FALSE)
	} else{
		stop("\nsqp-->error: inequality function given without lower bounds\n", call. = FALSE)
	}
	if(!is.null(ineqUB)){
		if(length(ineqUB)!=n)
			stop("\nsqp-->error: inequality function returns vector of different length to
inequality upper bounds\n", call. = FALSE)
	} else{
		stop("\nsqp-->error: inequality function given without upper bounds\n", call. = FALSE)
	}
	if(any(ineqLB>ineqUB))
		stop("\nsqp-->error: ineqUB must be greater than ineqLB\n", call. = FALSE)
	# transform from:
	# ineqLB =< ineqfun(x) =< ineqUB
	# to:
	# ineqfun(x) >= 0
	# this is also consisent with the solve.QP function which will be
	# called by sqp
	.sqp_ineqLB = ineqLB
	.sqp_ineqUB = ineqUB
	assign(".sqp_ineqLB", .sqp_ineqLB, envir = .env)
	assign(".sqp_ineqUB", .sqp_ineqUB, envir = .env)
	.sqp_ineqfun = function(x, ...) { retval = fun(x, ...); c(retval - .sqp_ineqLB, .sqp_ineqUB - retval) }
	assign(".sqp_ineqfun", .sqp_ineqfun, envir = .env)
	return(val)
}

# check the jacobian of inequality
.cheqjacineq = function(pars, fun, .env,  ...)
{
	# must be a matrix -> nrows = no.inequalities, ncol = length(pars)
	val = fun(pars, ...)
	if(!is.matrix(val))
		stop("\nsqp-->error: Jacobian of Inequality must return a matrix type object\n", call. = FALSE)
	nd = dim(val)
	if(nd[2]!=length(pars))
		stop("\nsqp-->error: Jacobian of Inequality column dimension must be equal to length
of parameters\n", call. = FALSE)
	if(nd[1]!=length(.sq_ineqUB))
		stop("\nsqp-->error: Jacobian of Inequality row dimension must be equal to length
						of inequality bounds vector\n", call. = FALSE)
	# as in inequality function, transforms from a 2 sided inequality to a one sided inequality
	# (for the jacobian).
	.sqp_ineqjac = function(x, ...) { retval = fun(x, ...); rbind( - retval, retval ) }
	assign(".sqp_ineqjac", .sqp_ineqjac, envir = .env)
	return(val)
}

# check the jacobian of equality
.cheqjaceq = function(pars, fun, .env, ...)
{
	# must be a matrix -> nrows = no.equalities, ncol = length(pars)
	val = fun(pars, ...)
	if(!is.matrix(val))
		stop("\nsqp-->error: Jacobian of Equality must return a matrix type object\n", call. = FALSE)
	nd = dim(val)
	if(nd[2]!=length(pars))
		stop("\nsqp-->error: Jacobian of Equality column dimension must be equal to length
						of parameters\n", call. = FALSE)
	if(nd[1]!=length(.sqp_eqB))
		stop("\nsqp-->error: Jacobian of Equality row dimension must be equal to length
						of equality bounds vector\n", call. = FALSE)
	assign(".sqp_eqjac", fun, envir = .env)
	return(val)
}

# reporting function
.report = function(iter, qpiter, funv)
{
	cat(paste("\nouter_iter:\t ", iter, " qp_iter:\t ", qpiter, " objective:\t", funv, sep = ""))
}

# finite difference gradient
.fdgrad = function(pars, fun, ...)
{
	if(!is.null(fun)){
	
		y0 = fun(pars, ...)
		nx = length(pars)
		grd = rep(0, nx)
		deltax = sqrt(.eps)
		for(i in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			grd[i] = (fun(pars, ...) - y0) / deltax
			pars[i] = init
		}
	}
	else
	{
		grd = 0
	}
	return(grd)
}

# finite difference jacobian
.fdjac = function(pars, fun, ...)
{
	nx = length(pars)
	if(!is.null(fun))
	{
		y0 = fun(pars, ...)
		nf = length (y0)
		jac = matrix(0, nrow = nf, ncol= nx)
		deltax = sqrt (.eps)
		for(i  in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			jac[,i] = (fun(pars, ...) - y0) / deltax
			pars[i] = init
		}
	} else{
		jac = rep(0, nx)
	}
	return(jac)
}

.meritfn = function(pars, mu, funv = NULL, .env, ...)
{
	fun = get(".sqp_fun", envir = .env)
	eqfun = get(".sqp_eqfun", envir = .env)
	ineqfun = get(".sqp_ineqlbfun", envir = .env)
	idx = get(".sqp_index", envir = .env)
	if(idx[7]) eqv = eqfun(pars, ...) else eqv = NULL
	if(idx[4])
	{
		ineqv = ineqfun(pars, ...) 
		ineqv = ineqv[ineqv<0]
	} else {
		ineqv = NULL
	}
	con = c(eqv, ineqv)
	if(is.null(funv))
	{
		funv = fun(pars, ...)
		.sqp_nfun <<- .sqp_nfun+1
	}
	merit = funv
	tmp = .colnorm(con) / mu
	if(!is.null(tmp)) merit = merit + tmp
	return(list(merit = merit, funv = funv))
}

.emptygrad = function(pars, ...)
{
	matrix(0, nrow = 0, ncol = 1)
}

.emptyjac = function(pars, ...)
{
	matrix(0, nrow = 0, ncol = length(pars))
}

.emptyfun = function(pars, ...)
{
	NULL
}

.ineqlbfun = function(pars, .env, ...)
{
	LB = get(".sqp_LB", envir = .env)
	UB = get(".sqp_UB", envir = .env)
	.sqp_ineqfun = get(".sqp_ineqfun", envir = .env)
	res = c(pars - LB,  UB - pars)
	if(!is.null(.sqp_ineqfun)) res = c(.sqp_ineqfun(pars, ...), res)
	res
}

.ineqlbjac = function(pars, .env, ...)
{
	.sqp_ineqjac = get(".sqp_ineqjac", envir = .env)
	n = length(pars)
	res = rbind(diag(n), -diag(n))
	if(!is.null(.sqp_ineqjac)) res = rbind(.sqp_ineqjac(pars, ...), res)
	res
}

.zeros<-function(n=1,m=1)
{
	if(missing(m)) m=n
	sol<-matrix(0,nrow=n,ncol=m)
	return(sol)
}

.ones<-function(n=1,m=1)
{
	if(missing(m)) m=n
	sol<-matrix(1,nrow=n,ncol=m)
	return(sol)
}

.vnorm = function(x){
	sum((x)^2)^(1/2)
}

.rownorm = function(x)
{
	max(apply(as.matrix(x), 1, FUN=function(x) sum(abs(x)) ) )
}

.colnorm = function(x)
{
	max(apply(as.matrix(x), 2, FUN=function(x) sum(abs(x)) ) )
}

.solvecond<-function(x)
{
	z=svd(x)$d
	if(any(z==0)) ret=Inf else ret=max(z)/min(z)
	return(ret)
}

.eps = .Machine$double.eps