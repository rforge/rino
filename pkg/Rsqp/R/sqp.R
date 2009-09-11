#################################################################################
##
##   R package rsqp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009
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
		LB = NULL, UB = NULL,  control = list(maxiter = 100, tol = 1e-8), ...)
{
# main function (coming soon)
}


# Checks and other related functions
.checkpars = function(pars, LB, UB)
{
	if(is.null(pars))
		stop("\nsqp-->error: must supply starting parameters\n", call. = FALSE)
	if(!is.null(LB)){
		if(length(pars)!=length(LB))
			stop("\nsqp-->error: LB length not equal to parameter length\n", call. = FALSE)
	} else{
		LB = rep(-Inf, length(pars))
	}
	if(!is.null(UB)){
		if(length(pars)!=length(UB))
			stop("\nsqp-->error: UB length not equal to parameter length\n", call. = FALSE)
	} else{
		UB = rep(Inf, length(pars))
	}
	if(any(LB>UB))
		stop("\nsqp-->error: UB must be greater than LB\n", call. = FALSE)
	return(list(LB = LB, UB = UB))
}

.checkfun = function(pars, fun, ...)
{
	val = fun(pars, ...)
	if(length(val)!=1)
		stop("\nsqp-->error: objective function returns value of length greater than 1!\n", call. = FALSE)
	
	.sqp_fun <<- fun
	return(val)
}

.checkgrad = function(pars, fun, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(val)!=n)
		stop("\nsqp-->error: gradient vector length must be equal to length(pars)\n", call. = FALSE)
	
	.sqp_fungrad <<- fun
	return(val)
}

.checkhess = function(pars, fun, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(as.vector(val))!=(n*n))
		stop("\nsqp-->error: hessian must be of length length(pars) x length(pars)\n", call. = FALSE)
	
	.sqp_funhess <<- fun
	return(val)
}

.checkeq = function(pars, fun, eq, ...)
{
	n = length(eq)
	val = fun(pars, ...)
	if(length(val)!=n)
		stop("\nsqp-->error: equality function returns vector of different length
to equality value\n", call. = FALSE)

	.sqp_eqfun <<- fun
	return(val)
}

.checkineq = function(pars, fun, ineqLB, ineqUB, ...)
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
	
	.sqp_ineqfun <<- .ineqtransform(pars, fun, ineqLB, ineqUB, ...)
	return(val)
}

# reporting function
.report = function(iter, qp_iter, alpha, nfun, fun)
{
	cat(paste("\n",iter," ",qp_iter," ",alpha," ",nfun," ",obj," ",sep = ""))
}

# finite difference gradient
.fdgrad = function(pars, fun, ...)
{
	if(!is.null(fun))
	{
		y0 = fun(pars, ...)
		nx = length(pars)
		grd = rep(0,nx)
		deltax = sqrt(.eps)
		for(i in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			grd[i] = (fun(pars, ...) - y0) / deltax
			x[i] = init
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
			x[i] = init
		}
	} else{
		jac = rep(0, nx)
	}
	return(jac)
}

phiL1 = function(pars, fun, ineqfun, eqfun, mu, funv, ...)
{
	if(!is.null(eqfun)) eqv = eqfun(pars, ...) else eqv = NULL
	if(!is.null(eqfun))
	{
		ineqv = ineqfun(pars, ...) 
		idx = ineq[ineqv<0]
		ineqv = ineqv[idx]
	} else 
	{
		ineqv = NULL
	}
	con = c(eqv, ineqv)
	if(is.null(funv))
	{
		funv = fun(pars, ...)
		.nfun <<- .nfun+1
	}
	merit = funv
	tmp = .norm(con, 1) / mu
	if(!is.null(tpm)) merit = merit + tmp
	return(list(merit = merit, funv = funv))
}

.ineqtransform = function(pars, ineqfun, ineqLB, ineqUB, ...)
{
	# transform from:
	# ineqLB =< ineqfun(x) =< ineqUB
	# to:
	# ineqfun(x) >= 0
	# this is also consisent with the solve.QP function which will be
	# called by sqp
	val = ineqfun(pars, ...)
	retval = c(ineqUB - val, val - ineqLB)
	return(retval)
}