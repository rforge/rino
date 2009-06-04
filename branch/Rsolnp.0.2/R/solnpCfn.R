#################################################################################
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

.solnp <- function(pars, Jfun, Efun=NULL, EQ=NULL, Ifun=NULL, ILB=NULL, IUB=NULL, LB=NULL, UB=NULL, control=list(), ...)
{
	
	# start the timer
	tic = Sys.time()
	
	# get length of objective function parameters
	np  = length(pars)
	
	# lower parameter bounds - indicator
	# lpb[1]=1 means lower/upper bounds present
	# lpb[2]=1 means lower/upper bounds OR inequality bounds present
	lpb = c(0, 0)
	
	# parameter checks
	#-------------------------------------
	check1 = .checkpars(pars, LB, UB)
	
	if( !is.null(check1) )  return( check1 )
	
	if( !is.null(LB) || !is.null(UB) ) lpb[ 1 ] = 1
	
	# inequality checks
	#-------------------------------------
	check2 =.checkineq(pars, Ifun, ILB, IUB, ...)
	
	if( !is.null(check2$message) && !is.null(Ifun) ) return( check2$message )
	
	# no. inequalities
	nineq  = check2$nineq
	xineq0 = check2$xineq0
	if( lpb[ 1 ] == 1 || nineq > 0) lpb[ 2 ] = 1
	
	# equality checks
	#-------------------------------------
	check3 = .checkeq(pars, Efun, EQ, ...)
	
	if( !is.null(check3$message) ) return( check3$message )
	
	# no. on inequalities (neq)
	neq = check3$neq
	
	# parameter bounds (pb)
	# 1:nineq [Inequality] ; (nineq+1):(nineq+np)  [Parameter]
	pb  = rbind( cbind(ILB, IUB), cbind(LB, UB) )
	
	# check control list
	#-------------------------------------
	ctrl  = .checkcontrol( control )
	rho   = ctrl[[ 1 ]]
	maxit = ctrl[[ 2 ]]
	minit = ctrl[[ 3 ]]
	delta = ctrl[[ 4 ]]
	tol   = ctrl[[ 5 ]]
	trace = ctrl[[ 6 ]]
	
	# total constraints (tc) = no.inequality constraints + no.equality constraints
	tc = nineq + neq
	
	# initialize fn value and inequalities and set to NULL those not needed
	startf = Jfun(pars, ...)
	
	if( nineq>0 ) starti = Ifun(pars,...) else starti = NULL
	
	if( neq>0 ) starte = Efun(pars,...) - EQ else starte = NULL
	
	j  = startf
	jh = startf
	tt = 0 * .ones(3, 1)
	
	## Establish the objective function and its environment
	#-------------------------------------
	ofn <- quote(Jfun(.par, ...))
	if(!is.null(Efun)) efn <- quote(Efun(.par, ...)) else efn=NULL
	if(!is.null(Ifun)) ifn <- quote(Ifun(.par, ...)) else ifn=NULL
	.solenv <- new.env(parent = environment())
	assign(".par", pars, envir = .solenv)
	if( tc>0 ) {
		# lagrange multipliers (l)
		l = 0 * .ones(tc, 1)
		# constraint = [1:neq 1:nineq]
		constraint = c(starte, starti)
		
		if( nineq>0 ) {
			tmpv = cbind(constraint[ (neq + 1):tc ] - ILB, IUB - constraint[ (neq + 1):tc ] )
			testmin = apply( tmpv, 1, FUN = function( x ) min(x[ 1 ], x[ 2 ]) )
			if( all(testmin > 0) ) xineq0 = constraint[ (neq + 1):tc ]
			constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - xineq0
		}
		
		tt[ 2 ] = .vnorm(constraint)
		if( max(tt[ 2 ] - 10 * tol, nineq) <= 0 ) rho = 0
	} else{
		l = 0
	}
	
	p  = c(xineq0, pars)
	h  = diag(np + nineq)
	mu = np
	iteration = 0
	ob = c(startf, starte, starti)
	hess=matrix(0,np,np)
	assign("ob", ob, .solenv)
	assign("hess", hess, .solenv)
	assign("lgmult",l,.solenv)
	assign("mu",mu,.solenv)
	assign("constraint",constraint,.solenv)
	.Call(R_solnp, ofn, efn, ifn, ebounds, ibound, xbounds, control)
}
