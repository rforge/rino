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

# Based on the original solnp by Yinyu Ye
# http://www.stanford.edu/~yyye/Col.html


#----------------------------------------------------------------------------------
# The Function SOLNP solves nonlinear programs in standard form:
#
#        minimize              J(P)
#        subject to            EC(P)  =0
#                   IB(:,1)<=  IC(P)  <=IB(:,2)
#                   PB(:,1)<=    P    <=PB(:,2).
#where
#
#  J       : Cost objective scalar function
#  EC      : Equality constraint vector function
#  IC      : Inequality constraint vector function
#  P       : Decision parameter vector
#  IB, PB  : lower and upper bounds for IC and P.
#----------------------------------------------------------------------------------

# control list
#           RHO  : penalty parameter
#           MAJIT: maximum number of major iterations
#           MINIT: maximum number of minor iterations
#           DELTA: relative step size in forward difference evaluation
#           TOL  : tolerance on feasibility and optimality
# defaults RHO=1, MAJIT=10, MINIT=10, DELTA=1.0e-5, TOL=1.0e-4

solnp <- function(pars, Jfun, Efun=NULL, EQ=NULL, Ifun=NULL, ILB=NULL, IUB=NULL, LB=NULL, UB=NULL, control=list(), ...)
{
	# start timer
	tic = Sys.time()
	
	np  = length(pars)
	
	# lower parameter bounds - initialize
	lpb = c(0, 0)
	
	# do parameter checks
	check1 = .checkpars(pars, LB, UB)
	
	if( !is.null(check1) )  return( check1 )
	
	if( !is.null(LB) || !is.null(UB) ) lpb[ 1 ] = 1
	
	# do inequality checks
	check2 =.checkineq(pars, Ifun, ILB, IUB, ...)
	
	if( !is.null(check2$message) && !is.null(Ifun) ) return( check2$message )
	
	# no. inequalities
	nineq  = check2$nineq
	xineq0 = check2$xineq0
	if( lpb[ 1 ] == 1 || nineq > 0) lpb[ 2 ] = 1
	
	# check equality
	check3 = .checkeq(pars, Efun, EQ, ...)
	
	if( !is.null(check3$message) ) return( check3$message )
	
	# no. on inequalities (neq)
	neq = check3$neq
	
	# parameter bounds (pb)
	pb  = rbind( cbind(ILB, IUB), cbind(LB, UB) ) 
	
	# check control list
	ctrl  = .checkcontrol( control )
	rho   = ctrl[[ 1 ]]
	maxit = ctrl[[ 2 ]]
	minit = ctrl[[ 3 ]]
	delta = ctrl[[ 4 ]]
	tol   = ctrl[[ 5 ]]
	trace = ctrl[[ 6 ]]
	
	# total constraints (tc)
	tc = nineq + neq
	
	# initialize fn value and inequalities and set to NULL those not needed
	startf = Jfun(pars, ...)
	
	if( nineq>0 ) starti = Ifun(pars,...) else starti = NULL
	
	if( neq>0 ) starte = Efun(pars,...) - EQ else starte = NULL
	
	j  = startf
	jh = startf
	tt = 0 * .ones(3, 1)
	
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
	h  = diag(np+nineq)
	mu = np
	iteration = 0
	ob = c(startf, starte, starti)
	
	while( iteration < maxit ){
		iteration = iteration + 1
		ctrlv = c(rho, minit, delta, tol, neq, nineq, np, lpb)
		res   = .subnp(p, Jfun, Efun, EQ, Ifun, ILB, IUB, LB, UB, control = ctrlv, 
				yy = l, ob = ob,h = h,l = mu, ...)
		p  = res$p
		l  = res$y
		h  = res$h
		mu = res$l
		temp = p[ (nineq + 1):(nineq + np) ]
		Jval = Jfun(temp, ...)
		tempdf = cbind(temp, Jval)
		
		if( trace ){
			xtemp = round(as.numeric(temp), 4)
			cat( paste( "\nIter: ", iteration, " fn: ", Jval, " Pars: ", sep="" ), xtemp, "\n" )
		}
		
		if( neq>0 ) Eval = Efun(temp, ...)-EQ else Eval = NULL
		
		if( nineq>0 ) Ival = Ifun(temp, ...) else Ival = NULL
		
		ob = c(Jval, Eval, Ival)
		tt[ 1 ] = (j - ob[ 1 ]) / max(abs(ob[ 1 ]), 1)
		j=ob[ 1 ]
		
		if( tc > 0.5 ){
			constraint = ob[ 2:(tc + 1) ]
			
			if( nineq > 0.5 ){
				tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ],
				              pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )
				              
				if( min(tempv) > 0 ) {
					p[ 1:nineq ] = constraint[ (neq + 1):tc ]
				}
				
				constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - p[ 1:nineq ]
			}
			
			tt[ 3 ] = .vnorm(constraint)
			
			if( tt[ 3 ] < 10 * tol ) { 
				rho = 0
				mu  = min(mu, tol)
			}
			
			if( tt[ 3 ] < 5 * tt[ 2 ]) {
				rho=rho/5
			}
			
			if( tt[ 3 ] > 10 * tt[ 2 ]) {
				rho = 5 * max( rho, sqrt(tol) )
			}
			
			if( max( c( tol + tt[ 1 ], tt[ 2 ] - tt[ 3 ] ) ) <= 0 ) { 
				l = 0 * l
				h = diag( diag ( h ) )
			}

			tt[ 2 ] = tt[ 3 ]
		}
		
		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			maxit = iteration
		}
		
		jh = c(jh, j)
	}
	
	if( nineq > 0.5 ) {
		xineq0 = p[ 1:nineq ]
	}
	
	p = p[ (nineq + 1):(nineq + np) ]
	
	if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
		convergence = 0
		cat( paste( "\nSOLNP--> Completed in ", iteration, " iterations\n", sep="" ) )
	} else{
		convergence = 1
		cat( paste( "\nSOLNP--> Exiting after maximum number of iterations\n",
						"Tolerance not achieved\n", sep="" ) )
	}
	# end timer
	toc = Sys.time() - tic
	ans = list(par = p, convergence = convergence, values = jh, lagrange = l, hessian = h, xine = xineq0, elapsed = toc)
	return( ans ) 
}