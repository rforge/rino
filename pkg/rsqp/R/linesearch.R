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


.constrlinesearch = function(x, p, lambda = NULL, funv, .env, ...)
{
	## Choose parameters
	##
	## eta in the range (0, 0.5)
	## tau in the range (0, 1)
	.sqp_index = get(".sqp_index", envir=.env)
	fun = get(".sqp_fun", envir=.env)
	gradfun = get(".sqp_gradfun", envir=.env)
	eqfun = get(".sqp_eqfun", envir=.env)
	ineqfun = get(".sqp_ineqlbfun", envir=.env)
	eta = 0.1
	tau = 0.7
	deltabar = sqrt (.eps)
	if (is.null(lambda)) mu = 1 / deltabar else mu = 1 / (.rownorm(lambda) + deltabar)
	alpha = 1
	gradv = gradfun(x, ...)
	if(.sqp_index[7]) eqv = eqfun(x, ...) else eqv = NULL
	tmp = .meritfn(x, mu, funv, .env,  ...)
	phixmu = tmp$merit
	funv = tmp$funv
	Dphixmu = t(gradv)%*%p
	ineqlbv = ineqfun(x, ...)
	## only those elements of d corresponding
	## to violated constraints should be included.
	idx = which(ineqlbv < 0)
	if(length(idx) == 0) didx = NULL else didx = ineqlbv[idx]
	if(!is.null(didx) || !is.null(eqv)) {
		xt = -.colnorm(c(eqv, didx)) / mu
		Dphixmu = Dphixmu + xt
	}
	while (1){
		tmp = .meritfn(x + alpha * p, mu, funv = NULL, .env, ...)
		p1 = tmp$merit
		funv = tmp$funv
		p2 = phixmu + eta*alpha%*%Dphixmu
		if (p1 > p2){
			tau_alpha = 0.9 * tau
			alpha = tau_alpha * alpha
		} else{
			break()
		}
	}
	xnew = x + alpha * p
	return(list(xnew = xnew, alpha = alpha, funv = funv))
}