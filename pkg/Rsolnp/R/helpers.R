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


.eps=.Machine$double.eps

.eqfun=function(Efun, pars, EQ,...)
{
	Efun(pars,...)-EQ
}

.subnpmsg<-function(m){
	g1=c("SOLNP-->")
	m1=paste("\n",g1,"Redundant constraints were found. Poor\n",
			g1,"intermediate results may result.Suggest that you\n",
			g1,"remove redundant constraints and re-OPTIMIZE\n",sep="")
	m2=paste("\n",g1,"The linearized problem has no feasible\n",
			g1,"solution.  The problem may not be feasible.\n",sep="")
	m3=paste("\n",g1,"Minor optimization routine did not converge in the \n",
			g1,"specified number of minor iterations.  You may need\n",
			g1,"to increase the number of minor iterations.        \n",sep="")
	ans=switch(m,
			m1=m1,
			m2=m2,
			m3=m3)
	cat(ans)
}
.checkpars<-function(pars,LB,UB){
	np=length(pars)
	message=NULL
	if((!is.null(LB) && np!=length(LB)) || (!is.null(UB) && np!=length(UB))){
		message=paste("lower or upper bound length nor equal to parameter length...exiting",sep="")
	}
	if(!is.null(LB) && any(pars<LB)) message=rbind(message,paste("starting parameters less than lower bounds",sep=""))
	if(!is.null(UB) && any(pars>UB)) message=rbind(message,paste("starting parameters greater than lower bounds",sep=""))
	if(!is.null(LB) && !is.null(UB) && any(LB>UB)){
		message=rbind(message,paste("lower bounds greater than upper bounds!",sep=""))
	}
	if((!is.null(LB) && any(is.na(LB))) || (!is.null(UB) && any(is.na(UB)))){
		message=rbind(message,paste("NAs in lower or upper bound values",sep=""))
	}
	if(any(is.na(pars))) message=rbind(message,paste("NAs in parameter values",sep=""))
	return(message)
}

.checkineq<-function(pars, Ifun,ILB,IUB,...){
	message=NULL
	nineq=0
	xineq0=NULL
	if(!is.function(Ifun)){
		message=rbind(message,paste("Ifun does not appear to be a function",sep=""))
	} else{
		if(!is.null(ILB) && !is.null(IUB) && length(IUB)==length(ILB)){
			testn=length(Ifun(pars,...))
			if(testn!=length(IUB)){
				message=rbind(message,paste("inequality function returns vector of different length
										than given bounds",sep=""))
			} else{
				if(any(ILB>IUB)){
					message=rbind(message,paste("lower inequality bounds greater than upper inequality bounds!",sep=""))
				} else{
					nineq=length(IUB)
					xineq0=(ILB+IUB)/2
				}}
		} else if(!is.null(ILB) && !is.null(IUB) && length(IUB)!=length(ILB)){
			message=rbind(message,paste("length of upper and lower inequality bounds do not match",sep=""))
		} else if(is.null(ILB) || is.null(IUB)){
			message=rbind(message,paste("both upper and lower inequality bounds must be supplied",sep=""))
		}
	}
	return(list(message=message,nineq=nineq,xineq0=xineq0))
}

.checkeq<-function(pars, Efun, EQ,...)
{
	message=NULL
	neq=0
	if(!is.null(Efun)){
		testn=length(Efun(pars,...))
		if(testn!=length(EQ)){
			message=rbind(message, paste("equality function returns vector of different 
									length than given bounds",sep=""))
		} else{
			neq=length(EQ)
		}
	}
	return(list(message=message,neq=neq))
}

.checkcontrol<-function(control){
	# parameters check is now case independent
	ans = list()
	params=unlist(control)
	npar = tolower(names(unlist(control)))
	names(params) = npar
	if(any(substr(npar, 1, 3)=="rho")) ans$rho = as.numeric(params["rho"]) else ans$rho=1
	if(any(substr(npar, 1, 5)=="majit")) ans$majit = as.numeric(params["majit"]) else ans$majit=50
	if(any(substr(npar, 1, 5)=="minit")) ans$minit = as.numeric(params["minit"]) else ans$minit=50
	if(any(substr(npar, 1, 5)=="delta")) ans$delta = as.numeric(params["delta"]) else ans$delta=1.0e-8
	if(any(substr(npar, 1, 3)=="tol")) ans$tol = as.numeric(params["tol"]) else ans$tol=1.0e-6
	if(any(substr(npar, 1, 5)=="trace")) ans$trace = as.numeric(params["trace"]) else ans$trace=0
	return(ans)
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

.vnorm<-function(x){
	sum((x)^2)^(1/2)
}

.solvecond<-function(x)
{
	z=svd(x)$d
	if(any(z==0)) ret=Inf else ret=max(z)/min(z)
	return(ret)
}
