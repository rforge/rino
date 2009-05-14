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
	if(!is.null(control$RHO) && is.numeric(control$RHO)) RHO=control$RHO else RHO=1
	if(!is.null(control$MAJIT) && is.numeric(control$MAJIT)) MAJIT=control$MAJIT else MAJIT=10
	if(!is.null(control$MINIT) && is.numeric(control$MINIT)) MINIT=control$MINIT else MINIT=10
	if(!is.null(control$DELTA) && is.numeric(control$DELTA)) DELTA=control$DELTA else DELTA=1.0e-10
	if(!is.null(control$TOL) && is.numeric(control$TOL)) TOL=control$TOL else TOL=1.0e-14
	if(!is.null(control$TRACE) && (control$TRACE==1 || control$TRACE==0)) TRACE=control$TRACE else TRACE=0
	ans=list(RHO=RHO, MAJIT=MAJIT, MINIT=MINIT, DELTA=DELTA, TOL=TOL, TRACE=TRACE)
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
