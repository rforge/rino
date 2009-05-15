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

# Based on the original subnp for matlab by Yinyu Ye
# http://www.stanford.edu/~yyye/Col.html
.subnp=function(pars, Jfun, Efun=NULL, EQ=NULL, Ifun=NULL, ILB=NULL, IUB=NULL, LB=NULL,
		UB=NULL, control, yy, ob, h, l, ...)
{
	rho=control[1]
	maxit=control[2]
	delta=control[3]
	tol=control[4]
	neq=control[5]
	nineq=control[6]
	np=control[7]
	lpb=control[8:9]
	ch=1
	alp=c(0,0,0)
	nc=neq+nineq
	npic=np+nineq
	p0=pars
	pb=rbind(cbind(ILB,IUB),cbind(LB,UB))
	sob=numeric()
	ptt=matrix()
	sc=numeric()
	# make the scale for the cost, the equality constraints, the inequality
	# constraints, and the parameters
	#
	if( neq>0 ){
		scale=c(ob[1],.ones(neq,1)*max(abs(ob[2:(neq+1)])))
	} else{
		scale=1
	}
	if(lpb[2]==0){
		scale=c(scale,p0)
	} else{
		scale=c(scale,rep(1,length=length(p0)))
	}
	scale=apply(matrix(scale,ncol=1),1,FUN=function(x) min(max(abs(x),tol),1/tol))
	
	# scale the cost, the equality constraints, the inequality constraints, 
	# the parameters (inequality parameters AND actual parameters), 
	# and the parameter bounds if there are any
	# Also make sure the parameters are no larger than (1-tol) times their bounds
	#
	ob=ob/scale[1:(nc+1)]
	p0=p0/scale[(neq+2):(nc+np+1)]
	if(lpb[2]==1){
		if(lpb[1]==0){
			mm=nineq
		} else{
			mm=npic
		}
		pb=pb/cbind(scale[(neq+2):(neq+mm+1)],scale[(neq+2):(neq+mm+1)])
	}
	# scale the lagrange multipliers and the Hessian
	if(nc>0){
		yy=scale[2:(nc+1)]*yy/scale[1]
	}
	h=h*(scale[(neq+2):(nc+np+1)]%*%t(scale[(neq+2):(nc+np+1)]))/scale[1]
	j=ob[1]
	
	# problem here
	if(nineq>0){
		if(neq==0){
			a=cbind(-diag(nineq),matrix(0,ncol=np,nrow=nineq))
		} else{
			a=rbind(cbind(0*.ones(neq,nineq),matrix(0,ncol=np,nrow=neq)),cbind(-diag(nineq),matrix(0,ncol=np,nrow=nineq)))
		}
	} else{
		a=.zeros(neq,np)
	}
	
	g=0*.ones(npic,1)
	#
	if(nc>0){
		constraint=ob[2:(nc+1)]
		for(i in 1:np){
			p0[nineq+i]=p0[nineq+i]+delta
			Jval=Jfun(p0[(nineq+1):npic]*scale[(nc+2):(nc+np+1)],...)
			if(neq>0) Eval=Efun(p0[(nineq+1):npic]*scale[(nc+2):(nc+np+1)],...)-EQ else Eval=NULL
			if(nineq>0) Ival=Ifun(p0[(nineq+1):npic]*scale[(nc+2):(nc+np+1)],...) else Ival=NULL
			ob=c(Jval,Eval,Ival)
			ob=ob/scale[1:(nc+1)]
			g[nineq+i]=(ob[1]-j)/delta
			a[,nineq+i]=(ob[2:(nc+1)]-constraint)/delta
			p0[nineq+i]=p0[nineq+i]-delta
		}
		if(nineq>0){
			constraint[(neq+1):(neq+nineq)]=constraint[(neq+1):(neq+nineq)]-p0[1:nineq]
		}
		# need to create custom messages
		if(.solvecond(a)>1/.eps){
			.subnpmsg("m1")
		}
		b=a%*%p0-constraint
		#if nc>0.5,
		ch=-1
		alp[1]=tol-max(abs(constraint))
		if(alp[1]<=0){
			ch=1
			if(lpb[2]==0){
				p0=p0-t(a)%*%solve((a%*%t(a)),constraint)	
				alp[1]=1
			}
		}
		if(alp[1]<=0){
			p0[npic+1]=1
			a=cbind(a, -constraint)
			cx=cbind(.zeros(1,npic), 1)
			dx=.ones(npic+1,1)
			go=1 
			minit=0
			while(go>=tol){
				minit=minit+1
				gap=cbind(p0[1:mm]-pb[,1],pb[,2]-p0[1:mm])
				gap=t(apply(gap,1,FUN=function(x) sort(x)))
				dx[1:mm]=gap[,1]
				dx[npic+1]=p0[npic+1]
				if(lpb[1]==0){
					dx[(mm+1):npic]=max(c(dx[1:mm],100))*.ones(npic-mm,1)
				}
				y=qr.solve(t(a%*%diag(as.numeric(dx))),(dx*t(cx)))
				v=dx*(dx*(t(cx)-t(a)%*%y))
				if(v[npic+1]>0){
					z=p0[npic+1]/v[npic+1]
					for(i in 1:mm){
						if(v[i]<0){
							z=min(z,-(pb[i,2]-p0[i])/v[i])
						} else if(v[i]>0){ 
							z=min(z,(p0[i]-pb[i,1])/v[i]) 
						}
					}
					if(z>=p0[npic+1]/v[npic+1]){
						p0=p0-z*v
					} else{
						p0=p0-0.9*z*v 
					}
					go=p0[npic+1]
					if(minit >= 10){
						go=0 
					}
				} else{
					go=0
					minit=10
				}
			}
			if(minit>=10){
				.subnpmsg("m2")
			}
			a=matrix(a[,1:npic],ncol=npic)
			b=a%*%p0[1:npic]
		}
	}
	#
	p=p0[1:npic]
	y=0 
	if(ch>0){
		temppars=p[(nineq+1):npic]*scale[(nc+2):(nc+np+1)]
		startf=Jfun(temppars,...)
		if(nineq>0) starti=Ifun(temppars,...) else starti=NULL
		if(neq > 0) starte=Efun(temppars,...)-EQ else starte=NULL
		ob=c(startf,starte,starti)/scale[1:(nc+1)]
	}						
	j=ob[1]
	#
	if(nineq>0){
		ob[(neq+2):(nc+1)]=ob[(neq+2):(nc+1)]-p[1:nineq]
	}							
	if(nc>0){
		ob[2:(nc+1)]=ob[2:(nc+1)]-a%*%p+b
		j=ob[1]-t(yy)%*%ob[2:(nc+1)]+rho*.vnorm(ob[2:(nc+1)])^2
	}	
	minit=0
	while(minit<maxit){
		minit=minit+1
		if(ch>0){
			for(i in 1:np){
				p[nineq+i]=p[nineq+i]+delta
				temppars=p[(nineq+1):npic]*scale[(nc+2):(nc+np+1)]
				startf=Jfun(temppars,...)
				if(nineq>0) starti=Ifun(temppars,...) else starti=NULL
				if(neq > 0) starte=Efun(temppars,...)-EQ else starte=NULL
				obm=c(startf,starte,starti)/scale[1:(nc+1)]
				if(nineq>0){
					obm[(neq+2):(nc+1)]=obm[(neq+2):(nc+1)]-p[1:nineq]
				}
				if(nc>0){
					obm[2:(nc+1)]=obm[2:(nc+1)]-a%*%p+b
					obm=obm[1]-t(yy)%*%obm[2:(nc+1)]+rho*.vnorm(obm[2:(nc+1)])^2
				}
				g[nineq+i]=(obm-j)/delta
				p[nineq+i]=p[nineq+i]-delta
			}
			if(nineq>0){
				g[1:nineq]=0*yy[(neq+1):nc]	
			}
		}
		# problem here
		if(minit>1){
			yg=g-yg
			sx=p-sx
			sc[1]=t(sx)%*%h%*%sx
			sc[2]=t(sx)%*%yg
			if((sc[1]*sc[2])>0){
				sx=h%*%sx
				h=h-(sx%*%t(sx))/sc[1]+(yg%*%t(yg))/sc[2]
			}
		}
		dx=0.01*.ones(npic,1)
		if(lpb[2]>0.5){
			gap=cbind(p[1:mm]-pb[,1], pb[,2]-p[1:mm])
			gap=t(apply(gap,1,FUN=function(x) sort(x)))
			gap=gap[,1]+sqrt(.eps)*.ones(mm,1)
			dx[1:mm,1]=.ones(mm,1)/gap
			if(lpb[1]<=0){
				dx[(mm+1):npic,1]=min(c(dx[1:mm,1],0.01))*.ones(npic-mm,1)
			}
		}
		go=-1
		##########
		l=l/10
		while(go<=0){
			cz=chol(h+l*diag(as.numeric(dx*dx)))
			cz=solve(cz)
			yg=t(cz)%*%g
			if(nc==0){
				u=-cz%*%yg
			} else{
				y=qr.solve(t(cz)%*%t(a),yg)
				u=-cz%*%(yg-(t(cz)%*%t(a))%*%y)
			}
			p0=u[1:npic]+p
			if(lpb[2]==0){
				go=1
			} else{
				go=min(c(p0[1:mm]-pb[,1],pb[,2]-p0[1:mm]))
				l=3*l
			}
		}
		alp[1]=0
		ob1=ob
		ob2=ob1
		sob[1]=j
		sob[2]=j
		ptt=cbind(p, p)
		alp[3]=1.0
		ptt=cbind(ptt,p0)
		temppars=ptt[(nineq+1):npic,3]*scale[(nc+2):(nc+np+1)]
		startf=Jfun(temppars,...)
		if(nineq>0) starti=Ifun(temppars,...) else starti=NULL
		if(neq > 0) starte=Efun(temppars,...)-EQ else starte=NULL
		ob3=c(startf,starte,starti)/scale[1:(nc+1)]
		sob[3]=ob3[1]
		if(nineq>0){
			ob3[(neq+2):(nc+1)]=ob3[(neq+2):(nc+1)]-ptt[1:nineq,3]
		}
		if(nc>0){
			ob3[2:(nc+1)]=ob3[2:(nc+1)]-a%*%ptt[,3]+b
			sob[3]=ob3[1]-t(yy)%*%ob3[2:(nc+1)]+rho*.vnorm(ob3[2:(nc+1)])^2
		}
		go=1
		while(go>tol){
			alp[2]=(alp[1]+alp[3])/2
			ptt[,2]=(1-alp[2])*p+alp[2]*p0
			temppars=ptt[(nineq+1):npic,2]*scale[(nc+2):(nc+np+1)]
			startf=Jfun(temppars,...)
			if(nineq>0) starti=Ifun(temppars,...) else starti=NULL
			if(neq > 0) starte=Efun(temppars,...)-EQ else starte=NULL
			ob2=c(startf,starte,starti)/scale[1:(nc+1)]
			sob[2]=ob2[1]
			if(nineq>0){
				ob2[(neq+2):(nc+1)]=ob2[(neq+2):(nc+1)]-ptt[1:nineq,2]
			}
			if(nc>0){
				ob2[2:(nc+1)]=ob2[2:(nc+1)]-a%*%ptt[,2]+b
				sob[2]=ob2[1]-t(yy)%*%ob2[2:(nc+1)]+rho*.vnorm(ob2[2:(nc+1)])^2
			}
			obm=max(sob)
			if(obm<j){
				obn=min(sob)
				go=tol*(obm-obn)/(j-obm)
			}
			condif1=sob[2]>=sob[1]
			condif2=sob[1]<=sob[3] && sob[2]<sob[1]
			condif3=sob[2]<sob[1] && sob[1]>sob[3]
			if(condif1){
				sob[3]=sob[2]
				ob3=ob2
				alp[3]=alp[2]
				ptt[,3]=ptt[,2]
			}
			if(condif2){
				sob[3]=sob[2]
				ob3=ob2
				alp[3]=alp[2]
				ptt[,3]=ptt[,2]
			}
			if(condif3){
				sob[1]=sob[2]
				ob1=ob2
				alp[1]=alp[2]
				ptt[,1]=ptt[,2]
			}
			if(go>=tol){
				go=alp[3]-alp[1]
			}
		}
		sx=p
		yg=g
		ch=1
		obn=min(sob)
		if(j<=obn){
			maxit=minit
		}
		reduce=(j-obn)/(1+abs(j))
		if(reduce<tol){
			maxit=minit
		}
		condif1=sob[1]<sob[2]
		condif2=sob[3]<sob[2] && sob[1]>=sob[2]
		condif3=sob[1]>=sob[2] && sob[3]>=sob[2]
		if(condif1){
			j=sob[1]
			p=ptt[,1]
			ob=ob1
		}
		if(condif2){
			j=sob[3]
			p=ptt[,3]
			ob=ob3
		}
		if(condif3){
			j=sob[2]
			p=ptt[,2]
			ob=ob2
		}
	}
	p=p*scale[(neq+2):(nc+np+1)]  # unscale the parameter vector
	if(nc>0){
		y=scale[1]*y/scale[2:(nc+1)] # unscale the lagrange multipliers
	}
	h=scale[1]*h/(scale[(neq+2):(nc+np+1)]%*%t(scale[(neq+2):(nc+np+1)]))
	if(reduce>tol){
		.subnpmsg("m3")
	}
	ans=list(p=p,y=y,h=h,l=l)
	return(ans)
}