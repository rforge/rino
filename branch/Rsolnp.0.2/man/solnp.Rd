\name{solnp}
\Rdversion{1.1}
\alias{solnp}
\title{
Nonlinear optimization using augmented lagrange method.
}
\description{
The solnp function is based on the solver by Yinyu Ye which solves the general nonlinear
programming problem:\cr
minimize		\deqn{f(x)}\cr 
subject to		\deqn{g(x) = 0}\cr
			\deqn{l_h \le h(x) \le u_h}\cr
			\deqn{l_x \le x \le u_x}\cr
where, \eqn{f(x)},\eqn{g(x)} and \eqn{g(x)} are smooth functions.
}
\usage{
solnp(pars, Jfun, Efun = NULL, EQ = NULL, Ifun = NULL, ILB = NULL, IUB = NULL, 
LB = NULL, UB = NULL, control = list(), ...)
}
\arguments{
  \item{pars}{
The starting parameter vector.
}
  \item{Jfun}{
The main function which takes as first argument the parameter vector and returns
a single value.
}
  \item{Efun}{
(Optional) The equality constraint function returning the vector of evaluated equality 
constraints.
}
  \item{EQ}{
(Optioal) The equality constraints.
}
  \item{Ifun}{
(Optional) The inequality constraint function returning the vector of evaluated inequality 
constraints.
}
  \item{ILB}{
(Optional) The lower bound of the inequality constraints.
}
  \item{IUB}{
(Optional) The upper bound of the inequality constraints.
}
  \item{LB}{
(Optional) The lower bound on the parameters.
}
  \item{UB}{
(Optional) The upper bound on the parameters.
}
  \item{control}{
(Optional) The control list of optimization parameters. See below for details.
}
  \item{\dots}{
(Optional) Additional parameters passed to the main, equality or inequality functions
}
}
\details{
The solver belongs to the class of indirect solvers and implements the augmented Lagrange
multiplier method with an SQP interior algorithm.
}
\value{
A list containing the following values:
\item{par}{Optimal Parameters.}
\item{convergence }{Indicates whether the solver has converged (0) or not (1).}
\item{values}{Vector of function values during optimization with last one the
value at the optimal.}
\item{lagrange}{The vector of Lagrange multipliers.}
\item{hessian}{The Hessian at the optimal solution.}
\item{xineq}{The estimated optimal inequality vector of slack variables used for 
transforming the inequality into an equality constraint.}
\item{elapsed}{Time taken to compute solution.}
}
\section{Control}{
\describe{
\item{RHO}{Penalty parameter (default 1) }
\item{MAJIT}{Maximum number of major iterations (default 50)}
\item{MINIT}{Maximum number of minor iterations (default 50)}
\item{DELTA}{Relative step size in forward difference evaluation (default 1.0e-10)}
\item{TOL}{Tolerance on feasibility and optimality (default 1e-14)}
\item{TRACE}{The value of the objective function and the parameters is printed at every
major iteration (default 0)}
}}
\references{
Y.Ye,\emph{Interior algorithms for linear, quadratic, and linearly constrained non 
linear programming}, PhD Thesis, Department of EES Stanford University, Stanford CA.
}
\author{
Alexios Ghalanos and Stefan Theussl\cr
Y.Ye (original matlab version of solnp)
}
\note{
The control parameters TOL and DELTA are key in getting any possibility of successful
convergence, therefore it is suggested that the user change these appropriately to reflect
their problem specification.
}
\examples{
# From the original paper by Y.Ye
# see the unit tests for more....
#---------------------------------------------------------------------------------
# POWEL Problem
fn1=function(x)
{
	exp(x[1]*x[2]*x[3]*x[4]*x[5])
}

Jn1=function(x){
	z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
	z2=x[2]*x[3]-5*x[4]*x[5]
	z3=x[1]*x[1]*x[1]+x[2]*x[2]*x[2]
	return(c(z1,z2,z3))
}


pars1=c(-2,2,2,-1,-1)
ctrl=list(TOL=1e-6,DELTA=1e-8)
powell=solnp(pars1, Jfun=fn1, Efun=Jn1, EQ=c(10,0,-1), control=ctrl)
}
\keyword{optimize}