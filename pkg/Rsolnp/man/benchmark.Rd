\name{benchmark}
\Rdversion{1.1}
\alias{benchmark}
\title{
The Rsolnp Benchmark Problems Suite.
}
\description{
The function implements a set of benchmark problems against the MINOS solver of 
Murtagh and Saunders.
}
\usage{
benchmark(dat, id = "Powell", print = TRUE)
}
\arguments{
  \item{dat}{the benchmark data to be used.}
  \item{id}{
The name of the benchmark problem. A call to the function \code{\link{benchmarkids}}
will return the available benchmark problems.
}
  \item{print}{
Prints the output to screen with a short descriptions of the benchmark problem.
}
}
\details{
The benchmarks were run on an intel core 2 extreme x9000 cpu with 8GB of memory on windows
vista operating system. The MINOS solver was used via the tomlab interface to matlab. The
solnp benchmarks are available to be run in the unit test section of the package directory.
}
\value{
A data.frame containing the benchmark data.
}
\references{
W.Hock and K.Schittkowski, \emph{Test Examples for Nonlinear Programming Codes}, 
Lecture Notes in Economics and Mathematical Systems. Springer Verlag, 1981.\cr 
Y.Ye,\emph{Interior algorithms for linear, quadratic, and linearly constrained non 
linear programming}, PhD Thesis, Department of EES Stanford University, Stanford CA.\cr
B.A.Murtagh and M.A.Saunders, \emph{MINOS 5.5 User's Guide, Report SOL 83-20R}, Systems 
Optimization Laboratory, Stanford University (revised July 1998).
}
\author{
Alexios Ghalanos and Stefan Theussl\cr
Y.Ye (original matlab version of solnp)
}
\examples{
benchmarkids()
dat <- data("benchmarklist")
powell <- benchmark(get(dat) , id = "Powell", print = TRUE )
}
\keyword{optimize}