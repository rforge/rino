## R Test


test_bonmin <- function() {

	x = 1
	out <- .C('R_bonmin_solve', test=as.integer(x), PACKAGE='Rbonmin')
	out
}

