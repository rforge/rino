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

benchmarkids <- function()
{
return(c("Powell","Wright4","Wright9","Alkyla","Entropy","Box","RosenSuzuki","PortRachev",
"PortKappa"))
}

# We need:
# Optimal
# 1. Function Value
# 2. Parameters
# 3. Convergence code
# 4. (major) Iterations/Function Evaluations
# 5. Time taken

benchmark <- function( dat, id = "Powell", print = TRUE )
{
  if( !any(benchmarkids() == id[ 1L ]) )
    stop( "invalid benchmark id" )
 
  n1 <- which( id[ 1L ] == names(dat$solnp) )
  n2 <- which( id[ 1L ] == names(dat$minos) )
  d1 <- which( id[ 1L ] == names(dat$describe) )
  bench1 <- dat$solnp[[ n1 ]]
  bench2 <- dat$minos[[ n2 ]]
  bt <- .makebenchtable( bench1, bench2 )
  if( print ){
    cat( "\n" )
    cat( paste("Problem :", id, sep = "") )
    cat( "\n" )
    print( bt )
    cat( "\n" )
    cat( dat$describe[[ d1 ]], "\n" )
  }
  return( bt )
}



.makebenchtable <- function( bench1, bench2 )
{
  bt <- data.frame( solnp = rbind(round(bench1$fn, 5L),
                                  round(bench1$iter, 0L),
                                  round(bench1$exitflag, 0L),
                                  round(bench1$elapsed, 3L),
                                  matrix(round(bench1$pars, 5L), ncol = 1L)),
                   minos =  rbind(round(bench2$fn, 5L),
                                  round(bench2$iter, 0L),
                                  round(bench2$exitflag, 0L),
                                  round(bench2$elapsed, 3L),
                                  matrix(round(bench2$pars, 5L), ncol = 1L)) )
  rownames(bt) <- c("funcValue", "majorIter", "exitFlag", "time(sec)",
                    paste("par.", 1L:length(bench1$pars), sep = "") )
  return(bt)
}
