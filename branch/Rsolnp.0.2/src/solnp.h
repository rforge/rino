#ifndef SOLNP_H
#define SOLNP_H
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string>
#include <vector>
using namespace std;

#include "solfun.h"
#include "newmatap.h"
#include "newmat.h"
using namespace NEWMAT;

extern "C" {
// C versions of
SEXP solnp(SEXP , SEXP , SEXP , SEXP , SEXP, SEXP , SEXP ,SEXP , SEXP , SEXP ,
		SEXP , SEXP , SEXP  );
SEXP solnpeval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
}
#endif
