#ifndef CONDITIONAL_POISSON_BOOTSTRAP_BINARY_TABLES_R_PACKAGE
#define CONDITIONAL_POISSON_BOOTSTRAP_BINARY_TABLES_R_PACKAGE
#include <Rcpp.h>
namespace binaryTables
{
	SEXP conditionalPoissonBootstrap(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed);
}
#endif
