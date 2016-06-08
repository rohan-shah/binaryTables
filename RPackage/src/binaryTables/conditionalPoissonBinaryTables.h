#ifndef CONDITIONAL_POISSON_BINARY_TABLES_R_PACKAGE
#define CONDITIONAL_POISSON_BINARY_TABLES_R_PACKAGE
#include <Rcpp.h>
namespace binaryTables
{
	SEXP conditionalPoisson(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed);
}
#endif
