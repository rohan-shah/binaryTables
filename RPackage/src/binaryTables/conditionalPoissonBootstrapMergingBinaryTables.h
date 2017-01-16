#ifndef CONDITIONAL_POISSON_BOOTSTRAP_MERGING_BINARY_TABLES_R_PACKAGE
#define CONDITIONAL_POISSON_BOOTSTRAP_MERGING_BINARY_TABLES_R_PACKAGE
#include <Rcpp.h>
namespace binaryTables
{
	SEXP conditionalPoissonBootstrapMerging(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed, SEXP mergeFrequency);
}
#endif
