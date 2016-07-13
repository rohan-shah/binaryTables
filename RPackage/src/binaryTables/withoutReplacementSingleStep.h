#ifndef WITHOUT_REPLACEMENT_SINGLE_STEP_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#define WITHOUT_REPLACEMENT_SINGLE_STEP_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace binaryTables
{
	SEXP withoutReplacementSingleStep(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed);
}
#endif
