#ifndef WITHOUT_REPLACEMENT_MERGING_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#define WITHOUT_REPLACEMENT_MERGING_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace binaryTables
{
	SEXP withoutReplacementMerging(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed, SEXP mergeFrequency);
}
#endif
