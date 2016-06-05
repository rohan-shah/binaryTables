#ifndef EXHAUSTIVE_SEARCH_R_PACKAGE_HEADER_GUARD
#define EXHAUSTIVE_SEARCH_R_PACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace binaryTables
{
	SEXP exhaustiveSearch(SEXP rowSums, SEXP columnSums);
}
#endif
