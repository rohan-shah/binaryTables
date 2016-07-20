#ifndef WITHOUT_REPLACEMENT_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#define WITHOUT_REPLACEMENT_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace binaryTables
{
	SEXP withoutReplacement(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed, SEXP keepTables);
}
#endif
