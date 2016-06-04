#ifndef CRUDE_MC_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#define CRUDE_MC_BINARY_TABLES_R_PACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace binaryTables
{
	SEXP crudeMC(SEXP rowSums, SEXP columnSums, SEXP n, SEXP seed);
}
#endif
