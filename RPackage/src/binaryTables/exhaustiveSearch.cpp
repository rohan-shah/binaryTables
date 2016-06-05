#include "exhaustiveSearch.h"
#include "exhaustiveSearchImpl.h"
namespace binaryTables
{
	SEXP exhaustiveSearch(SEXP rowSums_sexp, SEXP columnSums_sexp)
	{
	BEGIN_RCPP
		std::vector<int> rowSums = Rcpp::as<std::vector<int> >(rowSums_sexp);
		std::vector<int> columnSums = Rcpp::as<std::vector<int> >(columnSums_sexp);
		problem problemObj(rowSums, columnSums);
		return Rcpp::wrap(exhaustiveSearch(problemObj));
	END_RCPP
	}
}
