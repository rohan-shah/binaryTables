#include "conditionalPoissonBootstrapBinaryTables.h"
#include "conditionalPoissonBootstrapImpl.h"
namespace binaryTables
{
	SEXP conditionalPoissonBootstrap(SEXP rowSums_sexp, SEXP columnSums_sexp, SEXP n_sexp, SEXP seed_sexp)
	{
	BEGIN_RCPP
		std::vector<int> rowSums = Rcpp::as<std::vector<int> >(rowSums_sexp);
		std::vector<int> columnSums = Rcpp::as<std::vector<int> >(columnSums_sexp);
		problem problemObj(rowSums, columnSums);
		conditionalPoissonBootstrapArgs args(problemObj);

		int n = Rcpp::as<int>(n_sexp);
		int seed = Rcpp::as<int>(seed_sexp);
		if(n < 1)
		{
			throw std::runtime_error("Input n must be a positive integer");
		}
		args.n = (std::size_t)n;
		args.randomSource.seed(seed);
		conditionalPoissonBootstrap(args);
		
		return Rcpp::wrap(args.estimate.str());
	END_RCPP
	}
}
