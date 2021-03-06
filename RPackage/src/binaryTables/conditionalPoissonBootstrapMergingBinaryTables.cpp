#include "conditionalPoissonBootstrapMergingBinaryTables.h"
#include "conditionalPoissonBootstrapMergingImpl.h"
namespace binaryTables
{
	SEXP conditionalPoissonBootstrapMerging(SEXP rowSums_sexp, SEXP columnSums_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP mergeFrequency_sexp)
	{
	BEGIN_RCPP
		std::vector<int> rowSums = Rcpp::as<std::vector<int> >(rowSums_sexp);
		std::vector<int> columnSums = Rcpp::as<std::vector<int> >(columnSums_sexp);
		int mergeFrequency = Rcpp::as<int>(mergeFrequency_sexp);
		problem problemObj(rowSums, columnSums);
		conditionalPoissonBootstrapMergingArgs args(problemObj);

		int n = Rcpp::as<int>(n_sexp);
		int seed = Rcpp::as<int>(seed_sexp);
		if(n < 1)
		{
			throw std::runtime_error("Input n must be a positive integer");
		}
		args.n = (std::size_t)n;
		args.randomSource.seed(seed);
		args.mergeFrequency = mergeFrequency;
		conditionalPoissonBootstrapMerging(args);
		
		return Rcpp::wrap(args.estimate.str());
	END_RCPP
	}
}
