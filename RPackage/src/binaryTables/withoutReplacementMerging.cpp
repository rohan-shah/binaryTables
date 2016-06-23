#include "withoutReplacementMerging.h"
#include "withoutReplacementMergingImpl.h"
namespace binaryTables
{
	SEXP withoutReplacementMerging(SEXP rowSums_sexp, SEXP columnSums_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP mergeFrequency_sexp)
	{
	BEGIN_RCPP
		std::vector<int> rowSums = Rcpp::as<std::vector<int> >(rowSums_sexp);
		std::vector<int> columnSums = Rcpp::as<std::vector<int> >(columnSums_sexp);

		if(std::find(rowSums.begin(), rowSums.end(), 0) != rowSums.end())
		{
			throw std::runtime_error("Input rowSums cannot contain 0");
		}
		if(std::find(columnSums.begin(), columnSums.end(), 0) != columnSums.end())
		{
			throw std::runtime_error("Input columnSums cannot contain 0");
		}

		std::size_t n = (std::size_t)Rcpp::as<int>(n_sexp);
		int seed = Rcpp::as<int>(seed_sexp);
		int mergeFrequency = Rcpp::as<int>(mergeFrequency_sexp);

		problem problemObj(rowSums, columnSums);
		withoutReplacementMergingArgs args(problemObj);
		args.n = n;
		args.randomSource.seed(seed);
		args.mergeFrequency = mergeFrequency;
		withoutReplacementMerging(args);
		std::string estimate_str = args.estimate.str();
		return Rcpp::wrap(estimate_str);
	END_RCPP
	}
}
