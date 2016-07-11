#include "withoutReplacement.h"
#include "withoutReplacementImpl.h"
namespace binaryTables
{
	SEXP withoutReplacement(SEXP rowSums_sexp, SEXP columnSums_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP keepTables_sexp)
	{
	BEGIN_RCPP
		std::vector<int> rowSums = Rcpp::as<std::vector<int> >(rowSums_sexp);
		std::vector<int> columnSums = Rcpp::as<std::vector<int> >(columnSums_sexp);
		int nRows = (int)rowSums.size(), nColumns = (int)columnSums.size();
		bool keepTables = Rcpp::as<bool>(keepTables_sexp);

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

		problem problemObj(rowSums, columnSums);
		withoutReplacementArgs args(problemObj);
		args.n = n;
		args.randomSource.seed(seed);
		args.keepTables = keepTables;
		withoutReplacement(args);
		std::string estimate_str = args.estimate.str();
		Rcpp::List retVal(2); 
		if(keepTables)
		{
			Rcpp::List tables(args.tables.size()/(nRows * nColumns));
			for(int i = 0; i < (int)args.tables.size()/(nRows * nColumns); i++)
			{
				Rcpp::IntegerMatrix table(nRows, nColumns);
				std::copy(args.tables.begin() + i * nRows * nColumns, args.tables.begin() + (i + 1) * nRows * nColumns, table.begin());
				tables(i) = table;
			}
			retVal(1) = tables;
		}
		retVal(0) = estimate_str;
		std::vector<std::string> retValNames = {"estimate", "tables"};
		retVal.names() = retValNames;
		return retVal;
	END_RCPP
	}
}
