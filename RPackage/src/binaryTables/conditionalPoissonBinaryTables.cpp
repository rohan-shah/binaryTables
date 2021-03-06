#include "conditionalPoissonBinaryTables.h"
#include "conditionalPoissonImpl.h"
namespace binaryTables
{
	SEXP conditionalPoissonR(SEXP rowSums_sexp, SEXP columnSums_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP keepTables_sexp)
	{
	BEGIN_RCPP
		std::vector<int> rowSums = Rcpp::as<std::vector<int> >(rowSums_sexp);
		std::vector<int> columnSums = Rcpp::as<std::vector<int> >(columnSums_sexp);
		int nRows = (int)rowSums.size(), nColumns = (int)columnSums.size();
		bool keepTables = Rcpp::as<bool>(keepTables_sexp);
		problem problemObj(rowSums, columnSums);
		conditionalPoissonArgs args(problemObj);

		int n = Rcpp::as<int>(n_sexp);
		int seed = Rcpp::as<int>(seed_sexp);
		if(n < 1)
		{
			throw std::runtime_error("Input n must be a positive integer");
		}
		args.n = (std::size_t)n;
		args.randomSource.seed(seed);
		args.keepTables = keepTables;
		conditionalPoissonImpl(args);
		
		Rcpp::List tables = Rcpp::List(0);
		Rcpp::CharacterVector tableWeights(0);
		if(keepTables)
		{
			tables = Rcpp::List(args.tables.size()/(nRows * nColumns));
			tableWeights = Rcpp::CharacterVector(tables.size());
			for(int i = 0; i < (int)args.tables.size()/(nRows * nColumns); i++)
			{
				Rcpp::IntegerMatrix table(nRows, nColumns);
				std::copy(args.tables.begin() + i * nRows * nColumns, args.tables.begin() + (i + 1) * nRows * nColumns, table.begin());
				tables(i) = table;
				tableWeights(i) = args.tableWeights[i].str();
			}
		}
		Rcpp::List retVal = Rcpp::List::create(Rcpp::Named("estimate") = args.estimate.str(), Rcpp::Named("varianceEstimate") = args.varianceEstimate.str(), Rcpp::Named("tables") = tables, Rcpp::Named("tableWeights") = tableWeights);
		return retVal;
	END_RCPP
	}
}
