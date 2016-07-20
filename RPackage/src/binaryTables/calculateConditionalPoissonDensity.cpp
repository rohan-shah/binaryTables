#include "calculateConditionalPoissonDensity.h"
#include "conditionalPoissonBase.h"
#include "samplingBase.h"
namespace binaryTables
{
	SEXP calculateConditionalPoissonDensity(SEXP binaryTables_sexp)
	{
		Rcpp::List binaryTables = Rcpp::as<Rcpp::List>(binaryTables_sexp);
		std::size_t nTables = binaryTables.size();
		Rcpp::IntegerMatrix firstTable = Rcpp::as<Rcpp::IntegerMatrix>(binaryTables[0]);
		int nRows = firstTable.nrow(), nColumns = firstTable.ncol();

		//Calculate row and column sums
		std::vector<int> rowSums(nRows, 0), columnSums(nColumns, 0);
		for(int row = 0; row < nRows; row++)
		{
			for(int column = 0; column < nColumns; column++)
			{
				if(firstTable(row, column) != 0 && firstTable(row, column) != 1)
				{
					throw std::runtime_error("Tables were not binary");
				}
				if(firstTable(row, column) == 1)
				{
					rowSums[row]++;
					columnSums[column]++;
				}
			}
		}
		std::vector<int> currentRowSums(nRows, 0), currentColumnSums(nColumns, 0);
		for(int tableCounter = 0; tableCounter < (int)nTables; tableCounter++)
		{
			Rcpp::IntegerMatrix currentTable = Rcpp::as<Rcpp::IntegerMatrix>(binaryTables[tableCounter]);
			if((int)currentTable.nrow() != nRows || (int)currentTable.ncol() != nColumns)
			{
				throw std::runtime_error("Tables had different dimensions");
			}
			std::fill(currentRowSums.begin(), currentRowSums.end(), 0);
			std::fill(currentColumnSums.begin(), currentColumnSums.end(), 0);
			for(int row = 0; row < nRows; row++)
			{
				for(int column = 0; column < nColumns; column++)
				{
					if(currentTable(row, column) != 0 && currentTable(row, column) != 1)
					{
					}
					if(currentTable(row, column) == 1)
					{
						currentRowSums[row]++;
						currentColumnSums[column]++;
					}
				}
			}
			if(memcmp(&(currentRowSums[0]), &(rowSums[0]), sizeof(int)*nRows))
			{
				throw std::runtime_error("Row sums were inconsistent");
			}
			if(memcmp(&(currentRowSums[0]), &(rowSums[0]), sizeof(int)*nRows))
			{
				throw std::runtime_error("Column sums were inconsistent");
			}
		}

		sampling::conditionalPoissonArgs samplingArgs;
		samplingArgs.deterministicInclusion.resize(nRows);
		samplingArgs.zeroWeights.resize(nRows);
		std::vector<int>& indices = samplingArgs.indices;
		std::vector<sampling::mpfr_class> densities(nTables);
		std::vector<std::string> densitiesAsStrings(nTables);
		for(int tableCounter = 0; tableCounter < (int)nTables; tableCounter++)
		{
			Rcpp::IntegerMatrix currentTable = Rcpp::as<Rcpp::IntegerMatrix>(binaryTables[tableCounter]);
			std::copy(rowSums.begin(), rowSums.end(), currentRowSums.begin());
			sampling::mpfr_class density = 1;
			for(int column = 0; column < nColumns; column++)
			{
				samplingArgs.n = columnSums[column];
				samplingArgs.weights.clear();
				for(int row = 0; row < nRows; row++)
				{
					samplingArgs.weights.push_back(currentRowSums[row] / (double)(nColumns - column));
				}
				std::fill(samplingArgs.deterministicInclusion.begin(), samplingArgs.deterministicInclusion.end(), false);
				std::fill(samplingArgs.zeroWeights.begin(), samplingArgs.zeroWeights.end(), false);
				indices.clear();
				int nDeterministic = 0, nZeroWeights = 0;
				sampling::samplingBase(samplingArgs.n, indices, samplingArgs.weights, samplingArgs.zeroWeights, samplingArgs.deterministicInclusion, nDeterministic, nZeroWeights);
				computeExponentialParameters(samplingArgs);
				calculateExpNormalisingConstants(samplingArgs);
				int deterministicCount = 0;
				for(int row = 0; row < nRows; row++)
				{
					if(currentTable(row, column))
					{
						if(!samplingArgs.deterministicInclusion[row])
						{
							density *= samplingArgs.expExponentialParameters[row];
						}
						else deterministicCount++;
						currentRowSums[row]--;
					}
				}
				if(columnSums[column] != deterministicCount)
				{
					density /= samplingArgs.expNormalisingConstant(0, columnSums[column] - deterministicCount - 1);
				}
			}
			densities[tableCounter] = density;
			densitiesAsStrings[tableCounter] = density.str();
		}
		return Rcpp::wrap(densitiesAsStrings);
	}
}
