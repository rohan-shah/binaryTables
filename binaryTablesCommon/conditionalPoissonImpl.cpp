#include "conditionalPoissonImpl.h"
#include "GayleRyserTest.h"
#include <boost/math/special_functions/binomial.hpp>
#include "conditionalPoisson/computeExponentialParameters.h"
namespace binaryTables
{
	void conditionalPoissonImpl(conditionalPoissonArgs& args)
	{
		std::size_t n = args.n;
		problem& problemObj = args.problemObj;
		const std::vector<int>& rowSums = problemObj.getRowSums();
		const std::vector<int>& columnSums = problemObj.getColumnSums();
		conditionalPoisson::conditionalPoissonSequentialArgs& samplingArgs = args.samplingArgs;

		std::size_t nColumns = columnSums.size(), nRows = rowSums.size();
		conditionalPoisson::exponentialPreCompute(samplingArgs.preComputation, nColumns);
		GayleRyserTestWorking working(false);

		std::vector<int>& indices = samplingArgs.indices;

		if(args.keepTables)
		{
			args.tables.clear();
			args.tables.reserve(n * nRows*nColumns);
		}
		std::vector<bool> currentTable(nRows*nColumns);

		std::vector<int> currentRowSums(nRows);
		std::vector<int> currentColumnSums(nColumns);
		mpfr_class secondMoment = 0;
		args.estimate = 0;
		for(std::size_t i = 0; i < n; i++)
		{
			mpfr_class density = 1;
			std::copy(rowSums.begin(), rowSums.end(), currentRowSums.begin());
			std::copy(columnSums.begin(), columnSums.end(), currentColumnSums.begin());
			std::fill(currentTable.begin(), currentTable.end(), false);
			for(std::size_t column = 0; column < nColumns; column++)
			{
				samplingArgs.n = currentColumnSums[column];
				conditionalPoisson::conditionalPoissonSequential(samplingArgs, args.randomSource, currentRowSums.begin(), currentRowSums.end(), nColumns - column);
				std::sort(indices.begin(), indices.end());
				int deterministicCount = 0;
				for(int j = 0; j < currentColumnSums[column]; j++)
				{
					if(!samplingArgs.deterministicInclusion[indices[j]])
					{
						density *= samplingArgs.expExponentialParameters[indices[j]];
					}
					else deterministicCount++;
					currentRowSums[indices[j]]--;
					currentTable[column*nRows + indices[j]] = true;
				}
				//If all the units are deterministic then we don't need to do anything. 
				if(currentColumnSums[column] != deterministicCount)
				{
					density /= samplingArgs.expNormalisingConstant(0, currentColumnSums[column] - deterministicCount - 1);
				}

				if(!GayleRyserTest(currentColumnSums, currentRowSums, column+1, working))
				{
					goto notValid;
				}
			}
			for(std::size_t rowCounter = 0; rowCounter < nRows; rowCounter++)
			{
				if(currentRowSums[rowCounter] > 0) goto notValid;
			}
			args.estimate += 1/density;
			secondMoment += 1/(density*density);
			if(args.keepTables)
			{
				args.tables.resize(args.tables.size() + nRows * nColumns);
				std::copy(currentTable.begin(), currentTable.end(), args.tables.end() - nRows * nColumns);
			}
notValid:
			;
		}
		secondMoment /= n;
		args.estimate /= n;
		args.varianceEstimate = secondMoment - args.estimate*args.estimate;
	}
}
