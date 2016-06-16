#include "conditionalPoissonImpl.h"
#include "GayleRyserTest.h"
#include <boost/math/special_functions/binomial.hpp>
namespace binaryTables
{
	void conditionalPoisson(conditionalPoissonArgs& args)
	{
		std::size_t n = args.n;
		problem& problemObj = args.problemObj;
		const std::vector<int>& rowSums = problemObj.getRowSums();
		const std::vector<int>& columnSums = problemObj.getColumnSums();
		sampling::conditionalPoissonDraftingArgs& samplingArgs = args.samplingArgs;
		std::vector<mpfr_class>& samplingWeights = samplingArgs.weights;

		std::size_t nColumns = columnSums.size(), nRows = rowSums.size();
		GayleRyserTestWorking working(true);

		std::vector<int>& indices = samplingArgs.indices;

		std::vector<int> currentRowSums(nRows);
		std::vector<int> currentColumnSums(nColumns);
		args.estimate = 0;
		for(std::size_t i = 0; i < n; i++)
		{
			mpfr_class density = 1;
			std::copy(rowSums.begin(), rowSums.end(), currentRowSums.begin());
			std::copy(columnSums.begin(), columnSums.end(), currentColumnSums.begin());
			for(std::size_t column = 0; column < nColumns; column++)
			{
				samplingArgs.n = currentColumnSums[column];
				samplingWeights.clear();
				for(std::size_t row = 0; row < nRows; row++)
				{
					samplingWeights.push_back(currentRowSums[row]);
				}
				conditionalPoissonDrafting(samplingArgs, args.randomSource);
				std::sort(indices.begin(), indices.end());
				int deterministicCount = 0;
				for(int j = 0; j < currentColumnSums[column]; j++)
				{
					if(currentRowSums[indices[j]] != (int)nColumns - (int)column && samplingArgs.deterministicInclusion[indices[j]])
					{
						throw std::runtime_error("No units can be deterministically selected except for the last column");
					}
					if(!samplingArgs.deterministicInclusion[indices[j]])
					{
						density *= samplingArgs.expExponentialParameters[indices[j]];
					}
					else deterministicCount++;
					currentRowSums[indices[j]]--;
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
notValid:
			;
		}
		args.estimate /= n;
		//args.estimate /= boost::multiprecision(mpfr_class(2), nRows*nColumns);
	}
}
