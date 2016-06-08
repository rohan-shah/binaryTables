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
		sampling::conditionalPoissonArgs& samplingArgs = args.samplingArgs;
		std::vector<mpfr_class>& inclusionProbabilities = args.inclusionProbablities;
		std::vector<mpfr_class>& samplingWeights = args.samplingWeights;
		std::vector<int> indices;

		std::size_t nColumns = columnSums.size(), nRows = rowSums.size();
		GayleRyserTestWorking working(true);

		samplingArgs.inclusionProbabilities = &inclusionProbabilities;
		samplingArgs.weights = &samplingWeights;
		std::vector<mpfr_class> rescaledWeights;
		samplingArgs.rescaledWeights = &rescaledWeights;
		samplingArgs.indices = &indices;

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
				samplingArgs.calculateInclusionProbabilities = true;
				conditionalPoisson(samplingArgs, args.randomSource);
				std::sort(indices.begin(), indices.end());
				int deterministicCount = 0;
				for(int j = 0; j < currentColumnSums[column]; j++)
				{
					if(column != nColumns - 1 && samplingArgs.deterministicInclusion[indices[j]])
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