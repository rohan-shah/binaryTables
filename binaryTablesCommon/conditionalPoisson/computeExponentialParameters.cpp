#include "conditionalPoisson/computeExponentialParameters.h"
namespace binaryTables
{
	namespace conditionalPoisson
	{
		void exponentialPreCompute(computeExponentialData& data, int nColumns)
		{
			data.ratios.resize(nColumns+1, nColumns+1);
			data.logRatios.resize(nColumns+1, nColumns+1);
			//Row is number of columns which are full
			for(int i = 0; i < nColumns+1; i++)
			{
				//Column is number of columns remaining total.
				for(int j = i+1; j < nColumns+1; j++)
				{
					mpfr_class weight = (double)i / (double)j;
					data.ratios(i, j) = weight / (1 - weight);
					data.logRatios(i, j) = boost::multiprecision::log(data.ratios(i, j));
				}
			}
		};
		void computeExponentialParameters(conditionalPoissonSequentialArgs& args, int nColumnsRemaining, std::vector<int>::const_iterator rowSumsStart, std::vector<int>::const_iterator rowSumsEnd)
		{
			int nUnits = (int)std::distance(rowSumsStart, rowSumsEnd);

			args.exponentialParameters.resize(nUnits);
			args.expExponentialParameters.resize(nUnits);
			computeExponentialData& data = args.preComputation;
			mpfr_class sumExponentialParameters = 0;
			int excluded = 0;
			for(int i = 0; i < nUnits; i++)
			{
				if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
				{
					args.expExponentialParameters[i] = data.ratios(*(rowSumsStart + i), nColumnsRemaining);
					args.exponentialParameters[i] = data.logRatios(*(rowSumsStart + i), nColumnsRemaining);
					sumExponentialParameters += args.exponentialParameters[i];
				}
				else
				{
					excluded++;
					args.expExponentialParameters[i] = 0;
					args.exponentialParameters[i] = 0;
				}
			}
			//Rescale so the exponential parameters sum no zero
			mpfr_class toSubtract = sumExponentialParameters / (nUnits - excluded);
			mpfr_class expToSubtract = boost::multiprecision::exp(toSubtract);
			for(int i = 0; i < nUnits; i++)
			{
				if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
				{
					args.exponentialParameters[i] -= toSubtract;
					args.expExponentialParameters[i] /= expToSubtract;
				}
			}
		}
	}
}
