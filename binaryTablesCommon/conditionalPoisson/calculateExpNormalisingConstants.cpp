#include "conditionalPoisson/calculateExpNormalisingConstants.h"
namespace binaryTables
{
	namespace conditionalPoisson
	{
		void calculateExpNormalisingConstants(conditionalPoissonSequentialArgs& args)
		{
			int n = args.n;
			int nUnits = (int)args.zeroWeights.size();
			std::vector<bool>& zeroWeights = args.zeroWeights;
			std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
			std::vector<mpfr_class>& expExponentialParameters = args.expExponentialParameters;
			std::vector<mpfr_class>& exponentialParameters = args.exponentialParameters;
			boost::numeric::ublas::matrix<mpfr_class>& expNormalisingConstant = args.expNormalisingConstant;
			int nZeroWeights = 0, nDeterministic = 0;;
			for(int i = 0; i < nUnits; i++)
			{
				if(zeroWeights[i]) nZeroWeights++;
				if(deterministicInclusion[i]) nDeterministic++;
			}
			//We start by computing the normalising constants. First index is k, second is z. All indices in this loop are 1 indexed
			expNormalisingConstant.resize(nUnits - nZeroWeights - nDeterministic, n - nDeterministic);
			//This will skip over the *ignored* units (the ones that were deterministically included)
			int k = (int)expExponentialParameters.size();
			for(int unitIndex = nUnits - nDeterministic - nZeroWeights; unitIndex >= 1; unitIndex--)
			{
				while(zeroWeights[k-1] || deterministicInclusion[k-1]) k--;
				for(int z = 1; z <= std::min(n - nDeterministic, nUnits - nZeroWeights - nDeterministic - unitIndex+1); z++)
				{
					if(z == 1)
					{
						mpfr_class sum = 0;
						for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
						{
							if(!zeroWeights[unitIndex2-1] && !deterministicInclusion[unitIndex2-1]) sum += expExponentialParameters[unitIndex2-1];
						}
						expNormalisingConstant(unitIndex-1, z-1) = sum;
#ifndef NDEBUG
						assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
					}
					else if(z == nUnits - nZeroWeights - nDeterministic - unitIndex + 1)
					{
						mpfr_class sum = 0;
						for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
						{
							if(!zeroWeights[unitIndex2-1] && !deterministicInclusion[unitIndex2-1]) sum += exponentialParameters[unitIndex2-1];
						}
						expNormalisingConstant(unitIndex-1, z-1) = exp(sum);
#ifndef NDEBUG
						assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
					}
					else
					{
						expNormalisingConstant(unitIndex-1, z-1) = expExponentialParameters[k-1] * expNormalisingConstant(unitIndex, z-2) + expNormalisingConstant(unitIndex, z-1);
#ifndef NDEBUG
						assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
					}
				}
				k--;
			}
		}
	}
}
