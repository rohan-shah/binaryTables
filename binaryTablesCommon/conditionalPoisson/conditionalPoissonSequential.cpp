#include "conditionalPoisson/conditionalPoissonSequential.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "conditionalPoisson/computeExponentialParameters.h"
#include "conditionalPoisson/calculateExpNormalisingConstants.h"
namespace binaryTables
{
	namespace conditionalPoisson
	{
		void conditionalPoissonBase(conditionalPoissonSequentialArgs& args, std::vector<int>::const_iterator rowSumsBegin, std::vector<int>::const_iterator rowSumsEnd, int nRemainingColumns, int& nDeterministic, int& nZeroWeights)
		{
			std::vector<int>& indices = args.indices;
			std::vector<bool>& zeroWeights = args.zeroWeights;
			std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
			int nUnits = (int)std::distance(rowSumsBegin, rowSumsEnd);
			indices.clear();
			nZeroWeights = 0;
			nDeterministic = 0;
			
			zeroWeights.resize(nUnits);
			deterministicInclusion.resize(nUnits);
			std::fill(zeroWeights.begin(), zeroWeights.end(), false);
			std::fill(deterministicInclusion.begin(), deterministicInclusion.end(), false);

			std::vector<int> remaining;
			remaining.reserve(nUnits);
			for(int i = 0; i < nUnits; i++)
			{
				if(*(rowSumsBegin + i) == 0) 
				{
					nZeroWeights++;
					zeroWeights[i] = true;
				}
				else if(*(rowSumsBegin + i) == nRemainingColumns)
				{
					nDeterministic++;
					deterministicInclusion[i] = true;
					indices.push_back(i);
				}
				else remaining.push_back(i);
			}
		}
		void conditionalPoissonSequential(conditionalPoissonSequentialArgs& args, boost::mt19937& randomSource, std::vector<int>::const_iterator rowSumsBegin, std::vector<int>::const_iterator rowSumsEnd, int nRemainingColumns)
		{
			std::vector<int>& indices = args.indices;
			int nUnits = (int)std::distance(rowSumsBegin, rowSumsEnd);
			
			int nDeterministic, nZeroWeights;
			conditionalPoissonBase(args, rowSumsBegin, rowSumsEnd, nRemainingColumns, nDeterministic, nZeroWeights);

			computeExponentialParameters(args, nRemainingColumns, rowSumsBegin, rowSumsEnd); 
			calculateExpNormalisingConstants(args);
			if(indices.size() == args.n) return;
			int chosen = 0;
			int skipped = 0;
			for(int i = 0; i < nUnits; i++)
			{
				if(!args.zeroWeights[i] && !args.deterministicInclusion[i])
				{
					double parameter;
					if(args.n - nDeterministic - 1 - chosen == 0)
					{
						parameter = (args.expExponentialParameters[i] / args.expNormalisingConstant(i-skipped, args.n - nDeterministic - chosen - 1)).convert_to<double>();
					}
					else
					{
						parameter = (args.expExponentialParameters[i] * args.expNormalisingConstant(i+1-skipped, args.n - nDeterministic - 1 - chosen - 1) / args.expNormalisingConstant(i-skipped, args.n - nDeterministic - chosen - 1)).convert_to<double>();
					}
#ifndef NDEBUG
					if(parameter > 1) throw std::runtime_error("Internal error");
#endif
					boost::random::bernoulli_distribution<> bernoulli(parameter);
					if(bernoulli(randomSource))
					{
						indices.push_back(i);
						chosen++;
					}
				}
				else skipped++;
				if(chosen == (int)args.n - nDeterministic) break;
			}
		}
	}
}
