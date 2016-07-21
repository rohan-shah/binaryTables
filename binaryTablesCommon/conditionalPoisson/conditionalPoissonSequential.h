#ifndef CONDITIONAL_POISSON_DRAFTING_BINARY_TABLES_HEADER_GUARD
#define CONDITIONAL_POISSON_DRAFTING_BINARY_TABLES_HEADER_GUARD
#include <vector>
#include "includeMPFRBinaryTables.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace binaryTables
{
	namespace conditionalPoisson
	{
		struct computeExponentialData
		{
			boost::numeric::ublas::matrix<mpfr_class> ratios;
			boost::numeric::ublas::matrix<mpfr_class> logRatios;
		};
		struct conditionalPoissonSequentialArgs
		{
		public:
			conditionalPoissonSequentialArgs()
			{}
			//The vector of selected units
			std::vector<int> indices;
			//Number of units to select
			std::size_t n;

			std::vector<mpfr_class> exponentialParameters;
			std::vector<mpfr_class> expExponentialParameters;
			boost::numeric::ublas::matrix<mpfr_class> expNormalisingConstant;
			std::vector<bool> deterministicInclusion, zeroWeights;

			computeExponentialData preComputation;
		private:
			conditionalPoissonSequentialArgs(const conditionalPoissonSequentialArgs& other);
			conditionalPoissonSequentialArgs& operator=(const conditionalPoissonSequentialArgs& other);
		};
		void conditionalPoissonBase(conditionalPoissonSequentialArgs& args, std::vector<int>::const_iterator rowSumsBegin, std::vector<int>::const_iterator rowSumsEnd, int nRemainingColumns, int& nDeterministic, int& nZeroWeights);
		void conditionalPoissonSequential(conditionalPoissonSequentialArgs& args, boost::mt19937& randomSource, std::vector<int>::const_iterator rowSumsBegin, std::vector<int>::const_iterator rowSumsEnd, int nRemainingColumns);
	}
}
#endif
