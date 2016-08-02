#ifndef CONDITIONAL_POISSON_IMPL_HEADER_GUARD
#define CONDITIONAL_POISSON_IMPL_HEADER_GUARD
#include "problem.h"
#include "conditionalPoisson/conditionalPoissonSequential.h"
#include "includeMPFRBinaryTables.h"
#include <boost/random/mersenne_twister.hpp>
namespace binaryTables
{
	struct conditionalPoissonArgs
	{
	public:
		conditionalPoissonArgs(problem& problemObj)
			: problemObj(problemObj), keepTables(false)
		{}
		std::size_t n;
		problem& problemObj;
		mpfr_class estimate, varianceEstimate;
		conditionalPoisson::conditionalPoissonSequentialArgs samplingArgs;
		boost::mt19937 randomSource;
		bool keepTables;
		std::vector<bool> tables;
		std::vector<mpfr_class> tableWeights;
	};
	void conditionalPoissonImpl(conditionalPoissonArgs& args);
}
#endif
