#ifndef CONDITIONAL_POISSON_IMPL_HEADER_GUARD
#define CONDITIONAL_POISSON_IMPL_HEADER_GUARD
#include "problem.h"
#include "conditionalPoissonDrafting.h"
#include "conditionalPoissonSequential.h"
#include "includeMPFRBinaryTables.h"
#include <boost/random/mersenne_twister.hpp>
namespace binaryTables
{
	struct conditionalPoissonArgs
	{
	public:
		conditionalPoissonArgs(problem& problemObj)
			: problemObj(problemObj)
		{}
		std::size_t n;
		problem& problemObj;
		mpfr_class estimate, varianceEstimate;
		sampling::conditionalPoissonDraftingArgs samplingArgs;
		boost::mt19937 randomSource;
		bool keepTables;
		std::vector<bool> tables;
	};
	void conditionalPoisson(conditionalPoissonArgs& args);
}
#endif
