#ifndef CONDITIONAL_POISSON_IMPL_HEADER_GUARD
#define CONDITIONAL_POISSON_IMPL_HEADER_GUARD
#include "problem.h"
#include "conditionalPoissonDrafting.h"
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
		mpfr_class estimate;
		sampling::conditionalPoissonDraftingArgs samplingArgs;
		std::vector<mpfr_class> inclusionProbablities, samplingWeights;
		boost::mt19937 randomSource;
	private:

	};
	void conditionalPoisson(conditionalPoissonArgs& args);
}
#endif
