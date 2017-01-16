#ifndef CONDITIONAL_POISSON_BOOTSTRAP_MERGING_BINARY_TABLES_HEADER_GUARD
#define CONDITIONAL_POISSON_BOOTSTRAP_MERGING_BINARY_TABLES_HEADER_GUARD
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "conditionalPoissonBootstrapImpl.h"
namespace binaryTables
{
	struct conditionalPoissonBootstrapMergingArgs
	{
	public:
		conditionalPoissonBootstrapMergingArgs(problem& problemObj)
			: problemObj(problemObj) 
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		std::vector<conditionalPoissonBootstrapSample> samples, newSamples;
		mpfr_class estimate;
		boost::mt19937 randomSource;
		int mergeFrequency;
	};
	void conditionalPoissonBootstrapMerging(conditionalPoissonBootstrapMergingArgs& args);
}
#endif
