#ifndef CONDITIONAL_POISSON_BOOTSTRAP_BINARY_TABLES_HEADER_GUARD
#define CONDITIONAL_POISSON_BOOTSTRAP_BINARY_TABLES_HEADER_GUARD
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/mersenne_twister.hpp>
namespace binaryTables
{
	struct conditionalPoissonBootstrapSample
	{
		conditionalPoissonBootstrapSample(int columnSum, mpfr_class density, boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants, boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters)
			:columnSum(columnSum), density(density), expNormalisingConstants(expNormalisingConstants), expExponentialParameters(expExponentialParameters)
		{}
		conditionalPoissonBootstrapSample(conditionalPoissonBootstrapSample&& other);
		conditionalPoissonBootstrapSample& operator=(conditionalPoissonBootstrapSample&& other);
		conditionalPoissonBootstrapSample& operator=(const conditionalPoissonBootstrapSample& other);
		conditionalPoissonBootstrapSample(const conditionalPoissonBootstrapSample& other);
		int columnSum;
		mpfr_class density;
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants;
		boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters;
		int skipped;
		int nRemainingZeros;
		int nRemainingDeterministic;
		std::vector<bool> deterministicInclusion;
	};
	struct conditionalPoissonBootstrapArgs
	{
	public:
		conditionalPoissonBootstrapArgs(problem& problemObj)
			: problemObj(problemObj) 
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		std::vector<conditionalPoissonBootstrapSample> samples, newSamples;
		mpfr_class estimate;
		boost::mt19937 randomSource;
	};
	void conditionalPoissonBootstrap(conditionalPoissonBootstrapArgs& args);
}
#endif
