#ifndef WITHOUT_REPLACEMENT_MERGING_WITH_VARIANCE_BINARY_TABLES_HEADER_GUARD
#define WITHOUT_REPLACEMENT_MERGING_WITH_VARIANCE_BINARY_TABLES_HEADER_GUARD
#include "sampford.h"
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/shared_ptr.hpp>
#include "conditionalPoissonSequential.h"
namespace binaryTables
{
	struct withoutReplacementMergingWithVarianceSample
	{
		withoutReplacementMergingWithVarianceSample()
		{}
		/*withoutReplacementMergingWithVarianceSample(int columnSum, mpfr_class importanceDensity, mpfr_class weight, int trueDensity, int totalRemaining, boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants, boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters)
			:columnSum(columnSum), trueDensity(trueDensity), importanceDensity(importanceDensity), weight(weight), totalRemaining(totalRemaining), expNormalisingConstants(expNormalisingConstants), expExponentialParameters(expExponentialParameters)
		{}*/
		withoutReplacementMergingWithVarianceSample(withoutReplacementMergingWithVarianceSample&& other);
		withoutReplacementMergingWithVarianceSample& operator=(withoutReplacementMergingWithVarianceSample&& other);
		int columnSum;
		double trueDensity;
		mpfr_class importanceDensity;
		mpfr_class weight;
		int totalRemaining;
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants;
		boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters;
		int skipped;
		int nRemainingZeros;
		int nRemainingDeterministic;
		int parentIndexWithinDesign;
		int indexWithinDesign;
		std::vector<bool> deterministicInclusion;
	};
	struct withoutReplacementMergingWithVarianceArgs
	{
	public:
		withoutReplacementMergingWithVarianceArgs(problem& problemObj)
			: problemObj(problemObj), samplingArgs(true)
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		sampling::conditionalPoissonSequentialArgs samplingArgs;
		std::vector<withoutReplacementMergingWithVarianceSample> samples, newSamples;
		mpfr_class estimate, varianceEstimate;
		boost::mt19937 randomSource;
	};
	void withoutReplacementMergingWithVariance(withoutReplacementMergingWithVarianceArgs& args);
}
#endif
