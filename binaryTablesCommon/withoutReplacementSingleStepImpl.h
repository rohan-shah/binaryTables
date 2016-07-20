#ifndef WITHOUT_REPLACEMENT_SINGLE_STEP_BINARY_TABLES_HEADER_GUARD
#define WITHOUT_REPLACEMENT_SINGLE_STEP_BINARY_TABLES_HEADER_GUARD
#include "sampford.h"
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "withoutReplacementImpl.h"
namespace binaryTables
{
	struct withoutReplacementSingleStepSample
	{
		withoutReplacementSingleStepSample(int columnSum, mpfr_class productInclusionProbabilities, int totalRemaining, boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants, boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters)
			:columnSum(columnSum), productInclusionProbabilities(productInclusionProbabilities), totalRemaining(totalRemaining), expNormalisingConstants(expNormalisingConstants), expExponentialParameters(expExponentialParameters)
		{}
		withoutReplacementSingleStepSample(withoutReplacementSingleStepSample&& other);
		withoutReplacementSingleStepSample& operator=(withoutReplacementSingleStepSample&& other);
		int columnSum;
		mpfr_class productInclusionProbabilities;
		int totalRemaining;
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants;
		boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters;
		int skipped;
		int nRemainingZeros;
		int nRemainingDeterministic;
		std::vector<bool> deterministicInclusion;
	};
	struct withoutReplacementSingleStepArgs
	{
	public:
		withoutReplacementSingleStepArgs(problem& problemObj)
			: problemObj(problemObj)
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		std::vector<withoutReplacementSingleStepSample> samples, newSamples;
		mpfr_class estimate;
		boost::mt19937 randomSource;
	};
	void withoutReplacementSingleStep(withoutReplacementSingleStepArgs& args);
}
#endif
