#ifndef WITHOUT_REPLACEMENT_BINARY_TABLES_HEADER_GUARD
#define WITHOUT_REPLACEMENT_BINARY_TABLES_HEADER_GUARD
#include "sampford.h"
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/numeric/ublas/matrix.hpp>
namespace binaryTables
{
	struct withoutReplacementSample
	{
		withoutReplacementSample(int columnSum, mpfr_class sizeVariable, mpfr_class productInclusionProbabilities, int totalRemaining, boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants, boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters)
			:columnSum(columnSum), sizeVariable(sizeVariable), productInclusionProbabilities(productInclusionProbabilities), totalRemaining(totalRemaining), expNormalisingConstants(expNormalisingConstants), expExponentialParameters(expExponentialParameters)
		{}
		withoutReplacementSample(withoutReplacementSample&& other);
		withoutReplacementSample& operator=(withoutReplacementSample&& other);
		int columnSum;
		mpfr_class sizeVariable;
		mpfr_class productInclusionProbabilities;
		int totalRemaining;
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants;
		boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters;
		int skipped;
		int nRemainingZeros;
		int nRemainingDeterministic;
	};
	struct withoutReplacementArgs
	{
	public:
		withoutReplacementArgs(problem& problemObj)
			: problemObj(problemObj)
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		std::vector<withoutReplacementSample> samples, newSamples;
		mpfr_class estimate;
		boost::mt19937 randomSource;
	};
	void withoutReplacement(withoutReplacementArgs& args);
}
#endif
