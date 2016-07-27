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
		withoutReplacementSample(int columnSum, mpfr_class sizeVariable, mpfr_class productInclusionProbabilities, boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants, boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters)
			:columnSum(columnSum), sizeVariable(sizeVariable), productInclusionProbabilities(productInclusionProbabilities), expNormalisingConstants(expNormalisingConstants), expExponentialParameters(expExponentialParameters)
		{}
		withoutReplacementSample(withoutReplacementSample&& other);
		withoutReplacementSample& operator=(withoutReplacementSample&& other);
		int columnSum;
		mpfr_class sizeVariable;
		mpfr_class productInclusionProbabilities;
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants;
		boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters;
		int skipped;
		int nRemainingZeros;
		int nRemainingDeterministic;
		std::vector<bool> table, deterministicInclusion;
	};
	struct withoutReplacementArgs
	{
	public:
		withoutReplacementArgs(problem& problemObj)
			: problemObj(problemObj), keepTables(false)
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		std::vector<withoutReplacementSample> samples, newSamples;
		mpfr_class estimate;
		boost::mt19937 randomSource;
		bool keepTables;
		std::vector<bool> tables;
	};
	void withoutReplacement(withoutReplacementArgs& args);
}
#endif
