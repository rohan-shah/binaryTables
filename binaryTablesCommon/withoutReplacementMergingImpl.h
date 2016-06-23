#ifndef WITHOUT_REPLACEMENT_MERGING_BINARY_TABLES_HEADER_GUARD
#define WITHOUT_REPLACEMENT_MERGING_BINARY_TABLES_HEADER_GUARD
#include "sampford.h"
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/numeric/ublas/matrix.hpp>
namespace binaryTables
{
	struct withoutReplacementMergingSample
	{
		withoutReplacementMergingSample(int columnSum, mpfr_class sizeVariable, mpfr_class weight, int totalRemaining, boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants, boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters)
			:columnSum(columnSum), sizeVariable(sizeVariable), weight(weight), totalRemaining(totalRemaining), expNormalisingConstants(expNormalisingConstants), expExponentialParameters(expExponentialParameters)
		{}
		withoutReplacementMergingSample(withoutReplacementMergingSample&& other);
		withoutReplacementMergingSample& operator=(withoutReplacementMergingSample&& other);
		int columnSum;
		mpfr_class sizeVariable;
		mpfr_class weight;
		int totalRemaining;
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > expNormalisingConstants;
		boost::shared_ptr<std::vector<mpfr_class> > expExponentialParameters;
		int skipped;
		int nRemainingZeros;
		int nRemainingDeterministic;
	};
	struct withoutReplacementMergingArgs
	{
	public:
		withoutReplacementMergingArgs(problem& problemObj)
			: problemObj(problemObj)
		{}
		std::size_t n;
		int mergeFrequency;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums;
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		std::vector<withoutReplacementMergingSample> samples, newSamples;
		mpfr_class estimate;
		boost::mt19937 randomSource;
	};
	void withoutReplacementMerging(withoutReplacementMergingArgs& args);
}
#endif
