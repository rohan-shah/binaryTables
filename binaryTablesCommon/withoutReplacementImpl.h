#ifndef WITHOUT_REPLACEMENT_BINARY_TABLES_HEADER_GUARD
#define WITHOUT_REPLACEMENT_BINARY_TABLES_HEADER_GUARD
#include "sampford.h"
#include "problem.h"
#include "includeMPFRBinaryTables.h"
namespace binaryTables
{
	struct withoutReplacementArgs
	{
	public:
		withoutReplacementArgs(problem& problemObj)
			: problemObj(problemObj)
		{}
		std::size_t n;
		problem& problemObj;
		std::vector<int> sampleRowSums, newSampleRowSums, columnSums, newColumnSums;
		std::vector<mpfr_class> sizeVariables, newSizeVariables;
		std::vector<mpfr_class> productInclusionProbabilities, newProductInclusionProbabilities;
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		mpfr_class estimate;
		boost::mt19937 randomSource;
	};
	void withoutReplacement(withoutReplacementArgs& args);
}
#endif
