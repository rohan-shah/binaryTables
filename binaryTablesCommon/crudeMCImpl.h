#ifndef CRUDE_MC_IMPL_HEADER_GUARD
#define CRUDE_MC_IMPL_HEADER_GUARD
#include "problem.h"
#include "includeMPFRBinaryTables.h"
#include <boost/random/mersenne_twister.hpp>
namespace binaryTables
{
	struct crudeMCArgs
	{
	public:
		crudeMCArgs(problem& problemObj)
			:problemObj(problemObj)
		{}
		problem& problemObj;
		std::size_t n;
		std::size_t valid;
		mpfr_class estimate;
		boost::mt19937 randomSource;
		bool GayleRyserTest;
	};
	void crudeMC(crudeMCArgs& args);
}
#endif
