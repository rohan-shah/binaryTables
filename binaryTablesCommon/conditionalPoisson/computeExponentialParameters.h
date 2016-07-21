#include <boost/numeric/ublas/matrix.hpp>
#include "includeMPFRBinaryTables.h"
#include "conditionalPoisson/conditionalPoissonSequential.h"
namespace binaryTables
{
	namespace conditionalPoisson
	{
		void exponentialPreCompute(computeExponentialData& data, int nColumns);
		void computeExponentialParameters(conditionalPoissonSequentialArgs& args, int nColumns, std::vector<int>::const_iterator rowSumsStart, std::vector<int>::const_iterator rowSumsEnd);
	}
}
