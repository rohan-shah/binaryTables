#ifndef GAYLE_RYSER_TEST_HEADER_GUARD
#define GAYLE_RYSER_TEST_HEADER_GUARD
#include <vector>
namespace binaryTables
{
	struct GayleRyserTestWorking
	{
		GayleRyserTestWorking(bool canReorder)
			:canReorder(canReorder)
		{}
		std::vector<int> table, sortedSums1, sortedSums2, starVector;
		bool canReorder;
	};
	bool GayleRyserTest(std::vector<int>& sums1, std::vector<int>& sums2, int sums1ToIgnore, GayleRyserTestWorking& working);
}
#endif
