#include "GayleRyserTest.h"
#include <algorithm>
namespace binaryTables
{
	bool GayleRyserTest(std::vector<int>& sums1, std::vector<int>& sums2, int sums1ToIgnore, GayleRyserTestWorking& working)
	{
		if(working.canReorder)
		{
			std::sort(sums1.begin() + sums1ToIgnore, sums1.end());
			std::sort(sums2.begin(), sums2.end());
		}
		else
		{
			working.sortedSums1 = sums1;
			working.sortedSums2 = sums2;
			std::sort(working.sortedSums1.begin() + sums1ToIgnore, working.sortedSums1.end());
			std::sort(working.sortedSums2.begin(), working.sortedSums2.end());
		}
		std::vector<int>& sortedSums1 = (working.canReorder ? sums1 : working.sortedSums1);
		std::vector<int>& sortedSums2 = (working.canReorder ? sums2 : working.sortedSums2);
		working.table.resize(sortedSums2.size()+1);
		std::fill(working.table.begin(), working.table.end(), 0);
		for(std::size_t i = sums1ToIgnore; i < sortedSums1.size(); i++)
		{
			working.table[sortedSums1[i]]++;
		}
		working.starVector.resize(sortedSums2.size());
		int accumulated = sortedSums1.size()-sums1ToIgnore;
		for(std::size_t i = 0; i < sortedSums2.size(); i++)
		{
			working.starVector[i] = accumulated;
			accumulated -= working.table[i+1];
		}
		std::sort(working.starVector.begin(), working.starVector.end());
		bool valid = true;
		int accumulatedStar = 0;
		int accumulatedSums = 0;
		for(std::size_t i = 0; i < working.starVector.size(); i++)
		{
			accumulatedStar += *(working.starVector.rbegin() + i);
			accumulatedSums += *(sortedSums2.rbegin() + i);
			if(accumulatedSums > accumulatedStar)
			{
				valid = false;
				break;
			}
		}
		if(accumulatedStar != accumulatedSums) valid = false;
		return valid;
	}
}
