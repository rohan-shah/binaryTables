#include "exhaustiveSearchImpl.h"
#include "GayleRyserTest.h"
#include <stdexcept>
namespace binaryTables
{
	std::size_t exhaustiveSearch(problem& problemObj)
	{
		const std::vector<int>& rowSums = problemObj.getRowSums();
		const std::vector<int>& columnSums = problemObj.getColumnSums();
		std::size_t nRows = rowSums.size(), nColumns = columnSums.size();
		std::size_t totalPossibilities = (1ULL << nRows * nColumns);
		if(totalPossibilities > (1ULL << 31))
		{
			throw std::runtime_error("Exhaustive enumeration can only be used if there are less than 2^32 possibilities");
		}
		std::size_t valid = 0;
		for(std::size_t i = 0; i < totalPossibilities; i++)
		{
			for(std::size_t row = 0; row < nRows; row++)
			{
				int sum = 0;
				for(std::size_t column = 0; column < nColumns; column++)
				{
					sum += (i & (1ULL << (row * nColumns + column))) > 0;
				}
				if(sum != rowSums[row]) goto notValid;
			}
			for(std::size_t column = 0; column < nColumns; column++)
			{
				int sum = 0;
				for(std::size_t row = 0; row < nRows; row++)
				{
					sum += (i & (1ULL << (row * nColumns + column))) > 0;
				}
				if(sum != columnSums[column]) goto notValid;
			}

			valid++;
		notValid:
			;
		}
		return valid;
	}
}
