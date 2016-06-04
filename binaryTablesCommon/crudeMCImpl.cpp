#include "crudeMCImpl.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "GayleRyserTest.h"
namespace binaryTables
{
	void crudeMC(crudeMCArgs& args)
	{
		problem& problemObj = args.problemObj;
		const std::vector<int>& rowSums = problemObj.getRowSums();
		const std::vector<int>& columnSums = problemObj.getColumnSums();
		if(std::find(rowSums.begin(), rowSums.end(), 0) != rowSums.end())
		{
			throw std::runtime_error("Input rowSums cannot contain 0");
		}
		if(std::find(columnSums.begin(), columnSums.end(), 0) != columnSums.end())
		{
			throw std::runtime_error("Input columnSums cannot contain 0");
		}

		std::size_t nRows = rowSums.size(), nColumns = columnSums.size();
		std::size_t n = args.n;
		args.valid = 0;
		std::vector<int> currentRowSums(nRows), currentColumnSums(nColumns);
		boost::random::bernoulli_distribution<> bernoulli(0.5);
		GayleRyserTestWorking working(true);
		for(std::size_t i = 0; i < n; i++)
		{
			std::copy(rowSums.begin(), rowSums.end(), currentRowSums.begin());
			std::copy(columnSums.begin(), columnSums.end(), currentColumnSums.begin());
			bool valid = true;
			for(std::size_t row = 0; row < nRows; row++)
			{
				for(std::size_t column = 0; column < nColumns; column++)
				{
					if(bernoulli(args.randomSource))
					{
						if(currentRowSums[row] == 0 || currentColumnSums[column] == 0) 
						{
							valid = false;
							goto endCurrent;
						}
						currentRowSums[row]--;
						currentColumnSums[column]--;
					}
				}
				if(currentRowSums[row] != 0)
				{
					valid = false;
					goto endCurrent;
				}
				if(args.GayleRyserTest && !GayleRyserTest(currentRowSums, currentColumnSums, row+1, working))
				{
					valid = false;
					goto endCurrent;
				}
			}
			for(std::size_t row = 0; row < nRows; row++) 
			{
				if(currentRowSums[row] != 0)
				{
					valid = false; goto endCurrent;
				}
			}
			for(std::size_t column = 0; column < nColumns; column++)
			{
				if(currentColumnSums[column] != 0)
				{
					valid = false; goto endCurrent;
				}
			}
endCurrent:
			if(valid) args.valid++;
		}
		mpfr_class probability = mpfr_class(args.valid)/mpfr_class(n);
		args.estimate = boost::multiprecision::pow(mpfr_class(2), rowSums.size() * columnSums.size()) * probability;
	}
}
