#include "withoutReplacementImpl.h"
namespace binaryTables
{
	void withoutReplacement(withoutReplacementArgs& args)
	{
		problem& problemObj = args.problemObj;
		const std::vector<int>& initialRowSums = problemObj.getRowSums();
		const std::vector<int>& initialColumnSums = problemObj.getColumnSums();
		std::size_t nRows = initialRowSums.size(), nColumns = initialColumnSums.size();
		std::size_t n = args.n;

		std::vector<int>& sampleRowSums = args.sampleRowSums;
		std::vector<int>& newSampleRowSums = args.newSampleRowSums;
		std::vector<int>& columnSums = args.columnSums;
		std::vector<int>& newColumnSums = args.newColumnSums;
		std::vector<mpfr_class>& sizeVariables = args.sizeVariables;
		std::vector<mpfr_class>& newSizeVariables = args.newSizeVariables;
		std::vector<mpfr_class>& productInclusionProbabilities = args.productInclusionProbabilities;
		std::vector<mpfr_class>& newProductInclusionProbabilities = args.newProductInclusionProbabilities;

		if(n < 2)
		{
			throw std::runtime_error("Input n must be at least 2");
		}
		if(nRows <= 1 || nColumns <= 1)
		{
			throw std::runtime_error("The number of columns and rows must be at least 2");
		}
		columnSums.resize(n);
		newColumnSums.resize(n);
		sampleRowSums.resize(nRows*n);
		newSampleRowSums.resize(nRows*n);
		sizeVariables.clear();
		newSizeVariables.clear();
		productInclusionProbabilities.clear();
		newProductInclusionProbabilities.clear();
		//In this case we can't have a 1 in the first entry
		if(initialRowSums[0] == 0)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			sizeVariables.push_back(1);
			productInclusionProbabilities.push_back(1);
			columnSums[0] = 0;
		}
		else if(initialRowSums[0] == (int)nColumns)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			sizeVariables.push_back(1);
			productInclusionProbabilities.push_back(1);
			sampleRowSums[0]--;
			columnSums[0] = 1;
		}
		else
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin() + nRows);
			sampleRowSums[0]--;
			//The first sample has a one in the first column
			columnSums[0] = 1;
			columnSums[1] = 0;
			int sumRowSums = 0;
			for(int i = 0; i < (int)nRows; i++) sumRowSums += initialRowSums[i];
			sizeVariables.push_back(initialRowSums[0] / (double)sumRowSums);
			sizeVariables.push_back(1 - initialRowSums[0] / (double)sumRowSums);
			productInclusionProbabilities.push_back(1);
			productInclusionProbabilities.push_back(1);
		}
		std::vector<int> choicesUp, choicesDown;
		int row = 1, column = 0;
		for(; column < (int)nColumns; column++)
		{
			for(; row < (int)nRows; row++)
			{
				choicesUp.clear();
				choicesDown.clear();
				for(std::size_t i = 0; i < productInclusionProbabilities.size(); i++)
				{
					if(sampleRowSums[i*nRows + row] == 0)
					{
						if(columnSums[i] + (int)nRows - row- 1 >= initialColumnSums[column] && columnSums[i] <= initialColumnSums[column]) choicesDown.push_back((int)i);
					}
					else if(sampleRowSums[i*nRows + row] == (int)nColumns - column)
					{
						if(columnSums[i]+1 <= initialColumnSums[column] && columnSums[i] + (int)nRows - row>= initialColumnSums[column]) choicesUp.push_back((int)i);
					}
					else
					{
						if(columnSums[i]+1 <= initialColumnSums[column] && columnSums[i] + (int)nRows - row >= initialColumnSums[column]) choicesUp.push_back((int)i);
						if(columnSums[i] + (int)nRows - row - 1 >= initialColumnSums[column] && columnSums[i] <= initialColumnSums[column]) choicesDown.push_back((int)i);
					}
				}
				newProductInclusionProbabilities.clear();
				newSizeVariables.clear();
				if(choicesUp.size() + choicesDown.size() <= n)
				{
					std::size_t outputCounter = 0;
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						std::copy(sampleRowSums.begin() + nRows * choicesUp[i], sampleRowSums.begin() + nRows * (choicesUp[i]+1), newSampleRowSums.begin() + outputCounter * nRows);
						int previousRowSumValue = newSampleRowSums[nRows * outputCounter + row];
						newSampleRowSums[nRows * outputCounter + row]--;
						int sumRowSums = 0;
						for(int j = 0; j < (int)nRows; j++) sumRowSums += sampleRowSums[nRows * choicesUp[i] + j];
						newSizeVariables.push_back(sizeVariables[choicesUp[i]]*previousRowSumValue / (double)sumRowSums);
						newProductInclusionProbabilities.push_back(productInclusionProbabilities[choicesUp[i]]);
						newColumnSums[outputCounter] = columnSums[choicesUp[i]]+1;
						outputCounter++;
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						std::copy(sampleRowSums.begin() + nRows * choicesDown[i], sampleRowSums.begin() + nRows * (choicesDown[i]+1), newSampleRowSums.begin() + outputCounter * nRows);
						int previousRowSumValue = newSampleRowSums[nRows * outputCounter + row];
						int sumRowSums = 0;
						for(int j = 0; j < (int)nRows; j++) sumRowSums += sampleRowSums[nRows * choicesDown[i] + j];
						newSizeVariables.push_back(sizeVariables[choicesDown[i]]*(1 - previousRowSumValue / (double)sumRowSums));
						newProductInclusionProbabilities.push_back(productInclusionProbabilities[choicesDown[i]]);
						newColumnSums[outputCounter] = columnSums[choicesDown[i]];
						outputCounter++;
					}
				}
				else
				{
					throw std::runtime_error("Not implemented yet");
				}
				productInclusionProbabilities.swap(newProductInclusionProbabilities);
				sizeVariables.swap(newSizeVariables);
				sampleRowSums.swap(newSampleRowSums);
				columnSums.swap(newColumnSums);
			}
			std::fill(columnSums.begin(), columnSums.end(), 0);
			row = 0;
		}
		args.estimate = 0;
		for(std::size_t i = 0; i < productInclusionProbabilities.size(); i++)
		{
			args.estimate += 1/productInclusionProbabilities[i];
		}
	}
}
