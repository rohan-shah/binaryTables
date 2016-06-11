#include "withoutReplacementImpl.h"
namespace binaryTables
{
	void withoutReplacement(withoutReplacementArgs& args)
	{
		problem& problemObj = args.problemObj;
		const std::vector<int>& initialRowSums = problemObj.getRowSums();
		const std::vector<int>& initialColumnSums = problemObj.getColumnSums();
		std::size_t nRows = initialRowSums.size(), nColumns = initialColumnSums.size();

		int totalOnes = 0;
		for(std::vector<int>::const_iterator i = initialRowSums.begin(); i != initialRowSums.end(); i++) totalOnes += *i;
		
		std::size_t n = args.n;

		std::vector<int>& sampleRowSums = args.sampleRowSums;
		std::vector<int>& newSampleRowSums = args.newSampleRowSums;
		std::vector<int>& columnSums = args.columnSums;
		std::vector<int>& newColumnSums = args.newColumnSums;
		std::vector<int>& totalRemaining = args.totalRemaining;
		std::vector<int>& newTotalRemaining = args.newTotalRemaining;
		std::vector<mpfr_class>& sizeVariables = args.sizeVariables;
		std::vector<mpfr_class>& newSizeVariables = args.newSizeVariables;
		std::vector<mpfr_class>& productInclusionProbabilities = args.productInclusionProbabilities;
		std::vector<mpfr_class>& newProductInclusionProbabilities = args.newProductInclusionProbabilities;

		std::vector<mpfr_class> sampfordInclusionProbabilities, sampfordRescaledWeights;
		std::vector<int> sampfordIndices;
		args.samplingArgs.inclusionProbabilities = &sampfordInclusionProbabilities;
		args.samplingArgs.n = n;
		args.samplingArgs.weights = &newSizeVariables;
		args.samplingArgs.rescaledWeights = &sampfordRescaledWeights;
		args.samplingArgs.indices = &sampfordIndices;

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
		totalRemaining.clear();
		newTotalRemaining.clear();
		//In this case we can't have a 1 in the first entry
		if(initialRowSums[0] == 0)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			sizeVariables.push_back(1);
			productInclusionProbabilities.push_back(1);
			columnSums[0] = 0;
			totalRemaining.push_back(totalOnes);
		}
		else if(initialRowSums[0] == (int)nColumns)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			sizeVariables.push_back(1);
			productInclusionProbabilities.push_back(1);
			sampleRowSums[0]--;
			columnSums[0] = 1;
			totalRemaining.push_back(totalOnes-1);
		}
		else
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin() + nRows);
			sampleRowSums[0]--;
			//The first sample has a one in the first column
			columnSums[0] = 1;
			columnSums[1] = 0;
			sizeVariables.push_back(initialRowSums[0] / (double)totalOnes);
			sizeVariables.push_back(1 - initialRowSums[0] / (double)totalOnes);
			productInclusionProbabilities.push_back(1);
			productInclusionProbabilities.push_back(1);

			totalRemaining.push_back(totalOnes-1);
			totalRemaining.push_back(totalOnes);
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
				newTotalRemaining.clear();
				if(choicesUp.size() + choicesDown.size() <= n)
				{
					std::size_t outputCounter = 0;
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						std::copy(sampleRowSums.begin() + nRows * choicesUp[i], sampleRowSums.begin() + nRows * (choicesUp[i]+1), newSampleRowSums.begin() + outputCounter * nRows);
						int previousRowSumValue = newSampleRowSums[nRows * outputCounter + row];
						newSampleRowSums[nRows * outputCounter + row]--;
						newSizeVariables.push_back(sizeVariables[choicesUp[i]]*previousRowSumValue / (double)totalRemaining[choicesUp[i]]);
						newProductInclusionProbabilities.push_back(productInclusionProbabilities[choicesUp[i]]);
						newColumnSums[outputCounter] = columnSums[choicesUp[i]]+1;
						newTotalRemaining.push_back(totalRemaining[choicesUp[i]]-1);
						outputCounter++;
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						std::copy(sampleRowSums.begin() + nRows * choicesDown[i], sampleRowSums.begin() + nRows * (choicesDown[i]+1), newSampleRowSums.begin() + outputCounter * nRows);
						int previousRowSumValue = newSampleRowSums[nRows * outputCounter + row];
						newSizeVariables.push_back(sizeVariables[choicesDown[i]]*(1 - previousRowSumValue / (double)totalRemaining[choicesDown[i]]));
						newProductInclusionProbabilities.push_back(productInclusionProbabilities[choicesDown[i]]);
						newColumnSums[outputCounter] = columnSums[choicesDown[i]];
						newTotalRemaining.push_back(totalRemaining[choicesDown[i]]);
						outputCounter++;
					}
					sizeVariables.swap(newSizeVariables);
				}
				else
				{
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						newSizeVariables.push_back(sizeVariables[choicesUp[i]]*sampleRowSums[nRows * choicesUp[i] + row] / (double)totalRemaining[choicesUp[i]]);
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						newSizeVariables.push_back(sizeVariables[choicesDown[i]]*(1 - sampleRowSums[nRows * choicesDown[i] + row] / (double)totalRemaining[choicesDown[i]]));
					}
					sampfordFromParetoNaive(args.samplingArgs, args.randomSource);
					sizeVariables.clear();
					newTotalRemaining.clear();
					for(std::size_t i = 0; i < n; i++)
					{
						if(sampfordIndices[i] >= (int)choicesUp.size())
						{
							int parent = choicesDown[sampfordIndices[i] - choicesUp.size()];
							std::copy(sampleRowSums.begin() + nRows * parent, sampleRowSums.begin() + nRows * (parent+1), newSampleRowSums.begin() + i * nRows);
							newProductInclusionProbabilities.push_back(productInclusionProbabilities[parent]*sampfordInclusionProbabilities[sampfordIndices[i]]);
							newColumnSums[i] = columnSums[parent];
							newTotalRemaining.push_back(totalRemaining[parent]);
							sizeVariables.push_back(newSizeVariables[sampfordIndices[i]]/sampfordInclusionProbabilities[sampfordIndices[i]]);
						}
						else
						{
							int parent = choicesUp[sampfordIndices[i]];
							std::copy(sampleRowSums.begin() + nRows * parent, sampleRowSums.begin() + nRows * (parent+1), newSampleRowSums.begin() + i * nRows);
							newSampleRowSums[i*nRows + row]--;
							newProductInclusionProbabilities.push_back(productInclusionProbabilities[parent]*sampfordInclusionProbabilities[sampfordIndices[i]]);
							newColumnSums[i] = columnSums[parent]+1;
							newTotalRemaining.push_back(totalRemaining[parent]-1);
							sizeVariables.push_back(newSizeVariables[sampfordIndices[i]]/sampfordInclusionProbabilities[sampfordIndices[i]]);
						}
					}
				}
				totalRemaining.swap(newTotalRemaining);
				productInclusionProbabilities.swap(newProductInclusionProbabilities);
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
