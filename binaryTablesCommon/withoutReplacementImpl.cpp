#include "withoutReplacementImpl.h"
#include "conditionalPoissonSequential.h"
#include "samplingBase.h"
namespace binaryTables
{
	void withoutReplacement(withoutReplacementArgs& args)
	{
		problem& problemObj = args.problemObj;
		const std::vector<int>& initialRowSums = problemObj.getRowSums();
		const std::vector<int>& initialColumnSums = problemObj.getColumnSums();
		std::size_t nRows = initialRowSums.size(), nColumns = initialColumnSums.size();

		//Count the total number of ones
		int totalOnes = 0;
		for(std::vector<int>::const_iterator i = initialRowSums.begin(); i != initialRowSums.end(); i++) totalOnes += *i;
		
		std::size_t n = args.n;

		//Extract data from args
		std::vector<int>& sampleRowSums = args.sampleRowSums;
		std::vector<int>& newSampleRowSums = args.newSampleRowSums;
		std::vector<withoutReplacementSample>& samples = args.samples;
		std::vector<withoutReplacementSample>& newSamples = args.newSamples;
		std::vector<mpfr_class> conditionalPoissonInclusionProbabilities;
		sampling::conditionalPoissonSequentialArgs cpSamplingArgs(true);
		cpSamplingArgs.n = initialColumnSums[0];

/*		if(n < 2)
		{
			throw std::runtime_error("Input n must be at least 2");
		}*/
		if(nRows <= 1 || nColumns <= 1)
		{
			throw std::runtime_error("The number of columns and rows must be at least 2");
		}
		//Set up data structures for the particles
		sampleRowSums.resize(nRows*std::max((std::size_t)2, n));
		newSampleRowSums.resize(nRows*std::max((std::size_t)2, n));
		samples.clear();
		newSamples.clear();
		//We are going to perform conditional Poisson sampling using the sequential method. Se we need the selection probabilities. 
		//This requires that we compute the exponentials of the normalising constants. 
		cpSamplingArgs.weights.clear();
		for(std::size_t i = 0; i < nRows; i++) cpSamplingArgs.weights.push_back(initialRowSums[i]);
		sampling::samplingBase(initialColumnSums[0], cpSamplingArgs.indices, cpSamplingArgs.weights, cpSamplingArgs.rescaledWeights, cpSamplingArgs.zeroWeights, cpSamplingArgs.deterministicInclusion, conditionalPoissonInclusionProbabilities);

		computeExponentialParameters(cpSamplingArgs);
		sampling::conditionalPoissonInclusionProbabilities(cpSamplingArgs, cpSamplingArgs.inclusionProbabilities);
		mpfr_class selectionProb;
		if(initialRowSums[0] == 0)
		{
			selectionProb = 0;
		}
		else if(initialColumnSums[0] == 1)
		{
			selectionProb = cpSamplingArgs.expExponentialParameters[0] / cpSamplingArgs.expNormalisingConstant(0, initialColumnSums[0] - 1);
		}
		else
		{
			selectionProb = cpSamplingArgs.expExponentialParameters[0] * cpSamplingArgs.expNormalisingConstant(1, initialColumnSums[0]- 2) / cpSamplingArgs.expNormalisingConstant(0, initialColumnSums[0] - 1);
		}
		//Swap out the expNormalisingConstant data into a shared ptr
		boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > initialExpNormalisingConstantData(new boost::numeric::ublas::matrix<mpfr_class>());
		(*initialExpNormalisingConstantData).swap(cpSamplingArgs.expNormalisingConstant);
		//Similarly for expExponentialParameters
		boost::shared_ptr<std::vector<mpfr_class> > initialExpExponentialParameters(new std::vector<mpfr_class>());
		(*initialExpExponentialParameters).swap(cpSamplingArgs.expExponentialParameters);
		//In this case we can't have a 1 in the first entry
		if(initialRowSums[0] == 0)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			samples.push_back(withoutReplacementSample(0, 1 - selectionProb, 1, totalOnes));
			samples[0].skipped = 1;
		}
		else if(initialRowSums[0] == (int)nColumns)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			samples.push_back(withoutReplacementSample(1, selectionProb, 1, totalOnes-1));
			sampleRowSums[0]--;
			samples[0].skipped = 1;
		}
		else
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin() + nRows);
			sampleRowSums[0]--;
			//The first sample has a one in the first column. To start off with, don't 
			samples.push_back(withoutReplacementSample(1, selectionProb, 1, totalOnes-1));
			samples.push_back(withoutReplacementSample(0, 1 - selectionProb, 1, totalOnes));
			samples[0].skipped = samples[1].skipped = 0;
		}
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			samples[i].expExponentialParameters = initialExpExponentialParameters;
			samples[i].expNormalisingConstants = initialExpNormalisingConstantData;
			//This is to do with skipping over the rows which have a row-sum of zero
			samples[i].nRemainingZeros = samples[i].nRemainingDeterministic = 0;
			for(int j = 0; j < (int)nRows; j++)
			{
				if(sampleRowSums[i * nRows + j] == 0) samples[i].nRemainingZeros++;
				if(sampleRowSums[i * nRows + j] == (int)nColumns) samples[i].nRemainingDeterministic++;
			}
		}

		std::vector<int> choicesUp, choicesDown;
		int row = 1, column = 0;
		for(; column < (int)nColumns; column++)
		{
			cpSamplingArgs.n = initialColumnSums[column];
			for(; row < (int)nRows; row++)
			{
				choicesUp.clear();
				choicesDown.clear();
				for(std::size_t i = 0; i < samples.size(); i++)
				{
					if(sampleRowSums[i*nRows + row] == 0)
					{
						if(samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column] && samples[i].columnSum <= initialColumnSums[column]) choicesDown.push_back((int)i);
					}
					else if(sampleRowSums[i*nRows + row] == (int)nColumns - column)
					{
						if(samples[i].columnSum + 1 <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros>= initialColumnSums[column]) choicesUp.push_back((int)i);
					}
					else
					{
						if(samples[i].columnSum + 1 <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros>= initialColumnSums[column]) choicesUp.push_back((int)i);
						if(samples[i].columnSum + (int)nRows - row - 1 - samples[i].nRemainingZeros >= initialColumnSums[column] && samples[i].columnSum <= initialColumnSums[column]) choicesDown.push_back((int)i);
					}
				}
				newSamples.clear();
				if(choicesUp.size() + choicesDown.size() <= n)
				{
					std::size_t outputCounter = 0;
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						int parentIndex = choicesUp[i];
						withoutReplacementSample& parentSample = samples[choicesUp[i]];
						std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + outputCounter * nRows);
						int newSkipped = parentSample.skipped, newDeterministic = parentSample.nRemainingDeterministic;
						if(sampleRowSums[nRows * parentIndex + row] == (int)nColumns - column)
						{
							selectionProb = 1;
							newSkipped++;
							newDeterministic--;
						}
						else if(parentSample.columnSum == initialColumnSums[column] - 1 - parentSample.nRemainingDeterministic)
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] / (*parentSample.expNormalisingConstants)(row-parentSample.skipped, 0);
						}
						else
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] * (*parentSample.expNormalisingConstants)(row+1 - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 2 - parentSample.nRemainingDeterministic) / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 1 - parentSample.nRemainingDeterministic);
						}
						newSampleRowSums[nRows * outputCounter + row]--;
						newSamples.push_back(withoutReplacementSample(parentSample.columnSum + 1, parentSample.sizeVariable * selectionProb, parentSample.productInclusionProbabilities, parentSample.totalRemaining - 1));
						newSamples.back().expExponentialParameters = parentSample.expExponentialParameters;
						newSamples.back().expNormalisingConstants = parentSample.expNormalisingConstants;
						newSamples.back().skipped = newSkipped;
						newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
						newSamples.back().nRemainingDeterministic = newDeterministic;
						outputCounter++;
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						int parentIndex = choicesDown[i];
						withoutReplacementSample& parentSample = samples[parentIndex];
						std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + outputCounter * nRows);
						int newSkipped = parentSample.skipped;
						if(parentSample.columnSum == initialColumnSums[column] - parentSample.nRemainingDeterministic || newSampleRowSums[outputCounter * nRows + row] == 0)
						{
							selectionProb = 0;
							newSkipped++;
						}
						else if(parentSample.columnSum == initialColumnSums[column] - 1 - parentSample.nRemainingDeterministic)
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, 0);
						}
						else
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] * (*parentSample.expNormalisingConstants)(row+1 - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 2 - parentSample.nRemainingDeterministic) / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 1 - parentSample.nRemainingDeterministic);
						}
						newSamples.push_back(withoutReplacementSample(parentSample.columnSum, parentSample.sizeVariable * (1 - selectionProb), parentSample.productInclusionProbabilities, parentSample.totalRemaining));
						newSamples.back().expExponentialParameters = parentSample.expExponentialParameters;
						newSamples.back().expNormalisingConstants = parentSample.expNormalisingConstants;
						newSamples.back().skipped = newSkipped;
						newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
						newSamples.back().nRemainingDeterministic = parentSample.nRemainingDeterministic;
						if(newSampleRowSums[outputCounter * nRows + row] == 0) newSamples.back().nRemainingZeros--;
						outputCounter++;
					}
				}
				else
				{
					args.samplingArgs.weights.clear();
					args.samplingArgs.n = n;
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						int parentIndex = choicesUp[i];
						withoutReplacementSample& parentSample = samples[parentIndex];
						if(sampleRowSums[nRows * parentIndex + row] == (int)nColumns - column)
						{
							selectionProb = 1;
						}
						else if(parentSample.columnSum == initialColumnSums[column]-1)
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, 0);
						}
						else
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] * (*parentSample.expNormalisingConstants)(row+1 - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 2 - parentSample.nRemainingDeterministic) / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 1 - parentSample.nRemainingDeterministic);
						}
#ifndef NDEBUG
						if(selectionProb.convert_to<double>() > 1) throw std::runtime_error("Internal error");
#endif
						args.samplingArgs.weights.push_back(parentSample.sizeVariable * selectionProb);
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						int parentIndex = choicesDown[i];
						withoutReplacementSample& parentSample = samples[parentIndex];
						if(parentSample.columnSum == initialColumnSums[column] - parentSample.nRemainingDeterministic|| sampleRowSums[parentIndex*nRows + row] == 0)
						{
							selectionProb = 0;
						}
						else if(parentSample.columnSum == initialColumnSums[column]-1 - parentSample.nRemainingDeterministic)
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, 0);
						}
						else
						{
							selectionProb = (*parentSample.expExponentialParameters)[row] * (*parentSample.expNormalisingConstants)(row+1 - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 2 - parentSample.nRemainingDeterministic) / (*parentSample.expNormalisingConstants)(row - parentSample.skipped, initialColumnSums[column] - parentSample.columnSum - 1 - parentSample.nRemainingDeterministic);
						}
#ifndef NDEBUG
						if(selectionProb.convert_to<double>() > 1) throw std::runtime_error("Internal error");
#endif
						args.samplingArgs.weights.push_back(parentSample.sizeVariable * (1 - selectionProb));
					}
					sampfordFromParetoNaive(args.samplingArgs, args.randomSource);
					for(std::size_t i = 0; i < n; i++)
					{
						int selected = args.samplingArgs.indices[i];
						mpfr_class& inclusionProbability = args.samplingArgs.rescaledWeights[selected];
						if(args.samplingArgs.indices[i] >= (int)choicesUp.size())
						{
							int parentIndex = choicesDown[selected - choicesUp.size()];
							withoutReplacementSample& parentSample = samples[parentIndex];
							std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + i * nRows);
							newSamples.push_back(withoutReplacementSample(parentSample.columnSum, args.samplingArgs.weights[selected] / inclusionProbability, parentSample.productInclusionProbabilities * inclusionProbability, parentSample.totalRemaining));
							newSamples.back().expExponentialParameters = parentSample.expExponentialParameters;
							newSamples.back().expNormalisingConstants = parentSample.expNormalisingConstants;
							newSamples.back().skipped = parentSample.skipped;
							newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
							newSamples.back().nRemainingDeterministic = parentSample.nRemainingDeterministic;
							if(sampleRowSums[parentIndex*nRows + row] == 0) newSamples.back().nRemainingZeros--;
							if(parentSample.columnSum == initialColumnSums[column] || sampleRowSums[parentIndex*nRows + row] == 0) newSamples.back().skipped++;
						}
						else
						{
							int parentIndex = choicesUp[selected];
							withoutReplacementSample& parentSample = samples[parentIndex];
							std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + i * nRows);
							newSampleRowSums[i*nRows + row]--;
							newSamples.push_back(withoutReplacementSample(parentSample.columnSum+1, args.samplingArgs.weights[selected] / inclusionProbability, parentSample.productInclusionProbabilities * inclusionProbability, parentSample.totalRemaining-1));
							newSamples.back().expExponentialParameters = parentSample.expExponentialParameters;
							newSamples.back().expNormalisingConstants = parentSample.expNormalisingConstants;
							newSamples.back().skipped = parentSample.skipped;
							newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
							newSamples.back().nRemainingDeterministic = parentSample.nRemainingDeterministic;
							if(sampleRowSums[nRows * parentIndex + row] == (int)nColumns - column)
							{
								newSamples.back().skipped++;
								newSamples.back().nRemainingDeterministic--;
							}
						}
					}
				}
				samples.swap(newSamples);
				sampleRowSums.swap(newSampleRowSums);
			}
			if(column != (int)nColumns - 1)
			{
				for(std::size_t i = 0; i < samples.size(); i++)
				{
					withoutReplacementSample& currentSample = samples[i];
					currentSample.columnSum = 0;
					//Reset the inclusion probabilities for the conditional Poisson sampling
					cpSamplingArgs.weights.clear();
					cpSamplingArgs.n = initialColumnSums[column+1];
					cpSamplingArgs.weights.insert(cpSamplingArgs.weights.begin(), sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows);
					sampling::samplingBase(cpSamplingArgs.n, cpSamplingArgs.indices, cpSamplingArgs.weights, cpSamplingArgs.rescaledWeights, cpSamplingArgs.zeroWeights, cpSamplingArgs.deterministicInclusion, conditionalPoissonInclusionProbabilities);
					computeExponentialParameters(cpSamplingArgs);
					sampling::conditionalPoissonInclusionProbabilities(cpSamplingArgs, cpSamplingArgs.inclusionProbabilities);

					currentSample.skipped = 0;
					currentSample.nRemainingZeros = currentSample.nRemainingDeterministic = 0;
					for(int j = 0; j < (int)nRows; j++)
					{
						if(sampleRowSums[i * nRows + j] == 0) currentSample.nRemainingZeros++;
						if(cpSamplingArgs.deterministicInclusion[j]) currentSample.nRemainingDeterministic++;
					}
					currentSample.expNormalisingConstants.reset(new boost::numeric::ublas::matrix<mpfr_class>());
					currentSample.expNormalisingConstants->swap(cpSamplingArgs.expNormalisingConstant);

					currentSample.expExponentialParameters.reset(new std::vector<mpfr_class>());
					currentSample.expExponentialParameters->swap(cpSamplingArgs.expExponentialParameters);
				}
			}
			row = 0;
		}
		args.estimate = 0;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			args.estimate += 1/samples[i].productInclusionProbabilities;
		}
	}
}
