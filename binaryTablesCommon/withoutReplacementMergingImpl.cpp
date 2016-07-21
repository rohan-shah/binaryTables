#include "withoutReplacementMergingImpl.h"
#include "conditionalPoissonSequential.h"
#include "samplingBase.h"
#include "GayleRyserTest.h"
#include "conditionalPoisson/conditionalPoissonSequential.h"
#include "conditionalPoisson/computeExponentialParameters.h"
#include "conditionalPoisson/calculateExpNormalisingConstants.h"
namespace binaryTables
{
	withoutReplacementMergingSample::withoutReplacementMergingSample(withoutReplacementMergingSample&& other)
		: columnSum(other.columnSum), sizeVariable(std::move(other.sizeVariable)), weight(std::move(other.weight)), totalRemaining(other.totalRemaining), expNormalisingConstants(other.expNormalisingConstants), expExponentialParameters(other.expExponentialParameters), skipped(other.skipped), nRemainingZeros(other.nRemainingZeros), nRemainingDeterministic(other.nRemainingDeterministic), deterministicInclusion(std::move(other.deterministicInclusion))
	{
	}
	withoutReplacementMergingSample& withoutReplacementMergingSample::operator=(withoutReplacementMergingSample&& other)
	{
		columnSum = other.columnSum;
		sizeVariable = std::move(other.sizeVariable);
		weight = std::move(other.weight);
		totalRemaining = other.totalRemaining;
		expNormalisingConstants = other.expNormalisingConstants;
		expExponentialParameters = other.expExponentialParameters;
		skipped = other.skipped;
		nRemainingZeros = other.nRemainingZeros;
		nRemainingDeterministic = other.nRemainingDeterministic;
		deterministicInclusion = std::move(other.deterministicInclusion);
		return *this;
	}
	void withoutReplacementMerging(withoutReplacementMergingArgs& args)
	{
		if(args.mergeFrequency <= 0)
		{
			throw std::runtime_error("Input mergeFrequency must be at least 1");
		}
		problem& problemObj = args.problemObj;
		const std::vector<int>& initialRowSums = problemObj.getRowSums();
		const std::vector<int>& initialColumnSums = problemObj.getColumnSums();
		std::size_t nRows = initialRowSums.size(), nColumns = initialColumnSums.size();

		//Count the total number of ones
		int totalOnes = 0;
		for(std::vector<int>::const_iterator i = initialRowSums.begin(); i != initialRowSums.end(); i++) totalOnes += *i;
		
		std::size_t n = args.n;

		GayleRyserTestWorking gayleRyserTestWorking(true);

		//Extract data from args
		std::vector<int>& sampleRowSums = args.sampleRowSums;
		std::vector<int>& newSampleRowSums = args.newSampleRowSums;
		std::vector<int> sortedSampleRowSums(nRows*std::max((std::size_t)2, n));
		std::vector<withoutReplacementMergingSample>& samples = args.samples;
		std::vector<withoutReplacementMergingSample>& newSamples = args.newSamples;
		std::vector<mpfr_class> conditionalPoissonInclusionProbabilities;
		conditionalPoisson::conditionalPoissonSequentialArgs cpSamplingArgs;
		cpSamplingArgs.n = initialColumnSums[0];
		conditionalPoisson::exponentialPreCompute(cpSamplingArgs.preComputation, nColumns);

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
		int nDeterministic = 0, nZeroWeights = 0;
		conditionalPoisson::conditionalPoissonBase(cpSamplingArgs, initialRowSums.begin(), initialRowSums.end(), nColumns, nDeterministic, nZeroWeights);

		conditionalPoisson::computeExponentialParameters(cpSamplingArgs, nColumns, initialRowSums.begin(), initialRowSums.end());
		conditionalPoisson::calculateExpNormalisingConstants(cpSamplingArgs);
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
			samples.push_back(withoutReplacementMergingSample(0, 1 - selectionProb, 1, totalOnes, initialExpNormalisingConstantData, initialExpExponentialParameters));
			samples[0].skipped = 1;
		}
		else if(initialRowSums[0] == (int)nColumns)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			samples.push_back(withoutReplacementMergingSample(1, selectionProb, 1, totalOnes-1, initialExpNormalisingConstantData, initialExpExponentialParameters));
			sampleRowSums[0]--;
			samples[0].skipped = 1;
		}
		else
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin() + nRows);
			sampleRowSums[0]--;
			//The first sample has a one in the first column. To start off with, don't 
			samples.push_back(withoutReplacementMergingSample(1, selectionProb, 1, totalOnes-1, initialExpNormalisingConstantData, initialExpExponentialParameters));
			samples.push_back(withoutReplacementMergingSample(0, 1 - selectionProb, 1, totalOnes, initialExpNormalisingConstantData, initialExpExponentialParameters));
			samples[0].skipped = samples[1].skipped = 0;
		}
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			//This is to do with skipping over the rows which have a row-sum of zero
			samples[i].nRemainingZeros = samples[i].nRemainingDeterministic = 0;
			samples[i].deterministicInclusion = cpSamplingArgs.deterministicInclusion;
			for(int j = 0; j < (int)nRows; j++)
			{
				if(sampleRowSums[i * nRows + j] == 0) samples[i].nRemainingZeros++;
				if(sampleRowSums[i * nRows + j] == (int)nColumns - (j == 0)) samples[i].nRemainingDeterministic++;
			}
		}

		int doneSteps = 1;
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
						if(samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column] && samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column]) choicesDown.push_back((int)i);
					}
					else if(sampleRowSums[i*nRows + row] == (int)nColumns - column)
					{
						if(samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column]) choicesUp.push_back((int)i);
					}
					else
					{
						if(samples[i].deterministicInclusion[row])
						{
							if(samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column]) choicesUp.push_back((int)i);
						}
						else
						{
							if(samples[i].columnSum + samples[i].nRemainingDeterministic + 1 <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column]) choicesUp.push_back((int)i);
							if(samples[i].columnSum + (int)nRows - row - 1 - samples[i].nRemainingZeros >= initialColumnSums[column] && samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column]) choicesDown.push_back((int)i);
						}
					}
				}
				newSamples.clear();
				if(choicesUp.size() + choicesDown.size() <= n)
				{
					std::size_t outputCounter = 0;
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						int parentIndex = choicesUp[i];
						withoutReplacementMergingSample& parentSample = samples[choicesUp[i]];
						std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + outputCounter * nRows);
						int newSkipped = parentSample.skipped, newDeterministic = parentSample.nRemainingDeterministic;
						if(sampleRowSums[nRows * parentIndex + row] == (int)nColumns - column || parentSample.deterministicInclusion[row])
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
						newSamples.push_back(withoutReplacementMergingSample(parentSample.columnSum + 1, parentSample.sizeVariable * selectionProb, parentSample.weight, parentSample.totalRemaining - 1, parentSample.expNormalisingConstants, parentSample.expExponentialParameters));
						newSamples.back().skipped = newSkipped;
						newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
						newSamples.back().nRemainingDeterministic = newDeterministic;
						newSamples.back().deterministicInclusion = parentSample.deterministicInclusion;
						outputCounter++;
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						int parentIndex = choicesDown[i];
						withoutReplacementMergingSample& parentSample = samples[parentIndex];
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
						newSamples.push_back(withoutReplacementMergingSample(parentSample.columnSum, parentSample.sizeVariable * (1 - selectionProb), parentSample.weight, parentSample.totalRemaining, parentSample.expNormalisingConstants, parentSample.expExponentialParameters));
						newSamples.back().skipped = newSkipped;
						newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
						newSamples.back().nRemainingDeterministic = parentSample.nRemainingDeterministic;
						if(newSampleRowSums[outputCounter * nRows + row] == 0) newSamples.back().nRemainingZeros--;
						newSamples.back().deterministicInclusion = parentSample.deterministicInclusion;
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
						withoutReplacementMergingSample& parentSample = samples[parentIndex];
						if(sampleRowSums[nRows * parentIndex + row] == (int)nColumns - column)
						{
							selectionProb = 1;
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
						args.samplingArgs.weights.push_back(parentSample.sizeVariable * selectionProb);
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						int parentIndex = choicesDown[i];
						withoutReplacementMergingSample& parentSample = samples[parentIndex];
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
							withoutReplacementMergingSample& parentSample = samples[parentIndex];
							std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + i * nRows);
							newSamples.push_back(withoutReplacementMergingSample(parentSample.columnSum, args.samplingArgs.weights[selected] / inclusionProbability, parentSample.weight / inclusionProbability, parentSample.totalRemaining, parentSample.expNormalisingConstants, parentSample.expExponentialParameters));
							newSamples.back().skipped = parentSample.skipped;
							newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
							newSamples.back().nRemainingDeterministic = parentSample.nRemainingDeterministic;
							if(sampleRowSums[parentIndex*nRows + row] == 0) newSamples.back().nRemainingZeros--;
							if(parentSample.columnSum == initialColumnSums[column] || sampleRowSums[parentIndex*nRows + row] == 0) newSamples.back().skipped++;
							newSamples.back().deterministicInclusion = parentSample.deterministicInclusion;
						}
						else
						{
							int parentIndex = choicesUp[selected];
							withoutReplacementMergingSample& parentSample = samples[parentIndex];
							std::copy(sampleRowSums.begin() + nRows * parentIndex, sampleRowSums.begin() + nRows * (parentIndex+1), newSampleRowSums.begin() + i * nRows);
							newSampleRowSums[i*nRows + row]--;
							newSamples.push_back(withoutReplacementMergingSample(parentSample.columnSum+1, args.samplingArgs.weights[selected] / inclusionProbability, parentSample.weight / inclusionProbability, parentSample.totalRemaining-1, parentSample.expNormalisingConstants, parentSample.expExponentialParameters));
							newSamples.back().skipped = parentSample.skipped;
							newSamples.back().nRemainingZeros = parentSample.nRemainingZeros;
							newSamples.back().nRemainingDeterministic = parentSample.nRemainingDeterministic;
							if(sampleRowSums[nRows * parentIndex + row] == (int)nColumns - column || parentSample.deterministicInclusion[row])
							{
								newSamples.back().skipped++;
								newSamples.back().nRemainingDeterministic--;
							}
							newSamples.back().deterministicInclusion = parentSample.deterministicInclusion;
						}
					}
				}
				samples.swap(newSamples);
				sampleRowSums.swap(newSampleRowSums);
				doneSteps++;
				if((doneSteps % args.mergeFrequency) == 0)
				{
					std::vector<bool> alreadySelected(samples.size(), false);
					newSamples.clear();
					//sort the sample row sums, to help us merge different samples. We need to make copies so we avoid screwing up the link between the poisson sampling data and the row sums. 
					std::copy(sampleRowSums.begin(), sampleRowSums.end(), sortedSampleRowSums.begin());
					for(int i = 0; i < (int)samples.size(); i++)
					{
						std::sort(sortedSampleRowSums.begin() + i * nRows, sortedSampleRowSums.begin() + i * nRows + row + 1);
						std::sort(sortedSampleRowSums.begin() + i * nRows + row + 1, sortedSampleRowSums.begin() + (i + 1) * nRows);
					}
					for(int i = 0; i < (int)samples.size(); i++)
					{
						if(!alreadySelected[i])
						{
							for(int j = i+1; j < (int)samples.size(); j++)
							{
								if(samples[i].columnSum == samples[j].columnSum && memcmp(&(sortedSampleRowSums[i*nRows]), &(sortedSampleRowSums[j*nRows]), sizeof(int)*nRows) == 0)
								{
									samples[i].weight += samples[j].weight;
									samples[i].sizeVariable += samples[j].sizeVariable;
									alreadySelected[j] = true;
								}
							}
							std::copy(sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows, newSampleRowSums.begin() + newSamples.size()*nRows);
							newSamples.push_back(std::move(samples[i]));
						}
					}
					samples.clear();
					newSamples.swap(samples);
					sampleRowSums.swap(newSampleRowSums);
				}
			}
			if(column != (int)nColumns - 1)
			{
				for(std::size_t i = 0; i < samples.size(); i++)
				{
					bool shouldRemove = false;
					//This throws an exception if there is no possible sample, in which case the sample is removed. 
					try
					{
						withoutReplacementMergingSample& currentSample = samples[i];

						//Use the Gayle Ryser test to check if we should continue
						gayleRyserTestWorking.sortedSums1.clear();
						gayleRyserTestWorking.sortedSums2.clear();
						gayleRyserTestWorking.sortedSums1.insert(gayleRyserTestWorking.sortedSums1.begin(), sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows);
						gayleRyserTestWorking.sortedSums2.insert(gayleRyserTestWorking.sortedSums2.begin(), initialColumnSums.begin(), initialColumnSums.end());
						if(!GayleRyserTest(gayleRyserTestWorking.sortedSums2, gayleRyserTestWorking.sortedSums1, column+1, gayleRyserTestWorking))
						{
							shouldRemove = true;
							goto testRemove;
						}
						currentSample.columnSum = 0;
						//Reset the inclusion probabilities for the conditional Poisson sampling
						cpSamplingArgs.n = initialColumnSums[column+1];
						currentSample.nRemainingDeterministic = 0;
						currentSample.nRemainingZeros = 0;
						conditionalPoisson::conditionalPoissonBase(cpSamplingArgs, sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows, nColumns - (column + 1), currentSample.nRemainingDeterministic, currentSample.nRemainingZeros);
						conditionalPoisson::computeExponentialParameters(cpSamplingArgs, nColumns - (column + 1), sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows);
						conditionalPoisson::calculateExpNormalisingConstants(cpSamplingArgs);

						currentSample.skipped = 0;
						currentSample.deterministicInclusion = cpSamplingArgs.deterministicInclusion;
						//Swap in data from cpSamplingArgs to the sample
						currentSample.expNormalisingConstants.reset(new boost::numeric::ublas::matrix<mpfr_class>());
						currentSample.expNormalisingConstants->swap(cpSamplingArgs.expNormalisingConstant);

						currentSample.expExponentialParameters.reset(new std::vector<mpfr_class>());
						currentSample.expExponentialParameters->swap(cpSamplingArgs.expExponentialParameters);
					}
					catch(...)
					{
						shouldRemove = true;
					}
testRemove:
					if(shouldRemove)
					{
						if(i != samples.size() - 1)
						{
							samples[i] = std::move(samples.back());
							std::copy(sampleRowSums.begin() + (samples.size()-1)*nRows, sampleRowSums.begin() + samples.size() * nRows, sampleRowSums.begin() + i*nRows);
							i--;
						}
						samples.pop_back();
					}
				}
			}
			row = 0;
		}
		args.estimate = 0;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			args.estimate += samples[i].weight;
		}
	}
}
