#include "conditionalPoissonBootstrapImpl.h"
#include "conditionalPoissonSequential.h"
#include "conditionalPoisson/conditionalPoissonSequential.h"
#include "conditionalPoisson/computeExponentialParameters.h"
#include "conditionalPoisson/calculateExpNormalisingConstants.h"
#include "samplingBase.h"
#include "GayleRyserTest.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
namespace binaryTables
{
	conditionalPoissonBootstrapSample::conditionalPoissonBootstrapSample(conditionalPoissonBootstrapSample&& other)
		: columnSum(other.columnSum), density(std::move(other.density)), expNormalisingConstants(other.expNormalisingConstants), expExponentialParameters(other.expExponentialParameters), skipped(other.skipped), nRemainingZeros(other.nRemainingZeros), nRemainingDeterministic(other.nRemainingDeterministic), deterministicInclusion(std::move(other.deterministicInclusion))
	{
	}
	conditionalPoissonBootstrapSample::conditionalPoissonBootstrapSample(const conditionalPoissonBootstrapSample& other)
		: columnSum(other.columnSum), density(other.density), expNormalisingConstants(other.expNormalisingConstants), expExponentialParameters(other.expExponentialParameters), skipped(other.skipped), nRemainingZeros(other.nRemainingZeros), nRemainingDeterministic(other.nRemainingDeterministic), deterministicInclusion(other.deterministicInclusion)
	{
	}
	conditionalPoissonBootstrapSample& conditionalPoissonBootstrapSample::operator=(conditionalPoissonBootstrapSample&& other)
	{
		columnSum = other.columnSum;
		density = std::move(other.density);
		expNormalisingConstants = other.expNormalisingConstants;
		expExponentialParameters = other.expExponentialParameters;
		skipped = other.skipped;
		nRemainingZeros = other.nRemainingZeros;
		nRemainingDeterministic = other.nRemainingDeterministic;
		deterministicInclusion = std::move(other.deterministicInclusion);
		return *this;
	}
	conditionalPoissonBootstrapSample& conditionalPoissonBootstrapSample::operator=(const conditionalPoissonBootstrapSample& other)
	{
		columnSum = other.columnSum;
		density = other.density;
		expNormalisingConstants = other.expNormalisingConstants;
		expExponentialParameters = other.expExponentialParameters;
		skipped = other.skipped;
		nRemainingZeros = other.nRemainingZeros;
		nRemainingDeterministic = other.nRemainingDeterministic;
		deterministicInclusion = other.deterministicInclusion;
		return *this;
	}
	void conditionalPoissonBootstrap(conditionalPoissonBootstrapArgs& args)
	{
		problem& problemObj = args.problemObj;
		const std::vector<int>& initialRowSums = problemObj.getRowSums();
		const std::vector<int>& initialColumnSums = problemObj.getColumnSums();
		std::size_t nRows = initialRowSums.size(), nColumns = initialColumnSums.size();

		std::size_t n = args.n;

		GayleRyserTestWorking gayleRyserTestWorking(true);

		//Extract data from args
		std::vector<int>& sampleRowSums = args.sampleRowSums;
		std::vector<int>& newSampleRowSums = args.newSampleRowSums;
		std::vector<conditionalPoissonBootstrapSample>& samples = args.samples;
		std::vector<conditionalPoissonBootstrapSample>& newSamples = args.newSamples;
		conditionalPoisson::conditionalPoissonSequentialArgs cpSamplingArgs;
		exponentialPreCompute(cpSamplingArgs.preComputation, nColumns);

		if(nRows <= 1 || nColumns <= 1)
		{
			throw std::runtime_error("The number of columns and rows must be at least 2");
		}
		//Set up data structures for the particles
		sampleRowSums.resize(nRows*n);
		newSampleRowSums.resize(nRows*n);
		samples.clear();
		newSamples.clear();
		//Set up n samples
		{
			int nDeterministic = 0, nZeroWeights = 0;
			cpSamplingArgs.n = initialColumnSums[0];
			conditionalPoisson::conditionalPoissonBase(cpSamplingArgs, initialRowSums.begin(), initialRowSums.end(), nColumns, nDeterministic, nZeroWeights);
			conditionalPoisson::computeExponentialParameters(cpSamplingArgs, nColumns, initialRowSums.begin(), initialRowSums.end());
			conditionalPoisson::calculateExpNormalisingConstants(cpSamplingArgs);

			boost::shared_ptr<std::vector<mpfr_class> > copiedExpExponentialParameters(new std::vector<mpfr_class>());
			(*copiedExpExponentialParameters).swap(cpSamplingArgs.expExponentialParameters);

			boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > copiedExpNormalisingConstant(new boost::numeric::ublas::matrix<mpfr_class>());
			(*copiedExpNormalisingConstant).swap(cpSamplingArgs.expNormalisingConstant);
			//set up the original sample
			conditionalPoissonBootstrapSample initialSample(0, 1, copiedExpNormalisingConstant, copiedExpExponentialParameters);
			initialSample.skipped = 0;
			initialSample.nRemainingDeterministic = initialSample.nRemainingZeros = 0;
			initialSample.deterministicInclusion = cpSamplingArgs.deterministicInclusion;
			for(int i = 0; i < (int)initialRowSums.size(); i++)
			{
				if(initialRowSums[i] == 0) initialSample.nRemainingZeros++;
				if(initialRowSums[i] == (int)nColumns) initialSample.nRemainingDeterministic++;
			}
			//copy the original sample
			samples.insert(samples.begin(), n, initialSample);
			//Make copies of the original row sums
			for(int i = 0; i < (int)n; i++)
			{
				std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin() + i*nRows);
			}
		}
		mpfr_class meanDensity = 1;
		std::vector<double> accumulated;
		std::vector<int> resampleFromIndices;
		for(int column = 0; column < (int)nColumns; column++)
		{
			for(int row = 0; row < (int)nRows; row++)
			{
				//Perform sampling
				mpfr_class sum = 0;
				accumulated.clear();
				resampleFromIndices.clear();
				for(int i = 0; i < (int)n; i++)
				{
					conditionalPoissonBootstrapSample& current = samples[i];
					if(current.deterministicInclusion[row])
					{
						sampleRowSums[i*nRows+row]--;
						current.skipped++;
						current.nRemainingDeterministic--;
						current.columnSum++;
						current.density = 1;
					}
					else if(initialColumnSums[column] == current.columnSum + (int)nRows - row)
					{
						sampleRowSums[i*nRows+row]--;
						current.columnSum++;
						current.density = 1;
					}
					else if(sampleRowSums[i*nRows+row] == 0)
					{
						current.skipped++;
						current.nRemainingZeros--;
						current.density = 1;
					}
					else if(current.columnSum + current.nRemainingDeterministic == initialColumnSums[column])
					{
						current.density = 1;
					}
					else
					{
						mpfr_class selectionProb;
						if(initialColumnSums[column] - current.columnSum - current.nRemainingDeterministic == 1)
						{
							selectionProb = (*current.expExponentialParameters)[row] / (*current.expNormalisingConstants)(row - current.skipped, initialColumnSums[column] - current.columnSum - 1 - current.nRemainingDeterministic);
						}
						else
						{
							selectionProb = (*current.expExponentialParameters)[row] * (*current.expNormalisingConstants)(row+1 - current.skipped, initialColumnSums[column] - current.columnSum - 2 - current.nRemainingDeterministic) / (*current.expNormalisingConstants)(row - current.skipped, initialColumnSums[column] - current.columnSum - 1 - current.nRemainingDeterministic);
						}
						boost::random::bernoulli_distribution<> randomChoice(selectionProb.convert_to<double>());
						if(randomChoice(args.randomSource))
						{
							current.density = selectionProb;
							sampleRowSums[i*nRows+row]--;
							current.columnSum++;
						}
						else
						{
							current.density = (1 - selectionProb);
						}
					}
					//If we've just finishing simulating a column, then the resampling should use the gale-ryser test. 
					if(row == (int)nRows - 1)
					{
						gayleRyserTestWorking.sortedSums1.clear();
						gayleRyserTestWorking.sortedSums2.clear();
						gayleRyserTestWorking.sortedSums1.insert(gayleRyserTestWorking.sortedSums1.begin(), sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows);
						gayleRyserTestWorking.sortedSums2.insert(gayleRyserTestWorking.sortedSums2.begin(), initialColumnSums.begin(), initialColumnSums.end());
						if(GayleRyserTest(gayleRyserTestWorking.sortedSums2, gayleRyserTestWorking.sortedSums1, column+1, gayleRyserTestWorking))
						{
							sum += 1/current.density;
							accumulated.push_back(sum.convert_to<double>());
							resampleFromIndices.push_back(i);
						}
					}
					else
					{
						sum += 1/current.density;
						accumulated.push_back(sum.convert_to<double>());
						resampleFromIndices.push_back(i);
					}
				}
				meanDensity *= sum/n;
				//Resample
				newSamples.clear();
				boost::random::uniform_real_distribution<> resampler(0, sum.convert_to<double>());
				for(int i = 0; i < (int)n; i++)
				{
					double value = resampler(args.randomSource);
					int index = std::distance(accumulated.begin(), std::lower_bound(accumulated.begin(), accumulated.end(), value));
					index = resampleFromIndices[index];
					newSamples.push_back(samples[index]);
					std::copy(sampleRowSums.begin() + index * nRows, sampleRowSums.begin() + (index + 1) * nRows, newSampleRowSums.begin() + i * nRows);
				}
				samples.swap(newSamples);
				sampleRowSums.swap(newSampleRowSums);
			}
			//Reset column sums to zero and update conditional poisson data. 
			for(int i = 0; i < (int)n; i++)
			{
				conditionalPoissonBootstrapSample& current = samples[i];
				current.columnSum = 0;
				current.nRemainingDeterministic = current.nRemainingZeros = 0;
				current.skipped = 0;
				cpSamplingArgs.n = initialColumnSums[column+1];
				conditionalPoisson::conditionalPoissonBase(cpSamplingArgs, sampleRowSums.begin() + i * nRows, sampleRowSums.begin() + (i + 1) * nRows, nColumns-column-1, current.nRemainingDeterministic, current.nRemainingZeros);
				conditionalPoisson::computeExponentialParameters(cpSamplingArgs, nColumns - column - 1, sampleRowSums.begin() + i * nRows, sampleRowSums.begin() + (i + 1) * nRows);
				conditionalPoisson::calculateExpNormalisingConstants(cpSamplingArgs);

				current.expNormalisingConstants.reset(new boost::numeric::ublas::matrix<mpfr_class>());
				current.expNormalisingConstants->swap(cpSamplingArgs.expNormalisingConstant);

				current.expExponentialParameters.reset(new std::vector<mpfr_class>());
				current.expExponentialParameters->swap(cpSamplingArgs.expExponentialParameters);

				current.deterministicInclusion = cpSamplingArgs.deterministicInclusion;
			}
		}
		args.estimate = meanDensity;
	}
}
