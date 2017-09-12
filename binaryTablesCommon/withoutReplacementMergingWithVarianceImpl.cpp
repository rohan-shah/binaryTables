#include "withoutReplacementMergingWithVarianceImpl.h"
#include "conditionalPoissonSequential.h"
#include "samplingBase.h"
#include "GayleRyserTest.h"
#include "conditionalPoisson/conditionalPoissonSequential.h"
#include "conditionalPoisson/computeExponentialParameters.h"
#include "conditionalPoisson/calculateExpNormalisingConstants.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/graphviz.hpp>
namespace binaryTables
{
	withoutReplacementMergingWithVarianceSample::withoutReplacementMergingWithVarianceSample(withoutReplacementMergingWithVarianceSample&& other)
		: columnSum(other.columnSum), trueDensity(other.trueDensity), importanceDensity(std::move(other.importanceDensity)), weight(std::move(other.weight)), totalRemaining(other.totalRemaining), expNormalisingConstants(other.expNormalisingConstants), expExponentialParameters(other.expExponentialParameters), skipped(other.skipped), nRemainingZeros(other.nRemainingZeros), nRemainingDeterministic(other.nRemainingDeterministic), parentIndexWithinDesign(other.parentIndexWithinDesign), indexWithinDesign(other.indexWithinDesign), deterministicInclusion(std::move(other.deterministicInclusion)) 
	{
	}
	withoutReplacementMergingWithVarianceSample& withoutReplacementMergingWithVarianceSample::operator=(withoutReplacementMergingWithVarianceSample&& other)
	{
		columnSum = other.columnSum;
		trueDensity = other.trueDensity;
		importanceDensity = std::move(other.importanceDensity);
		weight = std::move(other.weight);
		totalRemaining = other.totalRemaining;
		expNormalisingConstants = other.expNormalisingConstants;
		expExponentialParameters = other.expExponentialParameters;
		skipped = other.skipped;
		nRemainingZeros = other.nRemainingZeros;
		nRemainingDeterministic = other.nRemainingDeterministic;
		parentIndexWithinDesign = other.parentIndexWithinDesign;
		indexWithinDesign = other.indexWithinDesign;
		deterministicInclusion = std::move(other.deterministicInclusion);
		return *this;
	}
	struct varianceGraphVertex
	{
	public:
		varianceGraphVertex()
			: indexWithinDesign(-1), samplingStage(-1), indexWithinSelected(-1), trueDensity(-1), accumulatedMean(0), V(0)
		{}
		int indexWithinDesign;
		int samplingStage;
		int indexWithinSelected;
		double trueDensity;
		mutable mpfr_class accumulatedMean, V;
	};
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_name_t, varianceGraphVertex> > varianceGraph;
	class vertexPropertyWriter
	{
	public:
		vertexPropertyWriter(const varianceGraph& g)
			:g(g)
		{}
		template<typename vertexDesc> void operator()(std::ostream& out, const vertexDesc& v)
		{
			const varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, g, v);
			out << std::endl << "[trueDensity=\"" << vertexInfo.trueDensity << "\"";
			out << ",accumulatedMean=\"" << vertexInfo.accumulatedMean.convert_to<double>() << "\"";
			out << ",indexWithinDesign=\"" << vertexInfo.indexWithinDesign << "\"";
			out << ",samplingStage=\"" << vertexInfo.samplingStage << "\"";
			out << ",V=\"" << vertexInfo.V.convert_to<double>() << "\"";
			out << ",indexWithinSelected=\"" << vertexInfo.indexWithinSelected <<"\"]" << std::endl;
		}
		const varianceGraph& g;
	};
	void withoutReplacementMergingWithVariance(withoutReplacementMergingWithVarianceArgs& args)
	{
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
		std::vector<int> sortedSampleRowSums(2 * nRows*std::max((std::size_t)2, n));
		std::vector<withoutReplacementMergingWithVarianceSample>& samples = args.samples;
		std::vector<withoutReplacementMergingWithVarianceSample>& newSamples = args.newSamples;
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
		sampleRowSums.resize(2 * nRows*std::max((std::size_t)2, n));
		newSampleRowSums.resize(2*nRows*std::max((std::size_t)2, n));
		samples.clear();
		newSamples.clear();

		std::vector<std::vector<::sampling::mpfr_class> > allInclusionProbabilities;
		std::vector<boost::numeric::ublas::matrix<sampling::mpfr_class> > allSecondOrderInclusionProbabilities;

		varianceGraph varianceEstimationGraph;
		std::vector<std::vector<int> > graphVertices(nRows * nColumns);
		varianceGraph::vertex_descriptor rootVertex = boost::add_vertex(varianceEstimationGraph);
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

			withoutReplacementMergingWithVarianceSample newParticle;
			newParticle.columnSum = 0;
			newParticle.skipped = 1;
			newParticle.importanceDensity = 1 - selectionProb;
			newParticle.weight = 1;
			newParticle.trueDensity = 0.5;
			newParticle.totalRemaining = totalOnes;
			newParticle.expNormalisingConstants = initialExpNormalisingConstantData;
			newParticle.expExponentialParameters = initialExpExponentialParameters;
			newParticle.parentIndexWithinDesign = -1;
			newParticle.indexWithinDesign = 0;
			samples.emplace_back(std::move(newParticle));

			varianceGraph::vertex_descriptor firstVertex = boost::add_vertex(varianceEstimationGraph);
			varianceGraphVertex& firstVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, firstVertex);
			firstVertexInfo.samplingStage = 0;
			firstVertexInfo.indexWithinDesign = firstVertexInfo.indexWithinSelected = 0;
			firstVertexInfo.trueDensity = 0.5;
			graphVertices[0].push_back((int)firstVertex);
			boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);

			//Specify inclusion probabilities and second order inclusion probabilities
			{
				std::vector<::sampling::mpfr_class> inclusion(1, 1.0);
				allInclusionProbabilities.emplace_back(std::move(inclusion));
				boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(1, 1, 1);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
			}
		}
		//In this case we *must* have a 1 in the first entry
		else if(initialRowSums[0] == (int)nColumns)
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			withoutReplacementMergingWithVarianceSample newParticle;
			newParticle.columnSum = 1;
			newParticle.skipped = 1;
			newParticle.importanceDensity = selectionProb;
			newParticle.weight = 1;
			newParticle.trueDensity = 0.5;
			newParticle.totalRemaining = totalOnes - 1;
			newParticle.expNormalisingConstants = initialExpNormalisingConstantData;
			newParticle.expExponentialParameters = initialExpExponentialParameters;
			newParticle.parentIndexWithinDesign = -1;
			newParticle.indexWithinDesign = 0;
			samples.emplace_back(std::move(newParticle));
			sampleRowSums[0]--;

			varianceGraph::vertex_descriptor firstVertex = boost::add_vertex(varianceEstimationGraph);
			varianceGraphVertex& firstVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, firstVertex);
			firstVertexInfo.samplingStage = 0;
			firstVertexInfo.indexWithinDesign = firstVertexInfo.indexWithinSelected = 0;
			firstVertexInfo.trueDensity = 0.5;
			graphVertices[0].push_back((int)firstVertex);
			boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);

			//Specify inclusion probabilities and second order inclusion probabilities
			{
				std::vector<::sampling::mpfr_class> inclusion(1, 1.0);
				allInclusionProbabilities.emplace_back(std::move(inclusion));
				boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(1, 1, 1);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
			}
		}
		else
		{
			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin());
			sampleRowSums[0]--;
			//The first sample has a one in the first column.
			withoutReplacementMergingWithVarianceSample firstParticle;
			firstParticle.columnSum = 1;
			firstParticle.skipped = 0;
			firstParticle.importanceDensity = selectionProb;
			firstParticle.weight = 1;
			firstParticle.trueDensity = 0.5;
			firstParticle.totalRemaining = totalOnes - 1;
			firstParticle.expNormalisingConstants = initialExpNormalisingConstantData;
			firstParticle.expExponentialParameters = initialExpExponentialParameters;
			firstParticle.parentIndexWithinDesign = -1;
			firstParticle.indexWithinDesign = 0;
			samples.emplace_back(std::move(firstParticle));

			std::copy(initialRowSums.begin(), initialRowSums.end(), sampleRowSums.begin() + nRows);
			withoutReplacementMergingWithVarianceSample secondParticle;
			secondParticle.columnSum = 0;
			secondParticle.skipped = 0;
			secondParticle.importanceDensity = 1 - selectionProb;
			secondParticle.weight = 1;
			secondParticle.trueDensity = 0.5;
			secondParticle.totalRemaining = totalOnes;
			secondParticle.expNormalisingConstants = initialExpNormalisingConstantData;
			secondParticle.expExponentialParameters = initialExpExponentialParameters;
			secondParticle.parentIndexWithinDesign = -1;
			secondParticle.indexWithinDesign = 1;
			samples.emplace_back(std::move(secondParticle));

			varianceGraph::vertex_descriptor firstVertex = boost::add_vertex(varianceEstimationGraph);
			varianceGraphVertex& firstVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, firstVertex);
			firstVertexInfo.samplingStage = 0;
			firstVertexInfo.indexWithinDesign = firstVertexInfo.indexWithinSelected = 0;
			firstVertexInfo.trueDensity = 0.5;
			graphVertices[0].push_back((int)firstVertex);
			boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);

			varianceGraph::vertex_descriptor secondVertex = boost::add_vertex(varianceEstimationGraph);
			varianceGraphVertex& secondVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, secondVertex);
			secondVertexInfo.samplingStage = 0;
			secondVertexInfo.indexWithinDesign = secondVertexInfo.indexWithinSelected = 1;
			secondVertexInfo.trueDensity = 0.5;
			graphVertices[0].push_back((int)secondVertex);
			boost::add_edge(rootVertex, secondVertex, varianceEstimationGraph);

			//Specify inclusion probabilities and second order inclusion probabilities
			{
				std::vector<::sampling::mpfr_class> inclusion(2, 1.0);
				allInclusionProbabilities.emplace_back(std::move(inclusion));
				boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(2, 2, 1);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
			}
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
				std::vector<int>& currentGraphVertices = graphVertices[doneSteps];
				const std::vector<int>& previousGraphVertices = graphVertices[doneSteps - 1];
				//organise the vector of possibilities (choicesUp, choicesDown) for the next step. 
				for(std::size_t i = 0; i < samples.size(); i++)
				{
					bool upAllowed = true, downAllowed = true;
					if(row == (int)nRows - 1 && column != (int)nColumns - 1)
					{
						//Use the Gayle Ryser test to check if we should continue
						gayleRyserTestWorking.sortedSums1.clear();
						gayleRyserTestWorking.sortedSums2.clear();
						gayleRyserTestWorking.sortedSums1.insert(gayleRyserTestWorking.sortedSums1.begin(), sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows);
						gayleRyserTestWorking.sortedSums2.insert(gayleRyserTestWorking.sortedSums2.begin(), initialColumnSums.begin(), initialColumnSums.end());
						downAllowed = GayleRyserTest(gayleRyserTestWorking.sortedSums2, gayleRyserTestWorking.sortedSums1, column+1, gayleRyserTestWorking);

						if(sampleRowSums[i*nRows + row] > 0)
						{
							gayleRyserTestWorking.sortedSums1.clear();
							gayleRyserTestWorking.sortedSums2.clear();
							gayleRyserTestWorking.sortedSums1.insert(gayleRyserTestWorking.sortedSums1.begin(), sampleRowSums.begin() + i*nRows, sampleRowSums.begin() + (i+1)*nRows);
							gayleRyserTestWorking.sortedSums2.insert(gayleRyserTestWorking.sortedSums2.begin(), initialColumnSums.begin(), initialColumnSums.end());
							gayleRyserTestWorking.sortedSums1[row]--;
							upAllowed = GayleRyserTest(gayleRyserTestWorking.sortedSums2, gayleRyserTestWorking.sortedSums1, column+1, gayleRyserTestWorking);
						}
					}
					if(sampleRowSums[i*nRows + row] == 0)
					{
						if(downAllowed && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column] && samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column]) choicesDown.push_back((int)i);
					}
					else if(sampleRowSums[i*nRows + row] == (int)nColumns - column)
					{
						if(upAllowed && samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column]) choicesUp.push_back((int)i);
					}
					else
					{
						if(samples[i].deterministicInclusion[row])
						{
							if(upAllowed && samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column]) choicesUp.push_back((int)i);
						}
						else
						{
							if(upAllowed && samples[i].columnSum + samples[i].nRemainingDeterministic + 1 <= initialColumnSums[column] && samples[i].columnSum + (int)nRows - row - samples[i].nRemainingZeros >= initialColumnSums[column]) choicesUp.push_back((int)i);
							if(downAllowed && samples[i].columnSum + (int)nRows - row - 1 - samples[i].nRemainingZeros >= initialColumnSums[column] && samples[i].columnSum + samples[i].nRemainingDeterministic <= initialColumnSums[column]) choicesDown.push_back((int)i);
						}
					}
				}
				//Expand choicesUp and choicesDown into an actual set of particles
				{
					newSamples.clear();
					std::size_t outputCounter = 0;
					for(std::size_t i = 0; i < choicesUp.size(); i++)
					{
						int parentIndex = choicesUp[i];
						withoutReplacementMergingWithVarianceSample& parentSample = samples[choicesUp[i]];
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
						withoutReplacementMergingWithVarianceSample currentNewSample;
						currentNewSample.columnSum = parentSample.columnSum + 1;
						currentNewSample.trueDensity = parentSample.trueDensity;
						currentNewSample.importanceDensity = parentSample.importanceDensity * selectionProb;
						currentNewSample.weight = parentSample.weight;
						currentNewSample.totalRemaining = parentSample.totalRemaining - 1;
						currentNewSample.expExponentialParameters = parentSample.expExponentialParameters;
						currentNewSample.expNormalisingConstants = parentSample.expNormalisingConstants;
						currentNewSample.skipped = newSkipped;
						currentNewSample.nRemainingZeros = parentSample.nRemainingZeros;
						currentNewSample.nRemainingDeterministic = newDeterministic;
						currentNewSample.parentIndexWithinDesign = parentSample.indexWithinDesign;
						currentNewSample.deterministicInclusion = parentSample.deterministicInclusion;
						newSamples.emplace_back(std::move(currentNewSample));
						outputCounter++;
					}
					for(std::size_t i = 0; i < choicesDown.size(); i++)
					{
						int parentIndex = choicesDown[i];
						withoutReplacementMergingWithVarianceSample& parentSample = samples[parentIndex];
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
						withoutReplacementMergingWithVarianceSample currentNewSample;
						currentNewSample.columnSum = parentSample.columnSum;
						currentNewSample.trueDensity = parentSample.trueDensity;
						currentNewSample.importanceDensity = parentSample.importanceDensity * (1 - selectionProb);
						currentNewSample.weight = parentSample.weight;
						currentNewSample.totalRemaining = parentSample.totalRemaining - 1;
						currentNewSample.expExponentialParameters = parentSample.expExponentialParameters;
						currentNewSample.expNormalisingConstants = parentSample.expNormalisingConstants;
						currentNewSample.skipped = newSkipped;
						currentNewSample.nRemainingZeros = parentSample.nRemainingZeros;
						currentNewSample.nRemainingDeterministic = parentSample.nRemainingDeterministic;
						currentNewSample.parentIndexWithinDesign = parentSample.indexWithinDesign;
						if(newSampleRowSums[outputCounter * nRows + row] == 0) currentNewSample.nRemainingZeros--;
						currentNewSample.deterministicInclusion = parentSample.deterministicInclusion;
						newSamples.emplace_back(std::move(currentNewSample));
						outputCounter++;
					}
				}
				//Particle merging. 
				{
					std::vector<bool> alreadySelected(newSamples.size(), false);
					samples.clear();
					//sort the sample row sums, to help us merge different samples. We need to make copies so we avoid screwing up the link between the poisson sampling data and the row sums. 
					std::copy(newSampleRowSums.begin(), newSampleRowSums.end(), sortedSampleRowSums.begin());
					for(int i = 0; i < (int)newSamples.size(); i++)
					{
						std::sort(sortedSampleRowSums.begin() + i * nRows, sortedSampleRowSums.begin() + i * nRows + row + 1);
						std::sort(sortedSampleRowSums.begin() + i * nRows + row + 1, sortedSampleRowSums.begin() + (i + 1) * nRows);
					}
					for(int i = 0; i < (int)newSamples.size(); i++)
					{
						if(!alreadySelected[i])
						{
							varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
							for(int j = i+1; j < (int)newSamples.size(); j++)
							{
								if(newSamples[i].columnSum == newSamples[j].columnSum && memcmp(&(sortedSampleRowSums[i*nRows]), &(sortedSampleRowSums[j*nRows]), sizeof(int)*nRows) == 0)
								{
									newSamples[i].weight += newSamples[j].weight;
									newSamples[i].importanceDensity += newSamples[j].importanceDensity;
									newSamples[i].trueDensity += newSamples[j].trueDensity;
									alreadySelected[j] = true;
									boost::add_edge(previousGraphVertices[newSamples[j].parentIndexWithinDesign], newVertex, varianceEstimationGraph);
								}
							}
							std::copy(newSampleRowSums.begin() + i*nRows, newSampleRowSums.begin() + (i+1)*nRows, sampleRowSums.begin() + samples.size()*nRows);

							varianceGraphVertex& newVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
							newVertexInfo.samplingStage = doneSteps;
							newVertexInfo.indexWithinDesign = (int)samples.size();
							newVertexInfo.indexWithinSelected = -1;
							newVertexInfo.trueDensity = newSamples[i].trueDensity;
							currentGraphVertices.push_back((int)newVertex);
							boost::add_edge(previousGraphVertices[newSamples[i].parentIndexWithinDesign], newVertex, varianceEstimationGraph);
							samples.push_back(std::move(newSamples[i]));
						}
					}
				}
				if(samples.size() <= n)
				{
					for(int i = 0; i < (int)samples.size(); i++)
					{
						int vertex = currentGraphVertices[i];
						varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, vertex);
						vertexInfo.indexWithinSelected = i;
						samples[i].indexWithinDesign = i;
					}
					//Specify inclusion probabilities and second order inclusion probabilities
					{
						std::vector<::sampling::mpfr_class> inclusion(samples.size(), 1.0);
						allInclusionProbabilities.emplace_back(std::move(inclusion));
						boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(samples.size(), samples.size(), 1);
						allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
					}
				}
				else
				{
					args.samplingArgs.weights.clear();
					args.samplingArgs.n = n;
					newSamples.clear();
					mpfr_class sumImportanceDensities = 0;
					for(std::size_t i = 0; i < samples.size(); i++) sumImportanceDensities += samples[i].importanceDensity;
					for(std::size_t i = 0; i < samples.size(); i++)
					{
						args.samplingArgs.weights.push_back(boost::multiprecision::min(n * samples[i].importanceDensity / sumImportanceDensities, mpfr_class(1)));
					}
					conditionalPoissonSequential(args.samplingArgs, args.randomSource);
					for(std::size_t i = 0; i < n; i++)
					{
						int selected = args.samplingArgs.indices[i];
						mpfr_class& inclusionProbability = args.samplingArgs.inclusionProbabilities[selected];

						withoutReplacementMergingWithVarianceSample& parentSample = samples[selected];
						std::copy(sampleRowSums.begin() + nRows * selected, sampleRowSums.begin() + nRows * (selected + 1), newSampleRowSums.begin() + i * nRows);
						newSamples.emplace_back(std::move(parentSample));
						newSamples.back().weight /= inclusionProbability;
						newSamples.back().importanceDensity /= inclusionProbability;
						newSamples.back().indexWithinDesign = selected;

						int vertex = currentGraphVertices[selected];
						varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, vertex);
						vertexInfo.indexWithinSelected = i;
					}
					boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(samples.size(), samples.size(), 0);
					sampling::conditionalPoissonSecondOrderInclusionProbabilities(args.samplingArgs, args.samplingArgs.inclusionProbabilities, secondOrder);
					allInclusionProbabilities.emplace_back(std::move(args.samplingArgs.inclusionProbabilities));
					allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));

					newSamples.swap(samples);
					newSampleRowSums.swap(sampleRowSums);
				}
				doneSteps++;
			}
			if(column != (int)nColumns - 1)
			{
				for(std::size_t i = 0; i < samples.size(); i++)
				{
					withoutReplacementMergingWithVarianceSample& currentSample = samples[i];

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
			}
			row = 0;
		}
		args.estimate = 0;
		for(std::size_t i = 0; i < samples.size(); i++)
		{
			args.estimate += samples[i].weight;
		}
		{
			//Initialize the graph
			for(int i = 0; i < (int)graphVertices.back().size(); i++)
			{
				varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertices.back()[i]);
				if(vertexInfo.indexWithinSelected != -1) vertexInfo.accumulatedMean = 1;
			}
			std::vector<boost::numeric::ublas::matrix<sampling::mpfr_class> > allCovariances(nRows * nColumns);
			allCovariances[nRows * nColumns - 1].resize(graphVertices[nRows * nColumns - 1].size(), graphVertices[nRows * nColumns - 1].size(), false);
			
			for(int time = nRows * nColumns - 2; time >= 0; time--)
			{
				const std::vector<sampling::mpfr_class>& currentInclusionProbabilities = allInclusionProbabilities[time+1];
				const boost::numeric::ublas::matrix<sampling::mpfr_class>& currentSecondOrderInclusionProbabilities = allSecondOrderInclusionProbabilities[time+1];

				allCovariances[time].resize(graphVertices[time].size(), graphVertices[time].size(), false);
				boost::numeric::ublas::matrix<sampling::mpfr_class>& currentCovariance = allCovariances[time];
				const boost::numeric::ublas::matrix<sampling::mpfr_class>& previousCovariance = allCovariances[time+1];

				//First estimate the accumulated means
				for(int particleCounter = 0; particleCounter < (int)graphVertices[time].size(); particleCounter++)
				{
					int currentVertex = graphVertices[time][particleCounter];
					varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, currentVertex);
					vertexInfo.accumulatedMean = 0;
					varianceGraph::out_edge_iterator current, end;
					boost::tie(current, end) = boost::out_edges(currentVertex, varianceEstimationGraph);
					for(; current != end; current++)
					{
						int targetVertex = boost::target(*current, varianceEstimationGraph);
						varianceGraphVertex& targetVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex);
						vertexInfo.accumulatedMean += (targetVertexInfo.accumulatedMean) / currentInclusionProbabilities[targetVertexInfo.indexWithinDesign];
					}
				}
				//now actually do the covariance estimation. 
				for(int particleCounter1 = 0; particleCounter1 < (int)graphVertices[time].size(); particleCounter1++)
				{
					int graphVertex1 = graphVertices[time][particleCounter1];
					varianceGraphVertex& graphVertex1Info = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertex1);
					for(int particleCounter2 = 0; particleCounter2 < (int)graphVertices[time].size(); particleCounter2++)
					{
						int graphVertex2 = graphVertices[time][particleCounter2];
						varianceGraphVertex& graphVertex2Info = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertex2);
						sampling::mpfr_class& currentCovarianceValue = currentCovariance(particleCounter1, particleCounter2);
					
						varianceGraph::out_edge_iterator current1, end1, current2, end2;
						boost::tie(current1, end1) = boost::out_edges(graphVertex1, varianceEstimationGraph);
						for(; current1 != end1; current1++)
						{
							int targetVertex1 = boost::target(*current1, varianceEstimationGraph);
							varianceGraphVertex& targetVertexInfo1 = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex1);
							if(targetVertexInfo1.indexWithinSelected != -1)
							{
								boost::tie(current2, end2) = boost::out_edges(graphVertex2, varianceEstimationGraph);
								for(; current2 != end2; current2++)
								{
									int targetVertex2 = boost::target(*current2, varianceEstimationGraph);
									varianceGraphVertex& targetVertexInfo2 = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex2);
									if(targetVertexInfo2.indexWithinSelected != -1)
									{
										sampling::mpfr_class inclusionProduct = currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign] * currentInclusionProbabilities[targetVertexInfo2.indexWithinDesign];
										if(targetVertex1 == targetVertex2)
										{
											currentCovarianceValue += ((currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign] - inclusionProduct) * targetVertexInfo1.accumulatedMean * targetVertexInfo2.accumulatedMean * graphVertex1Info.trueDensity * graphVertex2Info.trueDensity / (currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign]*inclusionProduct)) + (previousCovariance(targetVertexInfo1.indexWithinSelected, targetVertexInfo2.indexWithinSelected) / (currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign])) * (graphVertex1Info.trueDensity / targetVertexInfo1.trueDensity) * (graphVertex2Info.trueDensity / targetVertexInfo2.trueDensity);
										}
										else
										{
											currentCovarianceValue += ((currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign) - inclusionProduct) * targetVertexInfo1.accumulatedMean * targetVertexInfo2.accumulatedMean * graphVertex1Info.trueDensity * graphVertex2Info.trueDensity / (currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign) * inclusionProduct)) + (previousCovariance(targetVertexInfo1.indexWithinSelected, targetVertexInfo2.indexWithinSelected) / (currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign))) * (graphVertex1Info.trueDensity / targetVertexInfo1.trueDensity) * (graphVertex2Info.trueDensity / targetVertexInfo2.trueDensity);
										}
									}
								}
							}
						}
						if(particleCounter1 == particleCounter2)
						{
							graphVertex1Info.V = currentCovarianceValue;
						}
					}
				}
				//If a variance is zero, the covariances *have* to be zero. 
				for(int particleCounter1 = 0; particleCounter1 < (int)graphVertices[time].size(); particleCounter1++)
				{
					if(currentCovariance(particleCounter1, particleCounter1) == 0)
					{
						for(int particleCounter2 = 0; particleCounter2 < (int)graphVertices[time].size(); particleCounter2++)
						{
							currentCovariance(particleCounter1, particleCounter2) = currentCovariance(particleCounter2, particleCounter1) = 0;
						}
					}
				}
			}
			mpfr_class totalFromGraph = 0;
			{
				varianceGraph::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(0, varianceEstimationGraph);
				for(; current != end; current++)
				{
					int targetVertex = boost::target(*current, varianceEstimationGraph);
					totalFromGraph += boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex).accumulatedMean;
				}
			}
			sampling::mpfr_class totalCovarianceFromGraph = allCovariances[0](0, 0) + allCovariances[0](1, 0) + allCovariances[0](0, 1) + allCovariances[0](1, 1);
			args.varianceEstimate = totalCovarianceFromGraph;
		}
		{
			std::ofstream file("graph.dot");
			vertexPropertyWriter vp(varianceEstimationGraph);
			boost::write_graphviz(file, varianceEstimationGraph, vp);
		}
	}
}
