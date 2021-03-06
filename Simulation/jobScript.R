source("./generateScenarios.R")
SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
cat("SCENARIO_INDEX=", SCENARIO_INDEX, "\n", sep="")

library(binaryTables)
outputFile <- file.path("results", scenarios[SCENARIO_INDEX, "file"])
tmpFile <- paste0(outputFile, ".tmp")
scenario <- scenarios[SCENARIO_INDEX, "scenario"]
replications <- scenarios[SCENARIO_INDEX, "replications"]
sampleSize <- scenarios[SCENARIO_INDEX, "sampleSize"]
method <- scenarios[SCENARIO_INDEX, "method"]

if(scenario == "exact12x12")
{
	columnSums <- rowSums <- rep(2, 12)
} else if(scenario == "50x50_1")
{
	rowSums <- c(10, 8, 11, 11, 13, 11, 10, 9, 7, 9, 10, 16, 11, 9, 12, 14, 12, 7, 9, 10, 10, 6, 11, 8, 9, 8, 14, 12, 5, 10, 10, 8, 7, 8, 10, 10, 14, 6, 10, 7, 13, 4, 6, 8, 9, 15, 11, 12, 10, 6)
	columnSums <- c(9, 6, 12, 11, 9, 8, 8, 11, 9, 11, 13, 7, 10, 8, 9, 7, 8, 3, 10, 11, 13, 7, 5, 11, 10, 9, 10, 13, 9, 9, 7, 7, 6, 8, 10, 12, 8, 12, 16, 12, 15, 12, 13, 13, 10, 7, 12, 13, 6, 11)
} else if(scenario == "pathological_100")
{
	m <- 100
	beta <- 0.6
	gamma <- 0.7
	n <- m + floor(beta*m) - floor(gamma*m)
	columnSums <- c(rep(1, n), floor(gamma*m))
	rowSums <- c(rep(1, m), floor(beta*m))
} else if(scenario == "pathological_200")
{
	m <- 200
	beta <- 0.6
	gamma <- 0.7
	n <- m + floor(beta*m) - floor(gamma*m)
	columnSums <- c(rep(1, n), floor(gamma*m))
	rowSums <- c(rep(1, m), floor(beta*m))
} else if(scenario == "pathological_300")
{
	m <- 300
	beta <- 0.6
	gamma <- 0.7
	n <- m + floor(beta*m) - floor(gamma*m)
	columnSums <- c(rep(1, n), floor(gamma*m))
	rowSums <- c(rep(1, m), floor(beta*m))
} else
{
	stop("Unknown scenario")
}
counter <- 1
if(file.exists(outputFile))
{
	load(outputFile)
	counter <- length(results) + 1
} else results <- list()
if(method == "CP")
{
	while(counter < replications + 1)
	{
		results[[counter]] <- conditionalPoisson(columnSums = columnSums, rowSums = rowSums, seed = counter + 100000L*SCENARIO_INDEX, n = sampleSize)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method == "WOR")
{
	while(counter < replications + 1)
	{
		results[[counter]] <- withoutReplacement(columnSums = columnSums, rowSums = rowSums, seed = counter + 100000L*SCENARIO_INDEX, n = sampleSize)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}
} else if(method %in% c("WOR-Merge1", "WOR-Merge50", "WOR-Merge12", "WOR-Merge100", "WOR-Merge200", "WOR-Merge300"))
{
	if(method == "WOR-Merge1") mergeFrequency <- 1
	else if(method == "WOR-Merge50") mergeFrequency <- 50
	else if(method == "WOR-Merge12") mergeFrequency <- 12
	else if(method == "WOR-Merge100") mergeFrequency <- 100
	else if(method == "WOR-Merge200") mergeFrequency <- 200
	else if(method == "WOR-Merge300") mergeFrequency <- 300
	else stop("Internal error")
	while(counter < replications + 1)
	{
		results[[counter]] <- withoutReplacementMerging(columnSums = columnSums, rowSums = rowSums, seed = counter + 100000L*SCENARIO_INDEX, n = sampleSize, mergeFrequency = mergeFrequency)
		save(results, file = tmpFile)
		file.rename(from = tmpFile, to = outputFile)
		counter <- counter + 1
	}

} else
{
	stop("Unknown method")
}
