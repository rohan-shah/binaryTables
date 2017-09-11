library(xtable)
source("./generateScenarios.R")
load("./summarised.RData")
scenarioLabels <- c("pathological_100", "pathological_200", "pathological_300", "50x50_1", "50x50_2")
scenarioCaptions <- list()
scenarioCaptions[["pathological_100"]] <- "Simulation results for counting the number of binary tables with $\\beta = 0.6, \\gamma = 0.7$ and $m = 100$"
scenarioCaptions[["pathological_200"]] <- "Simulation results for counting the number of binary tables with $\\beta = 0.6, \\gamma = 0.7$ and $m = 200$"
scenarioCaptions[["pathological_300"]] <- "Simulation results for counting the number of binary tables with $\\beta = 0.6, \\gamma = 0.7$ and $m = 300$"
allData <- cbind(scenarios[,c("method", "scenario", "sampleSize")], averageEstimates, wnrv, relativeErrors, stringsAsFactors=FALSE)
allData$averageEstimates <- as.numeric(allData$averageEstimates)
for(scenarioLabel in scenarioLabels)
{
	subsetted <- allData[allData[,"scenario"] == scenarioLabel,]
	subsetted <- subsetted[,c("method", "sampleSize", "averageEstimates", "wnrv", "relativeErrors")]
	colnames(subsetted)[1:5] <- c("Method", "$\\sampleSize$", "Average estimate", "WNRV", "RE")
	if(scenarioLabel %in% names(scenarioCaptions))
	{
		table <- xtable(subsetted, display = c("s", "s", "d", "e", "e", "e"), caption = scenarioCaptions[[scenarioLabel]])
	}
	else table <- xtable(subsetted, display = c("s", "s", "d", "e", "e", "e"))
	print(table, math.style.exponents = TRUE, include.rownames=FALSE, sanitize.colnames.function = identity)
}
