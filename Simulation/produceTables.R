library(xtable)
source("./generateScenarios.R")
load("./summarised.RData")
scenarioLabels <- c("exact12x12", "50x50_1", "50x50_2", "4_8", "2_10")
scenarioCaptions <- list()
scenarioCaptions[["exact12x12"]] <- "Simulation results for counting the number of $12 \\times 12$ binary tables with all row and column sums equal to 2"
scenarioCaptions[["4_8"]] <- "Simulation results for counting the number of $12 \\times 12$ binary tables with $\\mathbf r = \\mathbf c = \\(4,4,4,4,\\protect\\underbrace{8, \\dots, 8}_{\\text{8 times}}\\)$"
scenarioCaptions[["2_10"]] <- "Simulation results for counting the number of $12 \\times 12$ binary tables with $\\mathbf r = \\mathbf c = \\(2,2,\\protect\\underbrace{10, \\dots, 10}_{\\text{10 times}}\\)$"
allData <- cbind(scenarios[,c("method", "scenario", "sampleSize")], wnrv)
for(scenarioLabel in scenarioLabels)
{
	subsetted <- allData[allData[,"scenario"] == scenarioLabel,]
	subsetted <- subsetted[,c("method", "sampleSize", "wnrv")]
	colnames(subsetted)[1:3] <- c("Method", "$\\sampleSize$", "WNRV")
	if(scenarioLabel %in% names(scenarioCaptions))
	{
		table <- xtable(subsetted, display = c("s", "s", "d", "e"), caption = scenarioCaptions[[scenarioLabel]])
	}
	else table <- xtable(subsetted, display = c("s", "s", "d", "e"))
	print(table, math.style.exponents = TRUE, include.rownames=FALSE, sanitize.colnames.function = identity)
}
