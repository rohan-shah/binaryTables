source("./generateScenarios.R")
library(Rmpfr)
allResults <- list()
for(i in 1:nrow(scenarios))
{
	path <- file.path("results", scenarios[i, "file"])
	if(file.exists(path))
	{
		load(path)
		allResults[[i]] <- results
	}
}
secondsPerRun <- lapply(allResults, function(x) unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
averageSecondsPerRun <- unlist(lapply(secondsPerRun, function(x) if(is.null(x)) NA else mean(x)))
averageEstimatesFunc <- function(x)
{
	if(is.null(x))
	{
		return(mpfr(NA))
	}
	values <- do.call(c, lapply(x, function(y) y@estimate))
	values <- mpfr(values, prec = 3*getPrec(values))
	total <- sum(values)
	return(total/length(values))
}
averageEstimates <- do.call(c, lapply(allResults, averageEstimatesFunc))
varianceFunc <- function(x)
{
	if(is.null(x))
	{
		return(mpfr(NA))
	}
	#There is no var for Rmpfr objects unfortunately. 
	values <- do.call(c, lapply(x, function(y) y@estimate))
	values <- mpfr(values, prec = 3*getPrec(values))
	total <- sum(values)

	totalSquared <- sum(values*values)
	return(totalSquared / length(values) - (total / length(values))^2)
}
variances <- do.call(c, lapply(allResults, varianceFunc))
workNormalizedVariance <- formatMpfr(variances * averageSecondsPerRun, digits = 10)

relativeErrors <- sqrt(variances) / averageEstimates
wnrv <- as.numeric(variances * averageSecondsPerRun / (averageEstimates^2))

variances <- formatMpfr(variances, digits = 10)
relativeErrors <- formatMpfr(relativeErrors, digits = 10)
averageEstimates <- formatMpfr(averageEstimates, digits = 20)

save(allResults, averageEstimates, averageSecondsPerRun, secondsPerRun, variances, workNormalizedVariance, wnrv, relativeErrors, file = "summarised.RData")
