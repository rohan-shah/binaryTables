conditionalPoisson <- function(rowSums, columnSums, n, seed, keepTables = FALSE)
{
	start <- Sys.time()
	result <- .Call("conditionalPoisson", rowSums, columnSums, n, seed, keepTables, PACKAGE="binaryTables")
	end <- Sys.time()

	estimate <- mpfr(result$estimate)
	varianceEstimate <- mpfr(result$varianceEstimate)
	tables <- list()
	tableWeights <- mpfr("0")[0]
	if(keepTables)
	{
		tables <- result$tables
		tableWeights <- mpfr(result$tableWeights)
	}
	return(new("conditionalPoissonResult", start = start, end = end, estimate = estimate, varianceEstimate = varianceEstimate, n = as.integer(n), seed = as.integer(seed), rowSums = as.integer(rowSums), columnSums = as.integer(columnSums), call = match.call(), tables = tables, keepTables = keepTables, tableWeights = tableWeights))

}
