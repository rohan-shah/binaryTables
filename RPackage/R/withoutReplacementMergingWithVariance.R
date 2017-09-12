withoutReplacementMergingWithVariance <- function(rowSums, columnSums, n, seed)
{
	start <- Sys.time()
	result <- .Call("withoutReplacementMergingWithVariance", rowSums, columnSums, n, seed, PACKAGE="binaryTables")
	end <- Sys.time()
	estimate <- mpfr(result$estimate)
	varianceEstimate <- mpfr(result$varianceEstimate)
	return(new("withoutReplacementMergingWithVarianceResult", start = start, end = end, estimate = estimate, n = as.integer(n), seed = as.integer(seed), rowSums = as.integer(rowSums), columnSums = as.integer(columnSums), call = match.call(), variance = varianceEstimate))
}
