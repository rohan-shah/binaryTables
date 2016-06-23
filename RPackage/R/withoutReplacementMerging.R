withoutReplacementMerging <- function(rowSums, columnSums, n, seed, mergeFrequency)
{
	start <- Sys.time()
	result <- .Call("withoutReplacementMerging", rowSums, columnSums, n, seed, mergeFrequency, PACKAGE="binaryTables")
	end <- Sys.time()
	estimate <- mpfr(result)
	return(new("withoutReplacementMergingResult", mergeFrequency = as.integer(mergeFrequency), start = start, end = end, estimate = estimate, n = as.integer(n), seed = as.integer(seed), rowSums = as.integer(rowSums), columnSums = as.integer(columnSums), call = match.call()))
}
