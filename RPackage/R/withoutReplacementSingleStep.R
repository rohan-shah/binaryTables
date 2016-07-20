withoutReplacementSingleStep <- function(rowSums, columnSums, n, seed)
{
	start <- Sys.time()
	result <- .Call("withoutReplacementSingleStep", rowSums, columnSums, n, seed, PACKAGE="binaryTables")
	end <- Sys.time()
	return(new("withoutReplacementSingleStepResult", start = start, end = end, estimate = mpfr(result), n = as.integer(n), seed = as.integer(seed), rowSums = as.integer(rowSums), columnSums = as.integer(columnSums), call = match.call()))
}
