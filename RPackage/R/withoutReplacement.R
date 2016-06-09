withoutReplacement <- function(rowSums, columnSums, n, seed)
{
	rowSums <- rowSums[rowSums > 0]
	columnSums <- columnSums[columnSums > 0]
	result <- .Call("withoutReplacement", rowSums, columnSums, n, seed, PACKAGE="binaryTables")
	return(mpfr(result))
}
