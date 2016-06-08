conditionalPoisson <- function(rowSums, columnSums, n, seed)
{
	result <- .Call("conditionalPoisson", rowSums, columnSums, n, seed, PACKAGE="binaryTables")
	return(mpfr(result))
}
