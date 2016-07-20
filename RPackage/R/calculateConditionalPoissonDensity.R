calculateConditionalPoissonDensity <- function(tables)
{
	densities <- .Call("calculateConditionalPoissonDensity", tables, PACKAGE="binaryTables")
	return(formatMpfr(mpfr(densities), 20))
}
