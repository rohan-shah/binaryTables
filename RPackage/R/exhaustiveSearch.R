exhaustiveSearch <- function(rowSums, columnSums)
{
	.Call("exhaustiveSearch", rowSums, columnSums, PACKAGE="binaryTables")
}
