.onLoad <- function(libname, pkgname)
{
	library.dynam(package="binaryTables", chname="binaryTables", lib.loc = .libPaths())
}
