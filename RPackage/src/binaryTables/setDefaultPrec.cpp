#include "setDefaultPrec.h"
#include "includeMPFRBinaryTables.h"
namespace binaryTables
{
	SEXP setDefaultPrec(SEXP prec_sexp)
	{
	BEGIN_RCPP
		binaryTables::mpfr_class::default_precision(Rcpp::as<int>(prec_sexp));
	VOID_END_RCPP
		return R_NilValue;
	}
}
