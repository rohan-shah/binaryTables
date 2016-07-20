#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
#include "includeMPFRBinaryTables.h"
#include "crudeMC.h"
#include "exhaustiveSearch.h"
#include "conditionalPoissonBinaryTables.h"
#include "setDefaultPrec.h"
#include "withoutReplacement.h"
#include "withoutReplacementMerging.h"
#include "calculateConditionalPoissonDensity.h"
#include "withoutReplacementSingleStep.h"
extern "C" const char* package_name = "binaryTables";
R_CallMethodDef callMethods[] = 
{
	{"crudeMC", (DL_FUNC)&binaryTables::crudeMC, 4},
	{"exhaustiveSearch", (DL_FUNC)&binaryTables::exhaustiveSearch, 2}, 
	{"conditionalPoisson", (DL_FUNC)&binaryTables::conditionalPoisson, 5},
	{"withoutReplacement", (DL_FUNC)&binaryTables::withoutReplacement, 5},
	{"withoutReplacementMerging", (DL_FUNC)&binaryTables::withoutReplacementMerging, 5},
	{"setDefaultPrec", (DL_FUNC)&binaryTables::setDefaultPrec, 1},
	{"calculateConditionalPoissonDensity", (DL_FUNC)&binaryTables::calculateConditionalPoissonDensity, 1},
	{"withoutReplacementSingleStep", (DL_FUNC)&binaryTables::withoutReplacementSingleStep, 4},
	{NULL, NULL, 0}
};
RcppExport void R_init_binaryTables(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* packageCallMethods = callMethods;
	while(packageCallMethods->name != NULL) packageCallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, packageCallMethods);

	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
	init_Rcpp_cache();

	//Setup default precision
	binaryTables::mpfr_class::default_precision(1024);
}
