#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
#include "includeMPFRBinaryTables.h"
#include "crudeMC.h"
#include "exhaustiveSearch.h"
extern "C" const char* package_name = "binaryTables";
R_CallMethodDef callMethods[] = 
{
	{"crudeMC", (DL_FUNC)&binaryTables::crudeMC, 4},
	{"exhaustiveSearch", (DL_FUNC)&binaryTables::exhaustiveSearch, 2}, 
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
