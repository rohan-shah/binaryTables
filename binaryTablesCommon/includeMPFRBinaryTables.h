#ifndef INCLUDE_MPFR_BINARY_TABLES_HEADER_GUARD
#define INCLUDE_MPFR_BINARY_TABLES_HEADER_GUARD
//#include "mpreal.h"
//typedef mpfr::mpreal mpfr_class;
#include <boost/multiprecision/mpfr.hpp>
namespace binaryTables
{
	typedef boost::multiprecision::mpfr_float mpfr_class;
	typedef boost::multiprecision::mpz_int mpz_class;
}
#endif 
