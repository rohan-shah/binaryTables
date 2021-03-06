cmake_minimum_required(VERSION 3.1)
if(WIN32)
	#It seems adding a non-satsfied dependency like this only works on windows
	add_custom_target(copyText ALL DEPENDS nonExistentFilename)
else()
	add_custom_target(copyText ALL)
endif()

set(RFiles R/onLoad.R R/crudeMC.R R/setDefaultPrec.R R/exhaustiveSearch.R R/conditionalPoisson.R R/withoutReplacement.R R/classes.R R/withoutReplacementMerging.R R/calculateConditionalPoissonDensity.R R/withoutReplacementSingleStep.R R/conditionalPoissonBootstrap.R R/conditionalPoissonBootstrapMerging.R R/withoutReplacementMergingWithVariance.R)
set(RPackageFiles DESCRIPTION NAMESPACE ${RFiles})

#Find location of Rcpp_exceptions.R 
get_property(RCPP_INCLUDE_DIRECTORIES TARGET Rcpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
foreach(RCPP_INCLUDE_DIRECTORY IN LISTS RCPP_INCLUDE_DIRECTORIES)
	if(EXISTS ${RCPP_INCLUDE_DIRECTORY}/Rcpp_exceptions.R)
		set(RCPP_EXCEPTIONS ${RCPP_INCLUDE_DIRECTORY}/Rcpp_exceptions.R)
	endif()
endforeach()

#Copy R package files
if(WIN32)
	foreach(BUILD_TYPE IN LISTS CMAKE_CONFIGURATION_TYPES)
		foreach(PACKAGE_FILE IN LISTS RPackageFiles)
			add_custom_command(TARGET copyText POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_SOURCE_DIR}/RPackage/${PACKAGE_FILE}" "${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/${PACKAGE_FILE}")
		endforeach()
		add_custom_command(TARGET copyText POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy "${RCPP_EXCEPTIONS}" "${CMAKE_BINARY_DIR}/RPackage/${BUILD_TYPE}/R/Rcpp_exceptions.R")
	endforeach()
else()
	foreach(PACKAGE_FILE IN LISTS RPackageFiles)
		add_custom_command(TARGET copyText POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy "${CMAKE_SOURCE_DIR}/RPackage/${PACKAGE_FILE}" "${CMAKE_BINARY_DIR}/RPackage/${PACKAGE_FILE}")
	endforeach()
	add_custom_command(TARGET copyText POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy "${RCPP_EXCEPTIONS}" "${CMAKE_BINARY_DIR}/RPackage//R/Rcpp_exceptions.R")
endif()
