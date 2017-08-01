/*
 * File:   programOptions.h
 * Author: Gilles Bellot
 * Date:   03/11/2015 - Dortmund - Germany
 *
 * Description: program options, uses boost
 *
 * History: 	- 13/10/16: added CUDA option
 * 				- 16/01/17: fixed a small oversight with precision settings
 * 				- 21/01/17: added option to enable verbose mode
 * 				- 13/03/17: changed verbose to int
 *
 * ToDo: everything :(
 */

#pragma once

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes util
#include "expected.h"

namespace util
{
// DEFINITIONS //////////////////////////////////////////////////////////////////////////
class ProgramOptions
{
private:
	Expected<void> createCommandLineParameters(int argc, char** argv);

public:
	ProgramOptions(int argc, char** argv);
	~ProgramOptions();

	// general options
	int verboseLevel;		// verbose mode - slower, shows additional information on the standard output stream
	bool debugMode;			// debug mode - slower, shows debug information and computes internal tests to ensure correctness of algorithms

	// cuda
	bool useCUDA;			// true iff we want to use CUDA algorithms

	// output
	int outputPrec;			// boost output precision
	int mpfrPrec;			// precision to use for mpfr

	// isom autom
	int automAlgo;			// 0: Plesken-Souvignier ; 1: Regev-Oded
	int isomAlgo;			// 0: Plesken-Souvignier ; 1: Regev-Oded
};
}
