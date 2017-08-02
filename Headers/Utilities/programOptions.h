#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	programOptions.h
 * @date	3/11/2015 - Dortmund - Germany
 *
 * @brief	Program options created from command line parameters with the help of BOOST.
 *
 * @section History
 *
 *  - 13/03/17: changed verbose to int
 *  - 21/01/17: added option to enable verbose mode
 *  - 16/01/17: fixed a small oversight with precision settings
 *  - 13/10/16: added CUDA option
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes util
#include "expected.h"

namespace util
{
// DEFINITIONS //////////////////////////////////////////////////////////////////////////
/*! @brief The ProgramOptions define the behaviour of TwoNeighbours.
 *
 * @section Description
 *
 * The program options are created from the command line parameters by the BOOST program options library.
 *
 * The following options are currently available:
 *  - verboseLevel: an int to control the verbose level of the entire application
 *  - debugMode: a boolean; if true, the application runs in debug mode, ensuring the correctness of all the algorithms. Obviously this is extremely slow.
 *  - useCUDA: a boolean which is true if and only if CUDA algorithms are to be prioritized
 *  - outputPrec: an int to set the output precision
 *  - mpfrPrec: an int to set the precision of mpfr related calculations
 *
 */
class ProgramOptions
{
private:
	Expected<void> createCommandLineParameters(int argc, char** argv);	// called from the default constructor, this function creates the program options from the passed command line parameters

public:
	ProgramOptions(int argc, char** argv);	// the constructor sets all members to their default values and then calls the createCommandLineParameters function
	~ProgramOptions();						// the destructor does what destructors do, it destroys!

	// general options
	int verboseLevel;		//!< This flag defines the verbose mode. The higher the number, the more information will be shown on the default output stream. Default: 0 (off).
	bool debugMode;			//!< True if and only if the debug mode is activated. Default: false.

	//!< The debug mode is obviously slower, shows debug information and computes internal tests to ensure correctness of algorithms.

	// cuda
	bool useCUDA;			//!< True if and only if CUDA algorithms should be prioritized.

	// output
	int outputPrec;			//!< Defines the output precision of BOOST.
	int mpfrPrec;			//!< Defines the precision to use for calculations involving mpfr.
};

}
