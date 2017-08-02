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
/** \addtogroup Utility
 *  @{
 */
// DEFINITIONS //////////////////////////////////////////////////////////////////////////
/*! \brief The ProgramOptions define the behaviour of TwoNeighbours.
 *
 */
class ProgramOptions
{
private:
	Expected<void> createCommandLineParameters(int argc, char** argv);	//!< Called from the default constructor, this function creates the program options from the passed command line parameters.

public:
	ProgramOptions(int argc, char** argv);	//!< The constructor sets all members to their default values and then calls the createCommandLineParameters function.
	~ProgramOptions();						//!< The destructor does what destructors do, it destroys!

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
/** @}*/
}
