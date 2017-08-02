#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	serviceLocator.h
 * @date	13/09/2016 - Lenningen - Luxembourg
 *
 * @brief	A service locator.
 *
 * @section Description
 *
 * Registers and provides services to the entire program.
 * See <a href="http://gameprogrammingpatterns.com/service-locator.html">Game Programming Patterns</a>, by R. Nystrom for further details.
 *
 * @section History
 *
 *  - 13/10/16: program options service added
 *  - 07/10/16: added file logging service
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes
#include "../../Headers/Utilities/log.h"
#include "../../Headers/Utilities/programOptions.h"


namespace util
{
// CLASSES //////////////////////////////////////////////////////////////////////////////
/*! \brief Registers and provides services to the entire program.
 *
 *  See <a href="http://gameprogrammingpatterns.com/service-locator.html">Game Programming Patterns</a>, by R. Nystrom for further details.
 */
class ServiceLocator
{
private:
	static std::shared_ptr<Logger<FileLogPolicy> > fileLogger;	//!< The file logger.
	static std::shared_ptr<ProgramOptions> progOpts;			//!< The program options.

public:
	/*!
	 *   @brief Returns the file logger.
	 *
	 *   @return Logger<FileLogPolicy>* Returns a pointer to the actual file logger, which can then be used to write messages to a log file.
	*/
	static Logger<FileLogPolicy>* getFileLogger() { return fileLogger.get(); }									// returns the file logger
	static void provideFileLoggingService(const std::shared_ptr<Logger<FileLogPolicy> > providedFileLogger);	// registers the file logging service

	/*!
	 *   @brief Returns the program options.
	 *
	 *   @return ProgramOptions* Returns a pointer to the actual program options.
	*/
	static ProgramOptions* getProgOpts() { return progOpts.get(); }												// gets the program options
	static void provideProgramOptions(const std::shared_ptr<ProgramOptions> providedProgramOptions);			// registers the program options as a service
};

}
