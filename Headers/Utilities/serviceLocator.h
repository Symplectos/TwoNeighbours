#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	serviceLocator.h
 * @date	13/09/2016 - Lenningen - Luxembourg
 *
 * @brief	A service locator.
 *
 * ### Description
 *
 * A service locator is a decoupling pattern to provide a global point of access to a service without coupling to the concrete class that implements it.<br>
 * See <a href="http://gameprogrammingpatterns.com/service-locator.html">Game Programming Patterns</a>, by R. Nystrom for further details.
 *
 * ### History
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
 *  At the moment, and this is clearly work in progress, when the TwoNeighbours starts, the service locator is created and the program options, defined in programOptions.h,
 *  are read from the console parameters and registered as a service. In addition a file logger, defined in log.h, is registered as a service.
 *  See <a href="http://gameprogrammingpatterns.com/service-locator.html">Game Programming Patterns</a>, by R. Nystrom for further details.
 *
 *  ## Example
 *  Here is an example of how to use the service locator class:
 *  ServiceLocator::getFileLogger()->print("Critical error! Not enough papers published!");
 *

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
