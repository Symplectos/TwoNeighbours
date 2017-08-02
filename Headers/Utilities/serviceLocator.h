#pragma once

/*************************************************************************************
* Author: Gilles Bellot
* Date:  13/09/2016 - Lenningen - Luxembourg
*
* Desc: locator service
*
* History: 	- 07/10/16: added file logging service
* 			- 13/10/16: program options service added
****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes
#include "../../Headers/Utilities/log.h"
#include "../../Headers/Utilities/programOptions.h"


namespace util
{
// CLASSES //////////////////////////////////////////////////////////////////////////////
/** \addtogroup Utility
 *  @{
 */
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
	// file logging
	static Logger<FileLogPolicy>* getFileLogger() { return fileLogger.get(); }									//!< Returns the file logger.
	static void provideFileLoggingService(const std::shared_ptr<Logger<FileLogPolicy> > providedFileLogger);	//!< Registers the file logging service.

	// program options
	static ProgramOptions* getProgOpts() { return progOpts.get(); }												//!< Gets the program options.
	static void provideProgramOptions(const std::shared_ptr<ProgramOptions> providedProgramOptions);			//!< Registers the program options as a service.
};
/** @}*/
}
