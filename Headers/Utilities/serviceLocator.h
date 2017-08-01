#pragma once

/***************************************************************************************
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

class ServiceLocator
{
private:
	static std::shared_ptr<Logger<FileLogPolicy> > fileLogger;	// the file logger
	static std::shared_ptr<ProgramOptions> progOpts;			// program options

public:
	// file logging
	static Logger<FileLogPolicy>* getFileLogger() { return fileLogger.get(); };
	static void provideFileLoggingService(const std::shared_ptr<Logger<FileLogPolicy> > providedFileLogger);

	// program options
	static ProgramOptions* getProgOpts() { return progOpts.get(); };
	static void provideProgramOptions(const std::shared_ptr<ProgramOptions> providedProgramOptions);
};

}
