#include "../../Headers/Utilities/serviceLocator.h"

namespace util
{

std::shared_ptr<Logger<FileLogPolicy> > ServiceLocator::fileLogger = NULL;
std::shared_ptr<ProgramOptions> ServiceLocator::progOpts = NULL;

// file logging
void ServiceLocator::provideFileLoggingService(std::shared_ptr<Logger<FileLogPolicy> > providedFileLogger)
{
	fileLogger = providedFileLogger;
}

// program options
void ServiceLocator::provideProgramOptions(std::shared_ptr<ProgramOptions> provProgOpts)
{
	progOpts = provProgOpts;
}

}
