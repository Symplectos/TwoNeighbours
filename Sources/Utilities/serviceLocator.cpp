/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	serviceLocator.cpp
 * @date	13/09/2016 - Lenningen - Luxembourg
 *
 * @brief	Implements the service locator defined in serviceLocator.h
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

#include "../../Headers/Utilities/serviceLocator.h"

namespace util
{

std::shared_ptr<Logger<FileLogPolicy> > ServiceLocator::fileLogger = NULL;
std::shared_ptr<ProgramOptions> ServiceLocator::progOpts = NULL;

/*!
 *   @brief Provides a file logging service to the entire program.
 *
 *   @param providedFileLogger A shared pointer to an instance of a FileLogger class.
 *   @return void
 */
void ServiceLocator::provideFileLoggingService(std::shared_ptr<Logger<FileLogPolicy> > providedFileLogger)
{
	fileLogger = providedFileLogger;
}

/*!
 *   @brief Provides program options which were created from command line parameters as a service to the entire program.
 *
 *   @param providedProgramOptions A shared pointer to an instance of ProgramOptions.
 *   @return void
 */
void ServiceLocator::provideProgramOptions(std::shared_ptr<ProgramOptions> providedProgramOptions)
{
	progOpts = providedProgramOptions;
}

}
