/*
 * File:   hn.cpp
 * Author: Gilles Bellot
 * Date:   23/07/2012 - Dortmund - Germany
 *
 * Description: TwoNeighbours is an attempt to rewrite TN and HN in modern C++
 *
 * History:	- 24/07/2012: program options added
 *			- 25/07/2012: mpfr precision can now be set with --precision arg
 *			- 22/01/2013: new timer added
 *			- 22/03/2013: mpfr is now thread safe
 *			- 25/04/2014: added crude console menu functionality
 *			- 02/11/2015: tabula rasa
 *			- 02/11/2015: added openMP support
 *			- 06/10/2016: tabula rasa
 *			- 07/10/2016: added file logger service
 *			- 13/10/2016: program options added
 *			- 01/08/2017: tabula rasa
 *			- 01/08/2017: Expected class for better error handling was added (expected.h)
 *			- 01/08/2017: added a service locator (seriveLocator.h)
 *			- 01/08/2017: file logger service was added (log.h)
 *			- 01/08/2017: program options were added (programOptions.h)
 *			- 01/08/2017: support for boost matrices was added (matrix.h)
 *			- 01/08/2017: linear algebra class added (linearAlgebra.h), most important algorithms: LLL and short vectors
 *			- 01/08/2017: general lattice class was added (lattice.h)
 *
 *
 * ToDo: everything :(
 */

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// multi-precision numbers
#include <gmp.h>											// the GNU MP libraries
#include <gmpxx.h>
#include <mpfr.h>											// the GNU MPFR library
#include "Headers/Multi-Precision/mpreal.h"					// GNU MPFR C++ wrapper by Pavel Holoborodko

// boost ublas
#include <boost/numeric/ublas/io.hpp>						// input and output for ublas objects
#include <boost/numeric/ublas/symmetric.hpp>				// symmetrical matrices

// CUDA includes
#include <cuda_runtime.h>									// cuda runtim library

// bell0bytes includes

// util
#include "Headers/Utilities/serviceLocator.h"				// service locator
#include "Headers/Utilities/programOptions.h"				// program options
#include "Headers/Utilities/expected.h"						// expected error handling

// mathematics
//#include "Headers/Mathematics/matrix.h"						// matrix helper
#include "Headers/Mathematics/linearAlgebra.h"				// linear algebra algorithms
#include "Headers/Mathematics/lattice.h"					// definition of lattices

// DEFINITIONS //////////////////////////////////////////////////////////////////////////

// classes
class TwoNeighbours
{
private:
	bool hasStarted;						// true iff the application was successfuly started
	bool hasFileLogger;						// true iff a file logger is available

	void printStartingLog();				// prints a starting log


public:
	// constructor
	TwoNeighbours();
	~TwoNeighbours();

	util::Expected<void> init(int argc, char** argv);		// initializes the application; creates and registers services
	void shutdown(util::Expected<void>* result = NULL);		// release memory and shut the application down
};

// FUNCTIONS ////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Main Function ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	TwoNeighbours tn;
	util::Expected<void> initialization = tn.init(argc, argv);
	if(initialization.wasSuccessful())
	{
		// run HN

		// read gram matrix of starting lattice
		boost::numeric::ublas::symmetric_matrix<mpz_class> A(1,1);
		std::cin >> A;

		// create lattice from gramian
		mathematics::Lattice L(&A);

		// print starting lattice
		L.print(util::ServiceLocator::getProgOpts()->verboseLevel);

		// shut down
		tn.shutdown();

		// return success
		return 0;
	}
	else
	{
		// shut down with error message
		tn.shutdown(&initialization);
		return -1;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Initialization ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
TwoNeighbours::TwoNeighbours() : hasStarted(false), hasFileLogger(false)
{

}

util::Expected<void> TwoNeighbours::init(int argc, char** argv)
{
	// register file logger service
	try
	{
		std::shared_ptr<util::Logger<util::FileLogPolicy> > fileLogger(new util::Logger<util::FileLogPolicy>("../Logs/TwoNeighbours.log"));
		util::ServiceLocator::provideFileLoggingService(fileLogger);
		util::ServiceLocator::getFileLogger()->setThreadName("mainThread");
	}
	catch(std::runtime_error&)
	{
		return std::runtime_error("Critical error: Unable to create the file logger!");
	}
	hasFileLogger=true;

	try
	{
		// create program options and register them as a service
		std::shared_ptr<util::ProgramOptions> progOpts(new util::ProgramOptions(argc, argv));
		util::ServiceLocator::provideProgramOptions(progOpts);
	}
	catch(std::runtime_error&)
	{
		return std::runtime_error("Critical error: Unable to parse the program options!");
	}

	// write starting activity to log file
	printStartingLog();

	// return success
	hasStarted = true;
	return {};
}

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Starting Log /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// write start log
void TwoNeighbours::printStartingLog()
{
	// get CPU info and write starting message
	std::stringstream cpuStream;
	cpuStream << "TwoNeighbours has started running on:\n";
	std::ifstream cpuInfo("/proc/cpuinfo");
	std::string readFile;

	// get number of cores
	int ncore = 0;
	while(getline(cpuInfo, readFile))
	{
		if(readFile.find("processor") != std::string::npos)
			ncore++;
	}

	// get back to the beginning of the file and read in the processor data
	cpuInfo.clear();
	cpuInfo.seekg(0, std::ios::beg);
	readFile.clear();
	while(getline(cpuInfo, readFile))
	{
		if(readFile.find("name") != std::string::npos)
		{
			cpuStream << "\n\t\t\t\t\t" << readFile.substr(13, std::string::npos-10);
			cpuStream << "\n\t\t\t\t\t\t Number of cores: " << ncore;
		}
		if(readFile.find("MHz") != std::string::npos)
			cpuStream << "\n\t\t\t\t\t\t CPU Clock: " << readFile.substr(11, std::string::npos-3) << " Mhz";
		if(readFile.find("cache size") != std::string::npos)
			cpuStream << "\n\t\t\t\t\t\t Cache size:" << readFile.substr(12, std::string::npos-3);
		if(readFile.find("management") != std::string::npos)
			break;
	}
	cpuStream << "\n";
	util::ServiceLocator::getFileLogger()->print<util::SeverityType::info>(std::stringstream(cpuStream.str()));

	// write status of loaded libraries
	std::stringstream streamlib;
	streamlib << "The following libraries were loaded successfully:\n";
	streamlib << "\n\t\t\t\t\t\t* GNU MP " << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL;
	streamlib << "\n\t\t\t\t\t\t* MPFR " << MPFR_VERSION_STRING << " (prec: " << mpfr::bits2digits(mpfr::mpreal::get_default_prec()) << ")";
	streamlib << "\n\t\t\t\t\t\t* Boost uBlas " << BOOST_LIB_VERSION;
	if(util::ServiceLocator::getProgOpts()->useCUDA)
		streamlib << "\n\t\t\t\t\t\t* CUDA " << CUDART_VERSION << " with cuBlas.\n";
	else
		streamlib << "\n";
	util::ServiceLocator::getFileLogger()->print<util::SeverityType::info>(std::stringstream(streamlib.str()));

	if(util::ServiceLocator::getProgOpts()->useCUDA)
	{
		// write gpu status
		std::stringstream gpuStream;

		const int kb = 1024;
		const int mb = kb * kb;
		gpuStream << "GPU Information:\n";

		int nDevices;
		cudaGetDeviceCount(&nDevices);
		for (int i = 0; i < nDevices; i++)
		{
			cudaDeviceProp props;
			cudaGetDeviceProperties(&props, i);
			gpuStream << "\n\t\t\t\t\t " << i << ": " << props.name << ": CUDA " << props.major << "." << props.minor;
			gpuStream << "\n\t\t\t\t\t\tGlobal memory:   " << props.totalGlobalMem / mb << "mb";
			gpuStream << "\n\t\t\t\t\t\tShared memory:   " << props.sharedMemPerBlock / kb << "kb";
			gpuStream << "\n\t\t\t\t\t\tConstant memory: " << props.totalConstMem / kb << "kb";
			gpuStream << "\n\t\t\t\t\t\tBlock registers: " << props.regsPerBlock;

			gpuStream << "\n\n\t\t\t\t\t\tWarp size           :   " << props.warpSize;
			gpuStream << "\n\t\t\t\t\t\tThreads per block   :   " << props.maxThreadsPerBlock;
			gpuStream << "\n\t\t\t\t\t\tMax block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]";
			gpuStream << "\n\t\t\t\t\t\tMax grid dimensions : [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]\n";
		}
		util::ServiceLocator::getFileLogger()->print<util::SeverityType::info>(std::stringstream(gpuStream.str()));
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Shutdown /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
TwoNeighbours::~TwoNeighbours()
{}

void TwoNeighbours::shutdown(util::Expected<void>* result)
{
	// check for error message
	if(result != NULL && !result->isValid())
	{
		// the application encountered an error -> try to clean up and to log the error message
		try
		{
			// clean up

			// throw the actual error message
			result->get();
		}
		catch(std::runtime_error& e)
		{
			if(hasFileLogger)
			{
				std::stringstream errorMessage;
				errorMessage << "TwoNeighbours was shut down: " << e.what();
				util::ServiceLocator::getFileLogger()->print<util::SeverityType::error>(std::stringstream(errorMessage.str()));
			}
			else
				std::cerr << "TwoNeighbours is shutting down: " << e.what() << "\n";
			return;
		}
	}

	// no error -> clean up and shut down normally
	util::ServiceLocator::getFileLogger()->print<util::SeverityType::info>("TwoNeighbours was successfully shut down!");
}

