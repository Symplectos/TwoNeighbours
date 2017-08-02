/*!
 *
 *  @mainpage
 *  TwoNeighbours is an attempt to rewrite the two neighbours program (TN) by B. Hemkemeier and the hermitian neighbours program (HN) by A. Schiemann in modern C++.
 *  @author		Gilles Bellot
 *  @version	0.0.0.1
 *  @pre		GMP, MPFR, BOOST and CUDA must be installed.
 *  @bug		No known bugs.
 *  @date		01/08/2017 - Dortmund - Germany
 *  @section 	History
 *
 * - 02/08/2017: changed comment structure to allow for automatic documentation generation by Doxygen
 * - 01/08/2017: general lattice class was added (lattice.h)
 * - 01/08/2017: linear algebra class added (linearAlgebra.h), most important algorithms: LLL and short vectors
 * - 01/08/2017: support for boost matrices was added (matrix.h)
 * - 01/08/2017: program options were added (programOptions.h)
 * - 01/08/2017: file logger service was added (log.h)
 * - 01/08/2017: added a service locator (serviceLocator.h)
 * - 01/08/2017: Expected class for better error handling was added (expected.h)
 * - 01/08/2017: tabula rasa
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 */

/*! @file tn.cpp
 *  @brief The main driver for two neighbours.
 *
 *  This file contains the TwoNeighbours class which drives the entire application.
 *
 *  @author Gilles Bellot
 *  @bug 	No known bugs.
 *
 */

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// boost ublas
#include <boost/numeric/ublas/io.hpp>						// input and output for ublas objects
#include <boost/numeric/ublas/symmetric.hpp>				// symmetrical matrices

// CUDA includes
#include <cuda_runtime.h>									// cuda runtim library

// bell0bytes util
#include "Headers/Utilities/serviceLocator.h"				// service locator
#include "Headers/Utilities/programOptions.h"				// program options
#include "Headers/Utilities/expected.h"						// expected error handling

// bell0bytes mathematics
#include "Headers/Mathematics/lattice.h"					// definition of lattices

// DEFINITIONS //////////////////////////////////////////////////////////////////////////

/*!
 * @brief Main driver for the two neighbours program.
 *
 * At initialisation, this class creates and registers all the services.<br>
 * This class also handles exceptions and tries to clean up and print the error message if something unexpected happens.
 */
class TwoNeighbours
{
private:
	bool hasStarted;						//!< True if and only if the application was successfully started.
	bool hasFileLogger;						//!< True if and only if a file logger is available.

	void printStartingLog();				// prints a starting log with CPU and GPU information

public:
	// constructor
	TwoNeighbours();						// the default constructor sets the two member booleans to false and calls the init function
	~TwoNeighbours();						// the destructor destroys

	util::Expected<void> init(int argc, char** argv);		// initialises the application; creates and registers services
	util::Expected<void> run();								// runs the two neighbour algorithm
	void shutdown(util::Expected<void>* result = NULL);		// releases memory and shuts the application down, reports error if the game was shut down by an error

};

// FUNCTIONS ////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Main Function ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
*   @brief The main entry point of the application.
*
*   @param argc The number of arguments given.
*   @param argv The actual arguments, used to create the program options.
*   @return int (0: no error; -1 error)
*/
int main(int argc, char** argv)
{
	TwoNeighbours tn;
	util::Expected<void> initialization = tn.init(argc, argv);
	if(initialization.wasSuccessful())
	{
		// run HN
		util::Expected<void> result = tn.run();

		// shut down
		tn.shutdown(&result);

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
////////////////////////////// Run the Algorithm/////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
*   @brief This function runs the actual program.
*
*   First, the input gram matrix is read from the console or a file.<br>
*   Then a lattice is created from the gramian and relevant properties of the lattice are computed.<br>
*   Once all that is done, the two neighbours algorithm is started.

*   @return An empty Expected<void>, or a filled one with an exception if an error occurred.
*/
util::Expected<void> TwoNeighbours::run()
{
	// read gram matrix of starting lattice
	boost::numeric::ublas::symmetric_matrix<mpz_class> A(1,1);
	std::cin >> A;

	// create lattice from gramian
	mathematics::Lattice L(&A);

	// print starting lattice
	L.print(util::ServiceLocator::getProgOpts()->verboseLevel);

	// start the actual algorithm

	// return success
	return { };
}

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Initialization ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
*   @brief The constructor of the TwoNeighbours class. The constructor sets both of its boolean members to false.
*/
TwoNeighbours::TwoNeighbours() : hasStarted(false), hasFileLogger(false)
{

}

/*!
*   @brief This function initialises the program.
*
*   It creates the file logger and the program options based on the given arguments.
*
*   @param argc The number of arguments given.
*   @param argv The actual arguments given via the console.
*   @return An empty Expected<void>, or a filled one with an exception if an error occurred.
*/
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
/*!
*   @brief This function used the file logger to print a starting log.
*
*   Information about the CPU, the GPU and the loaded libraries are printed.
*
*   @return void
*/
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
/*!
*   @brief The destructor destroys! But in this case, it simply does nothing.
*/
TwoNeighbours::~TwoNeighbours()
{}

/*!
*   @brief This function shuts down the application.
*
*   It tries to clean up everything and to free used memory.<br>
*   It also handles exceptions: If an error is encountered, the program still tries to clean up and prints the actual error message that lead to the abortion.
*
*   @param result A pointer to a void Expected which stores an exception if an error was encountered. If this is NULL, then everything went smoothly.
*   @return void
*/
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

