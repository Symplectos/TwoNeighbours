/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	programOptions.cpp
 * @date	3/11/2015 - Dortmund - Germany
 *
 * @brief	Implementation of the ProgramOptions class defined in programOptions.h
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

#include "../../Headers/Utilities/programOptions.h"
#include "../../Headers/Multi-Precision/mpreal.h"
#include "../../Headers/Utilities/serviceLocator.h"
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <stdio.h>
#include <iostream>
#include <mpfr.h>

namespace util
{
/*!
 *   @brief The constructor.
 *
 *   @section Description
 *
 *   The constructor sets all program options to their default value as follows:
 *    - verboseLevel: 0
 *    - debugMode: false
 *    - useCUDA: false
 *    - outputPrec: 10
 *    - mpfrPrec: 10
 *
 *	 It then creates the program options, based on the command line parameters given, by calling the createCommandLineParameters function.<br>
 *	 If an error is encountered, the constructor throws a std::runtime_error exception.
 *
 *   @param argc The number of arguments given.
 *   @param argv The actual arguments, used to create the program options.
 */
ProgramOptions::ProgramOptions(int argc, char** argv) : verboseLevel(0), debugMode(false), useCUDA(false), outputPrec(10), mpfrPrec(10)
{
	if(!this->createCommandLineParameters(argc, argv).wasSuccessful())
		throw std::runtime_error("Unable to create program options!");
}

/*!
 *   @brief The destructor actually has nothing to do. Shared pointer ftw!
 */
ProgramOptions::~ProgramOptions()
{ }

/*!
 *   @brief Creates the program options.
 *
 *   @section Description
 *
 *   Using boost::program_options, this function creates all the available options to be displayed by passing the "--help" parameter.
 *   At program start, the program options are created and set according to the specified command line parameters.
 *
 *   @param argc The number of arguments given.
 *   @param argv The actual arguments, used to create the program options.
 *   @return Expected<void> As always, the Expected is empty if everything went smoothly, else, a nasty exception is safely stored inside the Expected.
 */
Expected<void> ProgramOptions::createCommandLineParameters(int argc, char** argv)
{
	try
	{
		// declare supported options

		// general program options
		boost::program_options::options_description descMain("General options and settings");
		descMain.add_options()("help", "show help")
						("verboseLevel", boost::program_options::value<int>(), "sets arg (unsigned int) as verbose level - default: 0")
			    		("debug", "debug mode")
						("outputPrec", boost::program_options::value<int>(), "sets arg (int) as output precision - default: 10");

		boost::program_options::options_description descLib("Library options and settings");
		descLib.add_options()	("mpfrPrec", boost::program_options::value<int>(), "sets arg (int) as the precision for mpfr - default: 10")
								("cuda", "enable CUDA algorithms - default: false");

		// declare all options instance which will include all the visible options
		boost::program_options::options_description visible("All allowed options");
		visible.add(descMain);
		visible.add(descLib);

		// declare options instance which will also include hidden options
		boost::program_options::options_description all("All allowed options");
		all.add(descMain);
		all.add(descLib);

		boost::program_options::variables_map vm;
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, all), vm);
		boost::program_options::notify(vm);

		if(vm.count("help"))
		{
			std::cout << visible << std::endl;
			exit(1);
		}

		if(vm.count("cuda"))
			this->useCUDA = true;

		if(vm.count("debug"))
			this->debugMode = true;

		if(vm.count("verboseLevel"))
			this->verboseLevel = vm["verboseLevel"].as<int>();

		// set mpfr precision
		if(vm.count("mpfrPrec"))
		{
			mpfr::mpreal::set_default_prec(mpfr::digits2bits(vm["mpfrPrec"].as<int>()));
			this->mpfrPrec = vm["mpfrPrec"].as<int>();
		}
		else
			mpfr::mpreal::set_default_prec(mpfr::digits2bits(this->mpfrPrec));

		// set std::cout precision
		if (vm.count("outputPrec"))
		{
			std::cout.precision(vm["outputPrec"].as<int>());
			this->outputPrec = vm["outputPrec"].as<int>();
		}
		else
			std::cout.precision(this->outputPrec);
	}
	catch (...)
	{
		return std::runtime_error("Unable to create program options!");
	}

	// return suceess
	return {};
}
}
