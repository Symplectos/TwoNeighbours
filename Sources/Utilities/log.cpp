/************************************************************************************//**
* @author	Gilles Bellot
* @file		log.cpp
*
* @brief	Implements the file logger.
*
* @copyright	Gilles Bellot @ TU Dortmund
****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes
#include "../../Headers/Utilities/log.h"

namespace util
{
	// FUNCTIONS ////////////////////////////////////////////////////////////////////////////
	/*!
	 *   @brief Opens the file on the hard drive specified by the parameter filename.
	 *
	 *   @param filename A reference to a constant string specifying the desired file to open as output stream.
	 *   @return bool (true if and only if the output stream was opened successfully)
	 */
	bool FileLogPolicy::openOutputStream(const std::string& filename)
	{
		// try to open the file
		outputStream.open(filename.c_str(), std::ios_base::binary | std::ios_base::out);

		if(!outputStream.is_open())
			return false;

		// set output precision
		outputStream.precision(20);

		// return success
		return true;
	}

	/*!
	*   @brief Closes the file.
	*   @return void
	*/
	void FileLogPolicy::closeOutputStream()
	{
		outputStream.close();
	}

	/*!
	*   @brief Writes a string to the output stream.
	*
	*   @param msg A reference to a constant string to write to the output stream.
	*   @return void
	*/
	void FileLogPolicy::write(const std::string& msg)
	{
		// add the message to the stream
		outputStream << msg << std::endl;
	}
}
