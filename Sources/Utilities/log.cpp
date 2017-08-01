// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes
#include "../../Headers/Utilities/log.h"

namespace util
{
	// FUNCTIONS ////////////////////////////////////////////////////////////////////////////
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

	void FileLogPolicy::closeOutputStream()
	{
		outputStream.close();
	}

	void FileLogPolicy::write(const std::string& msg)
	{
		// add the message to the stream
		outputStream << msg << std::endl;
	}
}
