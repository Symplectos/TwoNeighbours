#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	log.h
 * @date	13/09/2016 - Lenningen - Luxembourg
 *
 * @brief	A thread-safe event logger.
 *
 * @section logDescription Description
 *
 * To be able to easier debug a program, a robust, yet lightweight, event logger, such that every class in the engine can provide a trace of its execution in a log file, is needed.
 * Obviously such a logging system must be very resilient, as like the captain of a sinking ship, it must stay "on board" until the very end. It should also be able to write out warnings of different severity levels (warning, debug, errors, ...) to various output channels.
 * Customizability is achieved by a purely abstract class called logging policy. A logging policy defines where messages will be printed to.
 * Now, for example, to create a file logger policy, it is enough to simply inherit from LogPolicyInterface and to specify a file on the hard drive to write out to.
 *
 * Check <a href="https://bell0bytes.eu/thread-safe-logger/">my personal website</a> for further details.
 *
 * @section logHistory History
 *
 *  - 02/07/2017: added overloaded print function to take a string
 *  - 01/07/2017: fixed a memory leak
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// c++ includes
#include <atomic>
#include <thread>
#include <mutex>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

// boost includes
#include <boost/date_time/posix_time/posix_time.hpp>

namespace util
{
// CLASSES //////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// LOG POLICIES //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 * @brief The LogPolicyInterface class is a virtual abstract class to act as an interface to open and close streams (and to write to them).
 */
class LogPolicyInterface
{
public:
	virtual ~LogPolicyInterface() = 0;

	/*!
	 *   @brief This virtual member must be overridden to handle the opening of the output stream.
	 *
	 *   @param name A reference to a constant string specifying the desired output stream (i.e. file location)
	 *   @return bool (true if and only if the output stream was opened successfully)
	 */
	virtual bool openOutputStream(const std::string& name) = 0;

	/*!
	 *   @brief This must be overridden to handle closing the output stream.
	 *   @return void
	 */
	virtual void closeOutputStream() = 0;

	/*!
	 *   @brief This must be overridden to specify how to write to the output stream.
	 *
	 *   @param msg A reference to a constant string to be written to the output stream.
	 *   @return void
	 */
	virtual void write(const std::string& msg) = 0;
};

/*!
 *   @brief The virtual empty default constructor.
 */
inline LogPolicyInterface::~LogPolicyInterface() {}

/*!
 * @brief The FileLogPolicy class is derived from the LogPolicyInterface and acts as an interface to write to files on the hard drive.
 *
 * Check <a href="https://bell0bytes.eu/thread-safe-logger/">my personal website</a> for further details.
 */
class FileLogPolicy : public LogPolicyInterface
{
private:
	std::ofstream outputStream;										//!< The output file stream.

public:
	/*!
	 *   @brief The default empty constructor.
	 */
	FileLogPolicy() : outputStream() {};

	/*!
	 *   @brief The default empty destructor.
	 */
	~FileLogPolicy() { };

	// member functions
	bool openOutputStream(const std::string& filename) override;
	void closeOutputStream() override;
	void write(const std::string& msg) override;
};

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MESSAGE TYPES ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

//! The severity types of the logged messages.
enum SeverityType
{ info = 0, debug, warning, error, config };

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// LOGGER ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// define the class
template<typename LogPolicy>
class Logger;

/*!
 *   @brief The logging daemon.
 *
 *   To output the contents of a stream buffer, the logger class uses a daemon: The running thread is locked and as long as the daemon is alive, it outputs the elements of the stream buffer.
 *   To start, the timed_mutex, timedMutex (a timed_mutex protects shared data from being simultaneously accessed by multiple threads) of the Logger class is locked by using the unique_lock with defer_lock, that is, the mutex is not immediately locked on construction, but it will be locked soon.
 *   The currently running thread is then put to sleep by the sleep_for function, which simply blocks the execution of the current thread for at least the specified duration.
 *   Once the thread is fast asleep, and there is actually data on the log buffer, attempts to lock the mutex are started. The mutex will be blocked for the length of its supposed slumber, or until the lock is acquired. If a lock can not be acquired at the moment, the thread is allowed to wake and continue safely on its journey (until captured again). If the lock succeeds, the content of the log buffer is written using the specified logging policy.
 *
 *   @param logger A constant pointer to a logger with a specified log policy.
 *   @return void
 */
template<typename LogPolicy>
void loggingDaemon(Logger<LogPolicy >* const logger)
{
	// dump the log data if present
	std::unique_lock<std::timed_mutex> lock(logger->writeMutex, std::defer_lock);	// lock the mutex
	do
	{
		std::this_thread::sleep_for(std::chrono::milliseconds{ 50 });
		if (logger->logBuffer.size())
		{
			if (!lock.try_lock_for(std::chrono::milliseconds{ 50 }))
				continue;
			for (auto& elem : logger->logBuffer)
				logger->policy.write(elem);											// write contents of the buffer
			logger->logBuffer.clear();
			lock.unlock();
		}
	} while (logger->isStillRunning.test_and_set() || logger->logBuffer.size());	// do this while while the logger is still active
}


/*!
 * @brief The actual logger class.
 *
 * This class is to be instantiated with a specific log policy.
 */
template<typename LogPolicy>
class Logger
{
private:
	unsigned int logLineNumber;							//!< Used to write line numbers.
	std::map<std::thread::id, std::string> threadName;	//!< A human readable name to each thread.
	LogPolicy policy;									//!< The log policy (i.e. write to file, ...).
	std::timed_mutex writeMutex;						//!< The mutual exclusive writer.
	std::vector<std::string> logBuffer;					//!< The log buffer.
	std::thread daemon;									//!< The logging daemon.
	std::atomic_flag isStillRunning{ ATOMIC_FLAG_INIT };//!< Lock-free boolean to check whether our daemon is still running or not.


public:
	// constructor and destructor
	Logger(const std::string& name);				// the constructor creates a logger to write to the specified file
	~Logger();										// the destructor destroys!

	void setThreadName(const std::string& name);	// sets a human-readable name for the current thread

	template<SeverityType severity>
	void print(const std::stringstream stream);		// prints a message (varies based on the severity level)
	template<SeverityType severity>
	void print(const std::string msg);				// prints a message (varies based on the severity level)

	/*!
	 *   @brief The logging daemon is obviously a friend of every logger.
	 *
	 *   @param logger A constant pointer to the logger that needs the daemon to print its output stream.
	 *   @return void
	 */
	template<typename Policy>
	friend void loggingDaemon(Logger<Policy>* const logger);
};

/*!
 *   @brief The constructor.
 *
 *	The constructor sets all the member variables to their default values and then opens the desired output file on the hard disk. Upon success, the daemon is started.
 *	If an error is encountered, the constructor throws a std::runtime_error exception.
 *
 *   @param name A reference to a constant string specifying the desired output stream on the hard disk.
 *
 *   @exception std::runtime_error thrown when the desired output file could not be opened.
 */
template<typename LogPolicy>
Logger<LogPolicy>::Logger(const std::string& name) : logLineNumber(0), threadName(), policy(), writeMutex(), logBuffer()
{
	if(policy.openOutputStream(name))
	{
		isStillRunning.test_and_set();		// mark the logging daemon as running
		daemon = std::move(std::thread{ loggingDaemon<LogPolicy>,this });
	}
	else
		throw std::runtime_error("Unable to open the log file!");
}

/*!
 *   @brief The destructor.
 *
 *   The destructor lets the daemon join the main thread, clears the thread names and the log buffer, before finally closing the output file.
 */
template<typename LogPolicy>
Logger<LogPolicy>::~Logger()
{
	// terminate the daemon by clearing the still running flag and letting it join to the main thread
	isStillRunning.clear();
	daemon.join();

	// clear the thread name map
	threadName.clear();
	std::map<std::thread::id, std::string>().swap(threadName);

	// clear the log vector
	logBuffer.clear();
	logBuffer.shrink_to_fit();

	// close the output stream
	policy.closeOutputStream();
}

/*!
 *   @brief Sets the name of the current thread to a human-readable form.
 *
 *   @param name A reference to a constant string specifying the name of the current thread.
 *   @return void
 */
template<typename LogPolicy>
void Logger<LogPolicy>::setThreadName(const std::string& name)
{
	threadName[std::this_thread::get_id()] = name;
}

/*!
 *   @brief Writes the content of a stringstream to the output file.
 *
 *   @param stream A constant stringstream to be added to the output file stream.
 *   @return void
 */
template<typename LogPolicy>
template<SeverityType severity>
void Logger<LogPolicy>::print(const std::stringstream stream)
{
	std::stringstream logStream;

	// all severity types but the config type allow custom formatting
	if(!(severity == SeverityType::config))
	{
		// header - line number and date
		logStream << logLineNumber++ << " < " << boost::posix_time::second_clock::local_time() << "\t";

		// write warning level
		switch (severity)
		{
		case SeverityType::info:
			logStream << "INFO:    ";
			break;
		case SeverityType::debug:
			logStream << "DEBUG:   ";
			break;
		case SeverityType::warning:
			logStream << "WARNING: ";
			break;
		case SeverityType::error:
			logStream << "ERROR:   ";
			break;
		};

		// write thread name
		logStream << threadName[std::this_thread::get_id()] << ":\t";
	}

	// write actual message
	logStream << stream.str();
	std::lock_guard<std::timed_mutex> lock(writeMutex);
	logBuffer.push_back(logStream.str());
}

/*!
 *   @brief Writes a string to the output file stream.
 *
 *   @param msg A constant string to be written to the output file stream.
 *   @return void
 */
template<typename LogPolicy>
template<SeverityType severity>
void Logger<LogPolicy>::print(const std::string msg)
{
	std::stringstream stream;
	stream << msg.c_str();
	this->print<severity>(std::stringstream(stream.str()));
}
}
