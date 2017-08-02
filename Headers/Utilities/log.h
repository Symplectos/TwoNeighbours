#pragma once

/****************************************************************************************
* Author:	Gilles Bellot
* Date:		13/09/2016 - Lenningen - Luxembourg
*
* Desc:		Event Logger
*
* History:	- 01/07/2017: fixed a memory leak
*			- 02/07/2017: added overloaded print function to take a string
*
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
/** \addtogroup Utility
 *  Utility classes are nice little helper classes to make life a tiny little bit easier.
 *  @{
 */
	// CLASSES //////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////// LOG POLICIES //////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	//! Virtual abstract class to act as an interface to open and close streams (and to write to them).

	//! Customizability is achieved by a purely abstract class called logging policy. A logging policy defines where messages will be printed to.
	//! Now, for example, to create a file logger policy, it is enough to simply inherit from LogPolicyInterface and to specify a file on the hard drive to write out to.
	//! Check <a href="https://bell0bytes.eu/thread-safe-logger/">my personal website</a> for further details.


	class LogPolicyInterface
	{
	public:
		virtual ~LogPolicyInterface() = 0;

		virtual bool openOutputStream(const std::string& name) = 0;
		virtual void closeOutputStream() = 0;
		virtual void write(const std::string& msg) = 0;
	};
	inline LogPolicyInterface::~LogPolicyInterface() {}

	//! Implementation of a policy to write to a file on the hard drive.
	class FileLogPolicy : public LogPolicyInterface
	{
	private:
		std::ofstream outputStream;

	public:
		FileLogPolicy() : outputStream() {};
		~FileLogPolicy() { };

		// member functions
		bool openOutputStream(const std::string& filename) override;	//!< Opens the file given by filename on the hard drive.
		void closeOutputStream() override;								//!< Closes the file.
		void write(const std::string& msg) override;					//!< Writes the message stored in msg to the output stream.
	};

	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// MESSAGE TYPES ////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	// message types
	enum SeverityType
	{ info = 0, debug, warning, error, config };

	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// LOGGER ///////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	// define the class
	template<typename LogPolicy>
	class Logger;

	//! The logging daemon.

	//! To output the contents of a stream buffer, the logger class uses a daemon: The running thread is locked and as long as the daemon is alive, it outputs the elements of the stream buffer.
	//! To start, the timed_mutex, timedMutex (a timed_mutex protects shared data from being simultaneously accessed by multiple threads) of the Logger class is locked by using the unique_lock with defer_lock, that is, the mutex is not immediately locked on construction, but it will be locked soon.
	//! The currently running thread is then put to sleep by the sleep_for function, which simply blocks the execution of the current thread for at least the specified duration.
	//! Once the thread is fast asleep, and there is actually data on the log buffer, attempts to lock the mutex are started. The mutex will be blocked for the length of its supposed slumber, or until the lock is acquired. If a lock can not be acquired at the moment, the thread is allowed to wake and continue safely on its journey (until captured again). If the lock succeeds, the content of the log buffer is written using the specified logging policy.
	template<typename LogPolicy>
	void loggingDaemon(Logger<LogPolicy >* logger)
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

	//! The actual logger class.

	//! This class is to be instantiated with a specific log policy.
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
		Logger(const std::string& name);				//!< The constructor creates a logger to write to the specified file.
		~Logger();										//!< The destructor destroys!

		void setThreadName(const std::string& name);	//!< Sets a human-readable name for the current thread.

		template<SeverityType severity>
		void print(std::stringstream stream);			//!< Prints a message (varies based on the severity level).
		template<SeverityType severity>
		void print(std::string msg);					//!< Prints a message (varies based on the severity level).

		template<typename Policy>
		friend void loggingDaemon(Logger<Policy>* logger);	//!< The logging deamon is obviously a friend of loggers.
	};
	/** @}*/
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


	template<typename LogPolicy>
	void Logger<LogPolicy>::setThreadName(const std::string& name)
	{
		threadName[std::this_thread::get_id()] = name;
	}

	template<typename LogPolicy>
	template<SeverityType severity>
	void Logger<LogPolicy>::print(std::stringstream stream)
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

	template<typename LogPolicy>
	template<SeverityType severity>
	void Logger<LogPolicy>::print(std::string msg)
	{
		std::stringstream stream;
		stream << msg.c_str();
		this->print<severity>(std::stringstream(stream.str()));
	}
}
