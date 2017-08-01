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
	// CLASSES //////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////// LOG POLICIES //////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	// virtual abstract class - interface to open and close streams (and to write to them)
	class LogPolicyInterface
	{
	public:
		virtual ~LogPolicyInterface() = 0;

		virtual bool openOutputStream(const std::string& name) = 0;
		virtual void closeOutputStream() = 0;
		virtual void write(const std::string& msg) = 0;
	};
	inline LogPolicyInterface::~LogPolicyInterface() {}

	// implementation of a policy to write to a file on the hard drive
	class FileLogPolicy : public LogPolicyInterface
	{
	private:
		std::ofstream outputStream;

	public:
		FileLogPolicy() : outputStream() {};
		~FileLogPolicy() { };

		// member functions
		bool openOutputStream(const std::string& filename) override;
		void closeOutputStream() override;
		void write(const std::string& msg) override;
	};

	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// MESSAGE TYPES ////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	// message types
	enum SeverityType
	{
		info = 0,
		debug,
		warning,
		error,
		config
	};

	/////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// LOGGER ///////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	// define the class
	template<typename LogPolicy>
	class Logger;

	// create the actual logging daemon
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

	// the actual logger class to be instantiated with a specific log policy
	template<typename LogPolicy>
	class Logger
	{
	private:
		unsigned int logLineNumber;							// used to write line numbers
		std::map<std::thread::id, std::string> threadName;	// give a human readable name to each thread
		LogPolicy policy;									// the log policy (i.e. write to file, ...)
		std::timed_mutex writeMutex;						// mutual exclusive writer
		std::vector<std::string> logBuffer;					// the log buffer
		std::thread daemon;									// the actual logging daemon
		std::atomic_flag isStillRunning{ ATOMIC_FLAG_INIT };// lock-free boolean to check whether our daemon is still running or not


	public:
		// constructor and destructor
		Logger(const std::string& name);
		~Logger();

		void setThreadName(const std::string& name);	// set human-readable name for the current thread

		template<SeverityType severity>
		void print(std::stringstream stream);			// print a message (varies based on the severity level)
		template<SeverityType severity>
		void print(std::string msg);

		template<typename Policy>
		friend void loggingDaemon(Logger<Policy>* logger);
	};

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
