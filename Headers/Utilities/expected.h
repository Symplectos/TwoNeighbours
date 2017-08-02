#pragma once

/****************************************************************************************
* Author:	Gilles Bellot
* Date:		01/07/2017 - Dortmund - Germany
*
* Desc:		functional error and exception handling
*			based on the talk: C++ and Beyond 2012: Andrei Alexandrescu - Systematic Error Handling in C++
*
****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// exception handling
#include <exception>
#include <stdexcept>

// atomic
#include <atomic>

namespace util
{
/** \addtogroup Utility
 *  @{
 */
	// CLASSES //////////////////////////////////////////////////////////////////////////////

	//! The Expected class for advanced error and exception handling.

	//! Based on the talk: C++ and Beyond 2012: Andrei Alexandrescu - Systematic Error Handling in C++.
	//! See <a href="https://bell0bytes.eu/expected/">my personal website</a> for further information.
	template<class T>
	class Expected
	{
	protected:
		union
		{
			T result;							//!< The expected result.
			std::exception_ptr spam;			//!< The error message if no result could be computed.
		};

		bool gotResult;							//!< True if and only if a valid result is available.
		Expected() : gotResult(false) {}		//!< A simple protected constructor.

	public:
		// constructors and destructor
		Expected(const T& r) : result(r), gotResult(true) {}		//!< Creates an Expected with a valid result.
		Expected(T&& r) : result(std::move(r)), gotResult(true) {}	//!< The move constructor when given a valid result.

		//! The copy constructor.
		Expected(const Expected& e) : gotResult(e.gotResult)
		{
			if (gotResult)
				new(&result) T(e.result);
			else
				new(&spam) std::exception_ptr(e.spam);
		}
		//! Move constructor when given another Expected.
		Expected(Expected&& e) : gotResult(e.gotResult)
		{
			if (gotResult)
				new(&result) T(std::move(e.result));
			else
				new(&spam) std::exception_ptr(std::move(e.spam));
		}
		~Expected() {}												//!< The destructor destroys!

		//! Swaps two Expected.
		void swap(Expected& e)
		{
			if (gotResult)
			{
				if (e.gotResult)
					std::swap(result, e.result);
				else
				{
					auto t = std::move(e.spam);
					new(&e.result) T(std::move(result));
					new(&spam) std::exception_ptr;
					std::swap(gotResult, e.gotResult);
				}
			}
			else
			{
				if (e.gotResult)
					e.swap(*this);
				else
					spam.swap(e.spam);
			}
		}

		// creating expect from exceptions
		template<typename E>
		Expected<T>(E const& e) : spam(std::make_exception_ptr(e)), gotResult(false) { }	//!< Creates an Expect from a standard exception.

		//! Creates an Expect from a standard exception.
		template<class E>
		static Expected<T> fromException(const E& exception)
		{
			if (typeid(exception) != typeid(E))
				throw std::invalid_argument("slicing detected!\n");
			return fromException(std::make_exception_ptr(exception));
		}
		//! Creates an Expect from a standard exception pointer.
		static Expected<T> fromException(std::exception_ptr p)
		{
			Expected<T> e;
			e.gotResult = false;
			new(&e.spam) std::exception_ptr(std::move(p));
			return e;
		}
		//! Creates an Expected from the current exception.
		static Expected<T> fromException()
		{
			return fromException(std::current_exception());
		}

		// operator overload
		//! Assignment operator.
		Expected<T>& operator=(const Expected<T>& e)
		{
			gotResult = e.gotResult;
			if (gotResult)
				new(&result) T(e.result);
			else
				new(&spam) std::exception_ptr(e.spam);
			return *this;
		}

		// getters
		bool isValid() const { return gotResult; }											//!< Returns true if and only if a valid result is available.
		bool wasSuccessful() const { return gotResult; }									//!< Returns true if and only if a valid result is available.

		//! Returns a valid result or rethrows the error message denying the result from being calculated.
		T& get()
		{
			if (!gotResult)
				std::rethrow_exception(spam);
			return result;
		}
		//! Const getter.
		const T& get() const
		{
			if (!gotResult)
				std::rethrow_exception(spam);
			return result;
		}

		// probe for exception
		//! Returns true if and only if an exception is stored in the Expect.
		template<class E>
		bool hasException() const
		{
			try
			{
				if (!gotResult)
					std::rethrow_exception(spam);
			}
			catch (const E& object)
			{
				(void)object;
				return true;
			}
			catch (...)
			{

			}
			return false;
		}

		friend class Expected<void>;														//!< The void Expected is a friend of all the other Expected.
	};

	//! The void Expected class for advanced error and exception handling.

	//! Based on the talk: C++ and Beyond 2012: Andrei Alexandrescu - Systematic Error Handling in C++.
	//! See <a href="https://bell0bytes.eu/expected/">my personal website</a> for further information.
	template<>
	class Expected<void>
	{
		std::exception_ptr spam;												//!< The actual error.

	public:
		// constructors and destructor
		//! The default constructor.
		template <typename E>
		Expected(E const& e) : spam(std::make_exception_ptr(e)) { }
		//! Copy constructor.
		template<typename T>
		Expected(const Expected<T>& e)
		{
			if (!e.gotResult)
				new(&spam) std::exception_ptr(e.spam);
		}

		Expected(Expected&& o) : spam(std::move(o.spam)) { }					//!< Move constructor.
		Expected() : spam() {}													//!< Empty constructor.

		// operator overload
		//! Assignment operator.
		Expected& operator=(const Expected& e)
		{
			if (!e.isValid())
				this->spam = e.spam;
			return *this;
		};

		// getters
		bool isValid() const { return !spam; }									//!< Returns true if and only if the computation was successful.
		bool wasSuccessful() const { return !spam; }							//!< Returns true if and only if the computation was successful.
		void get() const { if (!isValid()) std::rethrow_exception(spam); }		//!< Throws the error if the computation was not successful.
		void suppress() {}														//!< Suppresses error checking.
	};
	/** @}*/
}
