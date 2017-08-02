#pragma once

/************************************************************************************//**
* @author	Gilles Bellot
* @file		expected.h
* @date		01/07/2017 - Dortmund - Germany
*
* @brief	Functional error and exception handling.
*
* @section Description
* Based on the talk: C++ and Beyond 2012: Andrei Alexandrescu - Systematic Error Handling in C++.<br>
* See <a href="https://bell0bytes.eu/expected/">my personal website</a> for further details.
*
* @section History
*
* @version	1.0.2.0
* @bug 	No known bugs.
*
* @copyright	Gilles Bellot @ TU Dortmund
****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// exception handling
#include <exception>
#include <stdexcept>

// atomic
#include <atomic>

namespace util
{
	// CLASSES //////////////////////////////////////////////////////////////////////////////
	/*!
	 * @brief The Expected class for advanced error and exception handling.
	 *
	 * Based on the talk: C++ and Beyond 2012: Andrei Alexandrescu - Systematic Error Handling in C++.<br>
     * See <a href="https://bell0bytes.eu/expected/">my personal website</a> for further information.
	*/
	template<class T>
	class Expected
	{
	protected:
		/*!
		*   @union Nameless Union to define an Expected.
		*
		*   There is either a valid result, or an error message stored in the exception pointer.
		*/
		union
		{
			T result;							//!< The expected result.
			std::exception_ptr spam;			//!< The error message if no result could be computed.
		};										//!< The union.

		bool gotResult;							//!< True if and only if a valid result is available.

		/*!
		*   @brief Simple default empty constructor.
		*/
		Expected() : gotResult(true) {}			// a simple protected constructor

	public:
		// constructors and destructor

		/*!
		*   @brief Creates an Expected from a valid result.
		*
		*   @param r A result of type T, passed by reference.
		*/
		Expected(const T& r) : result(r), gotResult(true) {}

		/*!
		*   @brief "Moves" an Expected from a valid result.
		*
		*   @param r Address of a reference to a valid result of type T
		*/
		Expected(T&& r) : result(std::move(r)), gotResult(true) {}

		/*!
		*   @brief A copy constructor.
		*
		*   @param e A reference to the constant Expected to be copied.
		*/
		Expected(const Expected& e) : gotResult(e.gotResult)
		{
			if (gotResult)
				new(&result) T(e.result);
			else
				new(&spam) std::exception_ptr(e.spam);
		}

		/*!
		*   @brief Moves from another Expected.
		*
		*   @param e The address of a reference to an Expected.
		*/
		Expected(Expected&& e) : gotResult(e.gotResult)
		{
			if (gotResult)
				new(&result) T(std::move(e.result));
			else
				new(&spam) std::exception_ptr(std::move(e.spam));
		}

		/*!
		*   @brief An empty destructor.
		*/
		~Expected() {}

		/*!
		*   @brief Swaps two Expected.
		*
		*   @param e The reference to the Expected to swap with.
		*   @return void
		*/
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

		/*!
		*   @brief Creates an Expected from a standard exception.
		*
		*   @param e Reference to a constant exception of type E.
		*/
		template<typename E>
		Expected<T>(E const& e) : spam(std::make_exception_ptr(e)), gotResult(false) { }

		/*!
		*   @brief Creates an Expected from a standard exception.
		*
		*   @param exception A reference to a constant exception.
		*   @return Expected<T> The newly created Expected.
		*/
		template<class E>
		static Expected<T> fromException(const E& exception)
		{
			if (typeid(exception) != typeid(E))
				throw std::invalid_argument("slicing detected!\n");
			return fromException(std::make_exception_ptr(exception));
		}

		/*!
		*   @brief Creates an Expected from a standard exception pointer.
		*
		*   @param p A standard exception pointer.
		*   @return Expected<T> The newly created Expected.
		*/
		static Expected<T> fromException(std::exception_ptr p)
		{
			Expected<T> e;
			e.gotResult = false;
			new(&e.spam) std::exception_ptr(std::move(p));
			return e;
		}

		/*!
		*   @brief Creates an Expected from the current exception.
		*   @return Expected<T> The newly created Expected.
		*/
		static Expected<T> fromException()
		{
			return fromException(std::current_exception());
		}

		// operator overload

		/*!
		*   @brief Copy assignment.
		*
		*   @param e A reference to the constant Expected to be copied.
		*   @return Expected<T> The newly created Expected.
		*/
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

		/*!
		*   @brief Returns true if and only if a valid result is available.
		*
		*   @return bool
		*/
		bool isValid() const { return gotResult; }

		/*!
		*   @brief Returns true if and only if a valid result is available.
		*
		*   @return bool
		*/
		bool wasSuccessful() const { return gotResult; }

		/*!
		*   @brief Returns a valid result or re-throws the error that denied the result from being calculated.
		*
		*   @return T A valid result of type T, is possible.
		*/
		T& get()
		{
			if (!gotResult)
				std::rethrow_exception(spam);
			return result;
		}

		/*!
		*   @brief Returns a valid result or re-throws the error that denied the result from being calculated.
		*
		*   @return T A valid result of type T, is possible.
		*/
		const T& get() const
		{
			if (!gotResult)
				std::rethrow_exception(spam);
			return result;
		}

		// probe for exception

		/*!
		*   @brief Returns true if and only if an exception is stored in the Expect.
		*
		*   @return bool
		*/
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

	/*!
	 * @brief The void Expected class for advanced error and exception handling.
	 *
	 * Based on the talk: C++ and Beyond 2012: Andrei Alexandrescu - Systematic Error Handling in C++.<br>
	    * See <a href="https://bell0bytes.eu/expected/">my personal website</a> for further information.
	*/
	template<>
	class Expected<void>
	{
		std::exception_ptr spam;												//!< The actual error.

	public:
		// constructors and destructor

		/*!
		*   @brief Construct from a reference to a constant type E.
		*
		*   @param e The reference to the constant type E to be constructed from.
		*/
		template <typename E>
		Expected(E const& e) : spam(std::make_exception_ptr(e)) { }

		/*!
		*   @brief The copy constructor.
		*
		*   @param e The reference to the constant Expected to be constructed from,.
		*   @return void
		*/
		template<typename T>
		Expected(const Expected<T>& e)
		{
			if (!e.gotResult)
				new(&spam) std::exception_ptr(e.spam);
		}

		/*!
		*   @brief The move constructor.
		*
		*   @param e The address of a reference to a void Expected.
		*/
		Expected(Expected&& e) : spam(std::move(e.spam)) { }

		/*!
		*   @brief The empty constructor.
		*/
		Expected() : spam() {}

		// operator overload

		/*!
		*   @brief The assignment operator.
		*
		*   @param e A reference to a const Expected to be copied from.
		*   @return An Expected.
		*/
		Expected& operator=(const Expected& e)
		{
			if (!e.isValid())
				this->spam = e.spam;
			return *this;
		};

		// getters

		/*!
		*   @brief Returns true if and only if the computation was successful.
		*
		*   @return bool
		*/
		bool isValid() const { return !spam; }

		/*!
		*   @brief Returns true if and only if the computation was successful.
		*
		*   @return bool
		*/
		bool wasSuccessful() const { return !spam; }

		/*!
		*   @brief Re-throws the error if the computation was not successful.
		*
		*   @return void
		*/
		void get() const { if (!isValid()) std::rethrow_exception(spam); }

		/*!
		*   @brief Suppressed error checking.
		*
		*   @return void
		*/
		void suppress() {}
	};
}
