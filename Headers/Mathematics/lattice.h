#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	lattice.h
 * @date	09/01/2017 - Dortmund - Germany
 *
 * @brief	Algorithms to work with lattices over \f$\mathbb{Z}\f$.
 *
 * ### Bibliography
 *  - [SRLC] Scharlau, R. --- Lattices and Codes
 *	- [SRQF] Scharlau, R. --- Quadratic Forms
 *	- [HB] Hemkemeier, B. --- Algorithmische Konstruktion von Gittern
 *	- [HR] Haviv, I. and Regev, O. --- On the Lattice Isomorphism Problem
 *	- [KM] Kneser, M. --- Quadratische Formen
 *  - [HV1] Hemkemeier, B. and Vallentin, F. --- On the Decomposition of Lattices
 *  - [HV2] Hemkemeier, B. and Vallentin, F. --- Incremental Algorithms for Lattice Problems
 *
 * ### History
 * - 04/08/2017: added doxygen comments and consts
 * - 13/03/2017: fixed an error in decomp
 * - 10/03/2017: fixed an error in sucmin
 * - 06/02/2017: moved the automorphism algorithms to isomAutom.h
 * - 21/01/2017: fixed an error in the Regev-Haviv-algorithm (checking candidates for bijectivity)
 * - 20/01/2017: speed improvements to the naive automorphism algorithm
 * - 19/01/2017: added naive automorphism group algorithm for lattices with lambda_1 = lambda_n
 * - 19/01/2017: successive minima added
 * - 18/01/2017: added computation of the dual lattice
 * - 17/01/2017: added determinant computation (with LU-factorization)
 * - 09/01/2017: added constructor from gram matrices (LLL and short vectors) and copy constructor
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// c++ includes
#include <valarray>

// boost includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

// gmp includes
#include <gmp.h>
#include <gmpxx.h>
#include "../../Headers/Multi-Precision/mpreal.h"

// DEFINITIONS //////////////////////////////////////////////////////////////////////////

// forward definitions
class Autom;
class Isom;

/*!
 *   @brief This namespace defines mathematical functions.
 */
namespace mathematics
{

/*!
 * @brief A class to work with lattices.
 *
 * This class holds all information needed to work with lattices over \f$\mathbb{Z}\f$.
 * This is work in progress.
 */
class Lattice
{
private:
	// members
	boost::numeric::ublas::symmetric_matrix<mpz_class> gram;					//!< The Gram matrix of the lattice.
	boost::numeric::ublas::symmetric_matrix<mpz_class> dualGram;				//!< The dual of the lattice.
	std::vector<std::vector<std::valarray<mpz_class> > > decomp;				//!< The decomposition of the lattice (see [KM]).
	unsigned long nShVecs;														//!< The number of short vectors (up to highest diagonal entry of the Gram matrix).
	std::vector<std::pair<mpz_class, std::valarray<mpz_class> > > shVecs;		//!< The short vectors of the lattice (up to the highest diagonal entry of the Gram matrix).
	std::valarray<mpz_class> sucMin;											//!< The successive minima of the lattice.
	mpz_class determinant;														//!< The determinant of the Gram matrix (modulo squares, obviously).
	std::vector<boost::numeric::ublas::matrix<mpz_class> > automorphismGroup;	//!< The automorphism group of the lattice.

	// functions
	void dual();																// computes the gramian of the dual of the lattice
	void successiveMinima();													// computes the successive minima of the lattice
	void constructDecomposition();												// constructs an indecomposable decomposition of the lattice

public:
	// constructors and destructor
	Lattice(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const gram);	// creates lattice from gramian
	Lattice(const Lattice&);														// copy constructor
	~Lattice();


	unsigned long getNumberOfShortVectors() const {return this->nShVecs;}	//!< Gets the number of short vectors @return The amount of short vectors.
	std::vector<std::pair<mpz_class, std::valarray<mpz_class> > > getShortVectors() const {return this->shVecs;} //!< Gets the short vectors @return A vector of pairs of multi-precision integers and valarrays of the same type containg the short vectors and their lenghts.

	// elements
	bool contains (const std::valarray<mpz_class>* const v) const;				// returns true if and only if v is an element of the lattice

	// utilities
	void print(const unsigned int flags = 0) const;
};

}
