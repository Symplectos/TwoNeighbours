#pragma once

/*
 * File:   lattices.h
 * Author: Gilles Bellot
 * Date:   09/01/2017 - Dortmund - Germany
 *
 * Description: class to define lattces over Z
 *
 * Bibliography:	- [SRLC] Scharlau, R. --- Lattices and Codes
 *					- [SRQF] Scharlau, R. --- Quadratic Forms
 *					- [HB] Hemkemeier, B. --- Algorithmische Konstruktion von Gittern
 *					- [HR] Haviv, I. and Regev, O. --- On the Lattice Isomorphism Problem
 *					- [KM] Kneser, M. --- Quadratische Formen
 *
 * History:		- 09/01/2017: added constructor from gram matrices (LLL and short vectors) and copy constructor
 * 				- 17/01/2017: added determinant computation (with LU-factorization)
 * 				- 18/01/2017: added computation of the dual lattice
 * 				- 19/01/2017: successive minima added
 * 				- 19/01/2017: added naive automorphism group algorithm for lattices with lambda_1 = lambda_n
 * 				- 20/01/2017: speed improvements to the naive automorphism algorithm
 * 				- 21/01/2017: fixed an error in the Regev-Haviv-algorithm (checking candidates for bijectivity)
 * 				- 06/02/2017: moved the automorphism algorithms to isomAutom.h
 * 				- 10/03/2017: fixed an error in sucmin
 * 				- 13/03/2017: fixed an error in decomp
 *
 *
 * ToDo:	- clean up decomp and sucmin algorithms
 * 			- see magma for good output
 * 			- add better algorithm for automorphism groups and isometry testing
 *
*/

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

namespace mathematics
{

class Lattice
{
private:
	// members
	boost::numeric::ublas::symmetric_matrix<mpz_class> gram;					// gramian of the lattice
	boost::numeric::ublas::symmetric_matrix<mpz_class> dualGram;				// the dual of the lattice
	std::vector<std::vector<std::valarray<mpz_class> > > decomp;				// the decomposition of the lattice (see [KM])
	unsigned long nShVecs;														// number of short vectors (up to highest diagonal entry of the gramian)
	std::vector<std::pair<mpz_class, std::valarray<mpz_class> > > shVecs;		// short vectors of the lattice (up to the highest diagonal entry of the gramian)
	std::valarray<mpz_class> sucMin;											// successive minima of the lattice
	mpz_class determinant;														// determinant of the gram matrix (modulo squares, obviously)
	std::vector<boost::numeric::ublas::matrix<mpz_class> > automorphismGroup;	// the automorphism group of the lattice

	// functions
	void dual();																// computes the gramian of the dual of the lattice
	void successiveMinima();													// computes the successive minima of the lattice
	void constructDecomposition();												// constructs an indecomposable decomposition of the lattice

public:
	// constructors and destructor
	Lattice(const boost::numeric::ublas::symmetric_matrix<mpz_class>* gram);	// creates lattice from gramian
	Lattice(const Lattice&);													// copy constructor
	~Lattice();

	// getters and setters
	unsigned long getNumberOfShortVectors(){return this->nShVecs;};
	std::vector<std::pair<mpz_class, std::valarray<mpz_class> > > getShortVectors(){return this->shVecs;}

	// elements
	bool contains(std::valarray<mpz_class>* v);									// returns true if and only if v is an element of the lattice

	// utilities
	void print(unsigned int flags = 0);
};

}
