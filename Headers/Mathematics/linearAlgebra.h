#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	linearAlgebra.h
 * @date	14/10/2016 - Eichlinghofen - Germany
 *
 * @brief	General algorithms for linear algebra and lattices.
 *
 * ### Bibliography
 * - [CH] Cohen, H. --- A Course in Computational Algebraic Number Theory
 * - [HB] Hemkemeier, B. --- Algorithmische Konstruktion von Gittern
 *
 * ### History
 * - 08/08/2017: added comparison functions for boost vectors
 * - 04/08/2017: lllGram and lllBasis now use the Expected idiom for exception handling
 * - 04/08/2017: shortVectors now uses the Expected idiom for exception handling
 * - 04/08/2017: added doxygen style comments
 * - 02/08/2017: added const declarations
 * - 01/08/2017: cleaned some unused variables in the various LLL algorithms
 * - 10/03/2017: fixed an error in rank computation
 * - 10/03/2017: fixed an error in rankGaussBareiss
 * - 28/02/2017: fixed an error in gram schmidt
 * - 28/02/2017: fixed a stupid error in scalar product computation
 * - 28/02/2017: added computation of orthogonal complement
 * - 28/02/2017: added scalar product for mpq
 * - 28/02/2017: added naive gram-schmidt algorithm
 * - 28/02/2017: made basis extension algorithm public
 * - 28/02/2017: extended elementOfLattice algorithm to return mv
 * - 27/02/2017: added a serial algorithm with pruning to compute basis from a generating set, see [HB] chapter 3.5
 * - 27/02/2017: added naive algorithm for row reduced algorithm form
 * - 27/02/2017: added basis extension algorithm
 * - 27/02/2017: fixed an error in elementOfLattice
 * - 27/02/2017: deleted rankZ and echelonForm
 * - 21/02/2017: added elementOf, see [HC] 3.4.1
 * - 06/02/2017: added algorithm to compute a basis from a generating set based on [HB]
 * - 06/02/2017: added test to check whether a vector lies in the linear span of a given vector system
 * - 06/02/2017: added test for linear independence over Z
 * - 06/02/2017: added rank computation over the ring of integers
 * - 06/02/2017: added integral echelon form algorithm
 * - 21/01/2017: unrolled the lllDot functions
 * - 21/01/2017: adopted the det, inverse and rank functions from the MatrixHelper class
 * - 19/01/2017: unrolled the scalar product
 * - 19/01/2017: fixed an error in linear dependency over R algorithm - now uses gauss-bareiss
 * - 18/01/2017: added linear dependency test for vectors in Z^n over R
 * - 18/01/2017: added non-euclidean scalar products
 * - 18/01/2017: deleted QFA class -> cholesky added to LALA
 * - 17/01/2017: changed short vector algorithm to use mpz_class instead of longs (needs speed testing)
 * - 17/01/2017: added fudge variable to the short vectors algorithm (needs more testing) - set to 10^(-10) for now
 * - 17/01/2017: fixed iterations in mLLL
 * - 17/01/2017: fixed several asserts in dot product functions
 * - 13/01/2017: added mLLL for integral generating system (floating points used) - not tested yet
 * - 06/01/2017: added short vectors algorithm with Fincke-Pohst preprocessing
 * - 06/01/2017: cleaned LLL algorithms a little bit
 * - 05/01/2017: lllBasis added
 * - 05/01/2017: changed lll to lllGram
 * - 05/01/2017: added fixed point LLL
 * - 17/11/2016: changed cholesky to use valarrays
 * - 17/11/2016: changed lll to use valarrays (less aliasing)
 * - 17/11/2016: added rational cholesky decomposition
 * - 14/11/2016: fixed a loop error in the LLL algorithm
 * - 03/11/2016: added recursive LLL condition test
 * - 02/11/2016: added swap and red sub-agorithms for integral LLL
 *
 * @todo
 *  - write a more robust fp-LLL (see Schnorr-Euchner) --- see example in /Input/dim6 ...
 *  - add a BKZ algorithm
 *  - write a parallel basis construction algorithm for very large generating sets
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// c++
#include <valarray>

// gnu mp
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
#include "../../Headers/Multi-Precision/mpreal.h"

// boost ublas
#include <boost/numeric/ublas/matrix.hpp>

// bell0bytes util
#include "../../Headers/Utilities/expected.h"

namespace mathematics
{
// DEFINITIONS //////////////////////////////////////////////////////////////////////////

// (L)inear (A)lgebra and (L)attice (A)lgorithms
/*!
 * @brief This static class holds (L)inear (A)lgebra and (L)attice (A)lgorithms.
 *
 * ### Notes
 * - C++'s std::valarrays are usually used to encode vectors in the mathematical sense, whereas std::vectors are used as containers.
 * - Multi-precision integers and rationals are stored in instances of the mpz_class, resp. mpq_class.
 * - Multi-precision floats are handled by Pavel Holoborodko's mpfr::mpreal.
 */
class LALA
{
private:
	static constexpr double SHVEC_ELLIPSOID_EPSILON = 0.0000000001;		//!< A fudge variable used in the computation of short vectors.
	static constexpr float eps = std::numeric_limits<float>::epsilon();	//!< The difference between 1.0 and the next representable value of float type. Used as fudge variable in the fp LLL algorithm.

	// general functions
	static bool comparePairMPZValarray(const std::pair<mpz_class, std::valarray<mpz_class> >& a, const std::pair<mpz_class, std::valarray<mpz_class> >& b);

	// rational cholesky decomposition for positive definite forms
	static void choleskyDecomposition(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::triangular_matrix<mpq_class, boost::numeric::ublas::upper>* const r, std::valarray<mpq_class>* const d);

	// scalar products
	static mpfr::mpreal lllBasisDot(const boost::numeric::ublas::matrix<mpfr::mpreal>* const b, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct=nullptr);															// computes the scalar product of the i-th and j-th column of the matrix b
	static mpfr::mpreal lllBasisDot(const boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct = nullptr);	// computes the scalar product of the i-th column of b and j-th column of bStar
	static mpz_class lllBasisDot(const boost::numeric::ublas::matrix<mpz_class>* const b, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct=nullptr);																// computes the scalar product of the i-th and j-th column of the matrix b
	static mpfr::mpreal lllBasisDot(const boost::numeric::ublas::matrix<mpz_class>* const  b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct=nullptr);		// computes the scalar product of the i-th column of b and j-th column of bStar

	// floating point LLL
	static void lllBasisTestCondition(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, int* const k, const int kmax, const int dim);
	static void lllBasisRed(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const mu, const int k, const int l, const int dim);
	static void lllBasisSwap(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, const int k, const int kmax, const int dim);
	static void lllGramTestCondition(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const b, std::valarray<mpfr::mpreal>* const mu, int* const k, const int kmax, const int dim);
	static void lllGramRed(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const mu, const int k, const int l, const int dim);
	static void lllGramSwap(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const b, std::valarray<mpfr::mpreal>* const mu, const int k, const int kmax, const int dim);

	// integral LLL
	static void lllGramTestCondition(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpz_class>* const d, std::valarray<mpz_class>* const lambda, int* const k, const int kmax, const int dim);
	static void lllGramRed(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpz_class>* const d, std::valarray<mpz_class>* const lambda, const int k, const int l, const int dim);
	static void lllGramSwap(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpz_class>* const d, std::valarray<mpz_class>* const lambda, const int k, const int kmax, const int dim);

	// MLLL
	static void mlllBasisSwap(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, const int k, const int kmax, const int dim, const int nVectors);
	static void mlllBasisRed(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpfr::mpreal>* const mu, const int k, const int l, const int nVectors);
	static void mlllBasisTestCondition(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, int* const k, const int kmax, const int dim, const int nVevctors);

	// short vectors
	static void shortVectorsFinckePohst(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, boost::numeric::ublas::matrix<mpfr::mpreal>* const p, boost::numeric::ublas::matrix<mpfr::mpreal>* const q);	// Fincke-Pohst preprocessing, returns q
	static util::Expected<unsigned long> shortVectorsPrivate(const boost::numeric::ublas::matrix<mpfr::mpreal>* const q, const mpz_class c, std::vector<std::pair<mpz_class, std::valarray<mpz_class> > >* const X);	// short vectors algorithm

public:
	// compare functions
	static int compare(const boost::numeric::ublas::vector<mpz_class>* const v, const boost::numeric::ublas::vector<mpz_class>* const w);
	static int compare(const boost::numeric::ublas::vector<mpz_class>* const v, const std::valarray<mpz_class>* const w);

	// linear independence
	static bool linearlyIndependentR(const std::vector<std::valarray<mpz_class> >* const vectorSystem);	// returns true iff the vector system is linearly independent over R
	static bool linearlyIndependentZ(const std::vector<std::valarray<mpz_class> >* const vectorSystem);	// returns true iff the vector system is linearly independent over Z

	// element of
	static bool inLinearSpanZ(const std::vector<std::valarray<mpz_class> >* const b, const std::valarray<mpz_class>* const v);	// returns true iff v is an element of the vector space spanned by the column vectors in b
	static bool elementOfLattice(const boost::numeric::ublas::matrix<mpq_class>* const m, const std::valarray<mpz_class>* const v, const unsigned int dim, boost::numeric::ublas::vector<mpq_class>* const mv = nullptr); // returns true iff v is an element of the lattice given by the matrix m (m must be the inverse of a matrix of a lattice basis, possibly extended to a vector space basis), returns mv is so desired

	// scalar products
	static mpz_class scalarProduct(const std::valarray<mpz_class>* const x, const std::valarray<mpz_class>* const y, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct = nullptr);	// computes x^t * innerProduct * y
	static mpq_class scalarProduct(const std::valarray<mpq_class>* const x, const std::valarray<mpq_class>* const y, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct = nullptr);

	// gram schmidt
	static void gramSchmidt(std::vector<std::valarray<mpq_class> >* const basis, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g = nullptr);	// gram schmidt

	// rank of matrices
	static void echelonForm(boost::numeric::ublas::matrix<mpq_class>* const m, std::vector<unsigned int>* const basisExtension = nullptr);	// computes row-reduced echelon form of a given matrix, also outputs the non pivot columns, if so desired
	static void echelonFormConst(const boost::numeric::ublas::matrix<mpq_class>* const m, std::vector<unsigned int>* const basisExtension);	// computes basis extensions for partial basis given by the columns of m
	static unsigned long rank(const boost::numeric::ublas::matrix<mpfr::mpreal>* const m);	// computes the rank of m using gaussian elimination (over R)
	static unsigned long rank(const boost::numeric::ublas::matrix<mpq_class>* const m);		// computes the rank of m using gaussian elimination (over Q)
	static unsigned long rank(const boost::numeric::ublas::matrix<mpz_class>* const m);		// computes the rank of m using gaussian elimination (over Q)

	// inverse of matrices
	static void inverse(const boost::numeric::ublas::triangular_matrix<mpfr::mpreal, boost::numeric::ublas::upper>* const a, boost::numeric::ublas::triangular_matrix<mpfr::mpreal, boost::numeric::ublas::upper>* const b); // computes b such that a*b=1
	static bool inverse(const boost::numeric::ublas::matrix<mpfr::mpreal>* const a, boost::numeric::ublas::matrix<mpfr::mpreal>* const b);			// sets b = a^-1; uses LU-factorization
	static bool inverse(const boost::numeric::ublas::matrix<mpq_class>* const a, boost::numeric::ublas::matrix<mpq_class>* const b);

	// determinant of matrices
	static mpfr::mpreal det(const boost::numeric::ublas::matrix<mpfr::mpreal>* const m);					// computes the determinant of the matrix m using LU-factorization
	static mpfr::mpreal det(const boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const m);
	static mpz_class detGaussBareiss(const boost::numeric::ublas::matrix<mpz_class>* const m);				// computes the determinant of the matrix m using the Gauss-Bareiss algorithm
	static mpz_class detGaussBareiss(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const m);

	// basis extensions
	static void basisExtensionZ(boost::numeric::ublas::matrix<mpq_class>* const m);	// extends the columns of m to a vector space basis
	static void orthogonalComplement(const std::vector<std::valarray<mpq_class> >* const latticeBasis, std::vector<std::valarray<mpq_class> >* const orthCompl, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g = nullptr);	// computes the orthogonal complement of a lattice in a vector space

	// basis from generating set as described in [HB] 3.5
	static void basisFromGeneratingSystem(const std::vector<std::valarray<mpz_class> >* const gen, std::vector<std::valarray<mpz_class> >* const basis);

	// floating point LLL as described in [CH] Algorithm 2.6.3
	static util::Expected<void> lllBasis(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const h = nullptr);			// returns -1 iff given vector system was linearly dependent
	static util::Expected<void> lllGram(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h = nullptr);	// returns -1 iff given vector system was linearly dependent

	// fp-mLLL as described in [CH] Algorithm 2.6.8
	static unsigned int mlllBasis(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpz_class>* const h = nullptr, boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct = nullptr);		// returns the rank of the given vector system - computes basis from linearly dependent vectors

	// integral LLL as described in [CH] Algorithm 2.6.7
	static util::Expected<void> lllGram(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h = nullptr); 		// returns -1 iff given vector system was linearly dependent

	// short vectors of L, that is, vectors x \in L with N^2(x) <= c - with Fincke and Pohst preprocessing (LLL reduction)
	static util::Expected<unsigned long> shortVectors(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, const mpz_class c, std::vector<std::pair<mpz_class, std::valarray<mpz_class> > >* const X);
	};

}
