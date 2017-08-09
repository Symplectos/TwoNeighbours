#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	groupOperations.h
 * @date	07/08/2017 - Dortmund - Germany
 *
 * @brief	Algorithms for group operations.
 *
 * ### Bibliography
 *  - [SR] Scharlau, R. --- Algebra I
 *  - [SC] Sims, S. --- Computational methods in the study of permutation groups
 *  - [SA] Seress, Á. --- Permutation Group Algorithms. Cambridge University Press,
 *
 * ### History
 *  - 08/08/2017: operation added
 *
 * @version	0.0.0.1
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes util
#include "../../Headers/Utilities/expected.h"

// bell0bytes mathematics
#include "../../Headers/Mathematics/lattice.h"

// DEFINITIONS //////////////////////////////////////////////////////////////////////////
namespace mathematics
{
/*!
 * Static class for group operations.
 *
 * Let \f$V\f$ be a finite-dimensional vector space. This class defines algorithms to work with the operation of \f$GL(V)\f$ on a finite subset \f$U \subseteq V\f$,
 * more precisely, \f$GL(V)\f$ operates on a suitable set of short vectors of a lattice.
 *
 */
class GroupOperation
{
private:
	util::Expected<int> operate(const boost::numeric::ublas::matrix<mpz_class>* const A, const unsigned int i, const std::vector<unsigned int>* const indicesOfSuitableShortVectors, const Lattice* const L);

public:
};
}
