/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	groupOperations.cpp
 * @date	07/08/2017 - Dortmund - Germany
 *
 * @brief	Implementation of group operation algorithms.
 *
 * @version	0.0.0.1
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// bell0bytes util
#include "../../Headers/Utilities/serviceLocator.h"

// bell0bytes mathematics
#include "../../Headers/Mathematics/groupOperation.h"

// FUNCTIONS ////////////////////////////////////////////////////////////////////////////
namespace mathematics
{

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Action /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

/*!
 *   @brief Group operation.
 *
 *   @param A The element of \f$GL(V)\f$ operating on a subset of the vector space.
 *   @param i The index of the short vector under action.
 *   @param indicesOfSuitableShortVectors The integers in this vector specify the indices of the suitable short vectors in the lattice.
 *	 @param L Constant pointer to the constant lattice under consideration.
 *
 *   @return The index of the image of the vector under the action of A, or an exception if the image was not found in the set of short vectors.
 */
util::Expected<int> GroupOperation::operate(const boost::numeric::ublas::matrix<mpz_class>* const A, const unsigned int i, const std::vector<unsigned int>* const indicesOfSuitableShortVectors, const Lattice* const L)
{
	unsigned int dim = L->gram.size1();
	boost::numeric::ublas::vector<mpz_class> v(dim), w(dim);

	// get the actual vector from the list of short vectors of the lattice
	for(unsigned int j=0; j<dim; j++)
		v(j) = L->shVecs[indicesOfSuitableShortVectors->at(i)].second[j];

	// operate
	w = boost::numeric::ublas::prod(v, (*A));

	// find w in the list of suitable short vectors
	return L->find(&w, indicesOfSuitableShortVectors);
}


}
