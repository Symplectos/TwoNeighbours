#pragma once

/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	matrix.h
 * @date	13/10/2016 - Dortmund - Germany
 *
 * @brief	Wrapper class for BOOST uBLAS matrices.
 *
 * @section History
 * - 02/07/2017: added doxygen comments
 * - 01/08/2017: added const declarations
 * - 21/01/2017: moved rank, det and inverse function to the linear algebra class
 * - 20/01/2017: added comparison function
 * - 20/01/2017: added orthoginality check with inner product
 * - 19/01/2017: added orthogonality check
 * - 19/01/2017: added gauss bareiss for rank computation
 * - 19/01/2017: added gauss bareiss for determinant computation
 * - 18/01/2017: changed inverse function to const
 * - 18/01/2017: changed print functions to const
 * - 18/01/2017: added det function (with lu)
 * - 18/01/2017: clean up of printMatrix functions
 * - 18/01/2017: deleted printSymmetrixMatrix function
 * - 18/01/2017: deleted read function -> replaced by boosts own read function on std::cin
 * - 18/01/2017: deleted mult function -> replaced by boosts own prod function
 * - 18/01/2017: deleted makeUnit function -> replaced by boosts own assign function
 * - 17/01/2017: added function to delete zero columns from basis matrix
 * - 05/01/2017: deleted readMatrix
 * - 05/01/2017: fixed a template error
 * - 29/11/2016: fixed an error in trig matrix multiplication and inversion
 * - 25/11/2016: added matrix multiplication and inversion for triangular matrices
 * - 17/11/2016: added print support for triangular matrices and vectors
 * - 17/11/2016: added base change (add multiple, multiply) for triangular matrices
 * - 02/11/2016: added base change (swap) for gramians
 * - 17/10/2016: added base change (swap) for matrices
 * - 17/10/2016: added base change (add multiple) for gramians
 * - 15/10/2016: added base change (add multiple) for matrices
 * - 14/10/2016: added support for unit matrix
 * - 13/10/2016: added matrix reading (boost) and printing (human)
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// BOOST uBLAS
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

// BOOST utility
#include <boost/utility/swap.hpp>

// c++ stuff
#include <algorithm>
#include <initializer_list>

namespace mathematics
{
// DEFINITIONS //////////////////////////////////////////////////////////////////////////

/*!
 * @brief The MatrixHelper is a static templated class full of little helper function to facilitate the use of BOOST's uBLAS matrices.
 */
template<class T>
class MatrixHelper
{
private:

public:
	// elementary column operations
	static void colSwap(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j);				// swaps the i-th and j-th column of the matrix m
	static void colAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, T t);	// sets column_i = column_i + t * column_j
	static void colMultiply(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const T t);					// sets col_i *= t
	static std::valarray<T> getColumn(const boost::numeric::ublas::matrix<T>* const m, const unsigned int i);				// returns the i-th column of the matrix m
	static std::valarray<T> getColumn(const boost::numeric::ublas::symmetric_matrix<T>* const m, const unsigned int i);		// returns the i-th column of the matrix m

	// elementary row operations
	static void rowAddMultiple(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const unsigned int j, T t);	// sets row_i = row_i + t * row_j
	static void rowMultiply(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const T t);						// sets row i = t * row_i
	static void rowAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, const T t);										// sets row_i = row_i + t * row_j

	// base change for gram matrices
	static void baseSwap(boost::numeric::ublas::symmetric_matrix<T>* const g, unsigned int i, unsigned int j);							// swaps the i-th and j-th base vector
	static void baseAddMultiple(boost::numeric::ublas::symmetric_matrix<T>* const g, const unsigned int i, const unsigned int j, T t);	// sets: b_i = b_i + q * b_j

	// general utility functions for base or gram matrices
	static void baseDelZeroes(boost::numeric::ublas::matrix<T>* const b, const unsigned int rk);								// deletes zero columns of m
	static void setColumn(boost::numeric::ublas::matrix<T>* const b, const std::valarray<T>* const col, const unsigned int i);	// sets the i-th column of m to col

	// check properties
	static bool isOrthogonal(const boost::numeric::ublas::matrix<T>* const m, const boost::numeric::ublas::symmetric_matrix<T>* const g = nullptr);	// returns true if and only if m^t * m=id or m^t * g \cdot m=g

	// compare
	static bool areEqual(const boost::numeric::ublas::matrix<T>* const m1, const boost::numeric::ublas::matrix<T>* const m2);	// returns true if and only if both matrices are equal

	// input and output
	static void printMatrix(const boost::numeric::ublas::matrix<T>* const m);											// prints a matrix in the old TN format
	static void printMatrix(const boost::numeric::ublas::symmetric_matrix<T>* const m);									// prints a symmetric matrix in the old TN format
	static void printMatrix(const boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m); 	// prints an upper triangular matrix in the old TN format
	static void printDiagonalVector(const std::valarray<T>* const d);													// prints a diagonal vector
};
/** @}*/

// IMPLEMENTATION ///////////////////////////////////////////////////////////////////////

// base change for matrices
/*!
 *   @brief Swaps two columns of a matrix.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a matrix of type T.
 *   @param i Constant unsigned int specifying a column of the given matrix: \f$0 \leq i \leq n\f$.
 *   @param j Constant unsigned int specifying a column of the given matrix. \f$0 \leq j \leq n\f$.
 *
 *   Swaps the i-th and j-th column of the matrix m, where \f$n\f$ is the number of columns of m.<br>
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::colSwap(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j)
{
	assert(i<m->size2() && j<m->size2() && i>=0 && j>=0);

	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnI(*m, i);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnJ(*m, j);

	boost::swap(columnI, columnJ);
}

/*!
 *   @brief Adds a multiple of a column to another column.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a matrix of type T.
 *   @param i Constant unsigned int specifying a column of the given matrix: \f$0 \leq i \leq n\f$.
 *   @param j Constant unsigned int specifying a column of the given matrix: \f$0 \leq j \leq n\f$.
 *   @param t Constant element of type T.
 *
 *	 Let \f$n\f$ denote the number of columns of the matrix and \f$m_i\f$, \f$m_j\f$ the i-th, resp. j-th column of the matrix m.<br>
 *	 This function then sets \f$ m_i = m_i + t \cdot m_j\f$.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::colAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, const T t)
{
	assert(i < m->size2() && j < m->size2() && i>=0 && j>=0);

	// get the columns in question
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnI(*m, i);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnJ(*m, j);

	// add multiple
	for(unsigned int k=0; k<columnI.size(); k++)
		columnI(k)+=t*columnJ(k);
}

/*!
 *   @brief Multiplies a column of a matrix by a constant.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a matrix of type T.
 *   @param i Constant unsigned int specifying a column of the given matrix: \f$0 \leq i \leq n\f$.
 *   @param t Constant element of type T.
 *
 *	 Let \f$n\f$ denote the number of columns of the matrix and \f$m_i\f$ the i-th column of the matrix m.<br>
 *	 This function then sets \f$ m_i = t \cdot m_i\f$.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::colMultiply(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const T t)
{
	assert(i < m->size2() && i>=0);

	// get the column in question
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnI(*m, i);

	// multiply
	for(unsigned int k=0; k<columnI.size(); k++)
		columnI(k)*=t;
}

// base change for triangular matrices
/*!
 *   @brief Adds a multiple of a row to another row.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to an upper triangular matrix of type T.
 *   @param i Constant unsigned int specifying a row of the given matrix: \f$0 \leq i \leq n\f$.
 *   @param j Constant unsigned int specifying a row of the given matrix: \f$0 \leq j \leq n\f$.
 *   @param t Constant element of type T.
 *
 *	 Let \f$n\f$ denote the number of rows of the given matrix and \f$m_i\f$, \f$m_j\f$ the i-th, resp. j-th row of the matrix m.<br>
 *	 This function then sets \f$ m_i = m_i + t \cdot m_j\f$.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::rowAddMultiple(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const unsigned int j, const T t)
{
	assert(i < m->size1() && j < m->size1() && i>=0 && j>=0);
	unsigned int n = m->size1();

	// add multiple
	if(i==j)
		for(unsigned int k=i; k<n; k++)
			(*m)(i,k)+=t*(*m)(j,k);
	else if(i<j)
		// lower triangular, thus row i has "more" elements than row j
		for(unsigned int k=j; k<n; k++)
			(*m)(i,k)+=t*(*m)(j,k);
	else if(i>j)
		// lower triangular, thus row j has "more" elements than row i
		for(unsigned int k=i; k<n; k++)
			(*m)(i,k)+=t*(*m)(j,k);
}

/*!
 *   @brief Multiplies a row of a matrix by a constant.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to an upper triangular matrix of type T.
 *   @param i Constant unsigned int specifying a row of the given matrix: \f$0 \leq i \leq n\f$.
 *   @param t Constant element of type T.
 *
 *	 Let \f$n\f$ denote the number of rows of the matrix and \f$m_i\f$ the i-th row of the matrix m.<br>
 *	 This function then sets \f$ m_i = t \cdot m_i\f$.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::rowMultiply(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const T t)
{
	unsigned int n = m->size1();
	assert(i < n && i>=0);

	for(unsigned int k=i; k<n; k++)
		(*m)(i,k) *= t;
}


// base change for gram matrices
/*!
 *   @brief Swaps two basis vectors.
 *
 *	 @section Description
 *
 *   @param g Constant pointer to a symmetric matrix of type T.
 *   @param i Constant unsigned int specifying the i-th basis vector of the underlying module: \f$0 \leq i \leq n\f$.
 *   @param j Constant unsigned int specifying the j-th basis vector of the underlying module: \f$0 \leq j \leq n\f$.
 *
 *   Swaps the i-th and j-th basis vector, where \f$n\f$ is the dimension of the underlying module.<br>
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::baseSwap(boost::numeric::ublas::symmetric_matrix<T>* const g, unsigned int i, unsigned int j)
{
	unsigned int dim = g->size1();
	assert(i<dim && j<dim && i>=0 && j>=0);

	if (i > j)
		std::swap(i, j);
	else if (i==j)
		return;


	// swap diagonal entries
	std::swap((*g)(i,i), (*g)(j,j));

	// run through columns and swap elements
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::symmetric_matrix<T> > columnI(*g, i);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::symmetric_matrix<T> > columnJ(*g, j);

	// save original column
	std::vector<T> tempColI(dim);
	for(unsigned int k=0; k<dim; k++)
		tempColI[k]=columnI[k];

	// change elements
	for(unsigned int k=0; k<dim; k++)
		if(k!=i && k!=j)
		{
			columnI[k]=columnJ[k];
			columnJ[k]=tempColI[k];
		}
}

/*!
 *   @brief Adds a multiple of a basis vector to another basis vector.
 *
 *	 @section Description
 *
 *   @param g Constant pointer to a symmetric matrix of type T.
 *   @param i Constant unsigned int specifying the i-th basis vector of the underlying module: \f$0 \leq i \leq n\f$.
 *   @param j Constant unsigned int specifying the j-th basis vector of the underlying module: \f$0 \leq j \leq n\f$.
 *   @param t Constant element of type T.
 *
 *	 Let \f$n\f$ denote the dimension of the underlying module and \f$b_i\f$, \f$b_j\f$ the i-th, resp. j-th basis vector.<br>
 *	 This function then sets \f$ b_i = b_i + t \cdot b_j\f$.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::baseAddMultiple(boost::numeric::ublas::symmetric_matrix<T>* const g, const unsigned int i, const unsigned int j, T t)
{	assert(i < g->size1() && j < g->size1() && i>=0 && j>=0);

// temporarily save needed entries
T tempII = (*g)(i,i);
T tempIJ = (*g)(i,j);
T tempJJ = (*g)(j,j);

// change the i-th row and column
boost::numeric::ublas::matrix_row<boost::numeric::ublas::symmetric_matrix<T> > rowI(*g, i);
boost::numeric::ublas::matrix_row<boost::numeric::ublas::symmetric_matrix<T> > rowJ(*g, j);
boost::numeric::ublas::matrix_column<boost::numeric::ublas::symmetric_matrix<T> > colI(*g, i);
for(unsigned int k=0; k<g->size1(); k++)
	rowI(k) += rowJ(k) * t;

// change the diagonal entry
(*g)(i,i) = tempII + 2*t*tempIJ + t*t*tempJJ;
}

// row changes
/*!
 *   @brief Adds a multiple of a row to another row.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a matrix of type T.
 *   @param i Constant unsigned int specifying a row of the given matrix: \f$0 \leq i \leq n\f$.
 *   @param j Constant unsigned int specifying a row of the given matrix: \f$0 \leq j \leq n\f$.
 *   @param t Constant element of type T.
 *
 *	 Let \f$n\f$ denote the number of rows of the given matrix and \f$m_i\f$, \f$m_j\f$ the i-th, resp. j-th row of the matrix m.<br>
 *	 This function then sets \f$ m_i = m_i + t \cdot m_j\f$.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::rowAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, T t)
{
	assert(i < m->size1() && j < m->size1() && i>=0 && j>=0);

	// get columns
	boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<T> > rowI((*m), i);
	boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<T> > rowJ((*m), j);

	// add multiple
	for(unsigned int k=0; k<rowI.size(); k++)
		rowI[k] += t*rowJ[k];
}

// general utility functions for base and gram matrices
/*!
 *   @brief Deletes columns which are filled with only 0s.
 *
 *	 @section Description
 *
 *   @param b Constant pointer to a matrix of type T.
 *   @param rk Constant unsigned int specifying the rank of the matrix: \f$0 \leq i \leq n\f$.
 *
 *   Deleted the first n-rk columns of the given matrix m, where n is the number of columns of b.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::baseDelZeroes(boost::numeric::ublas::matrix<T>* const b, const unsigned int rk)
{
	// delete zero columns - expects the first n-rk columns to be zero!
	assert(rk <= b->size2() && rk>=0);

	if(rk == b->size2())
		return;

	// copy b
	unsigned int n = b->size1(), m = b->size2();
	boost::numeric::ublas::matrix<T> copyB((*b));

	// copy back non-zero columns
	b->resize(n, rk, false);
	for(unsigned int i=0; i<rk; i++)
		for(unsigned int j=0; j<n; j++)
			(*b)(j,i) = copyB(j,i+m-rk);
}

/*!
 *   @brief Sets a specific column of a matrix to a given column.
 *
 *	 @section Description
 *
 *   @param b Constant pointer to a matrix of type T.
 *   @param col Constant pointer to a constant valarray of type T defining the new column to be set.
 *   @param i Constant unsigned int specifying the column of the matrix to be changed: \f$0 \leq i \leq n\f$.
 *
 *   Sets the i-th column of the given matrix b to col. \f$n is the number of columns of b.\f$
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::setColumn(boost::numeric::ublas::matrix<T>* const b, const std::valarray<T>* const col, const unsigned int i)
{
	// set column i of b to col
	assert(i < b->size2() && i>=0);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > colI(*b, i);
	for(unsigned int j=0; j<b->size1(); j++)
		colI(j) = (*col)[j];
}

// check properties
/*!
 *   @brief Checks for orthogonality.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a constant matrix of type T.
 *   @param g Constant pointer to a constant symmetric matrix of type T.
 *
 *   @return Returns true if and only if \f$m^t \cdot m=id\f$ or \f$m^t \cdot g \cdot m=g\f$.
 */
template<class T>
bool MatrixHelper<T>::isOrthogonal(const boost::numeric::ublas::matrix<T>* const m, const boost::numeric::ublas::symmetric_matrix<T>* const g)
{
	if(g == nullptr)
	{
		boost::numeric::ublas::matrix<T> transpose((*m)), identity(m->size1(), m->size2());
		identity.assign(boost::numeric::ublas::identity_matrix<T> (m->size1()));
		transpose = boost::numeric::ublas::prod(transpose,(*m));

		for(unsigned int i=0; i<m->size1(); i++)
			for(unsigned int j=0; j<m->size2(); j++)
				if(transpose(i,j) != identity(i,j))
					return false;
		return true;
	}
	else
	{
		boost::numeric::ublas::matrix<T> transpose((*m)), res(m->size1(), m->size2());
		transpose = boost::numeric::ublas::trans((*m));
		res = boost::numeric::ublas::prod((*g), (*m));
		res = boost::numeric::ublas::prod(transpose, res);

		for(unsigned int i=0; i<m->size1(); i++)
			for(unsigned int j=0; j<m->size2(); j++)
				if(res(i,j) != (*g)(i,j))
					return false;
		return true;
	}
}

// getters
/*!
 *   @brief Retrieves a column of a matrix.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a constant matrix of type T.
 *   @param i Constant unsigned int specifying a column of the given matrix: \f$0 \leq i \leq n\f$, where \f$n\f$ denotes the number of columns of the given matrix.
 *
 *   @return The i-th column of the given matrix as a std::valarray<T>.
 */
template<class T>
std::valarray<T> MatrixHelper<T>::getColumn(const boost::numeric::ublas::matrix<T>* const m, const unsigned int i)
{
	assert(i<m->size2() && i>=0);

	// get column
	boost::numeric::ublas::matrix_column<const boost::numeric::ublas::matrix<T> > columnI(*m, i);

	// save entries in valarray
	std::valarray<T> col(columnI.size());
	for(unsigned int i=0; i<col.size(); i++)
		col[i] = columnI(i);

	return col;
}

/*!
 *   @brief Retrieves a column of a matrix.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a constant symmetric matrix of type T.
 *   @param i Constant unsigned int specifying a column of the given matrix: \f$0 \leq i \leq n\f$, where \f$n\f$ denotes the number of columns of the given matrix.
 *
 *   @return The i-th column of the given matrix as a std::valarray<T>.
 */
template<class T>
std::valarray<T> MatrixHelper<T>::getColumn(const boost::numeric::ublas::symmetric_matrix<T>* const m, const unsigned int i)
{
	assert(i<m->size2() && i>=0);

	// get column
	boost::numeric::ublas::matrix_column<const boost::numeric::ublas::symmetric_matrix<T> > columnI(*m, i);

	// save entries in valarray
	std::valarray<T> col(columnI.size());
	for(unsigned int i=0; i<col.size(); i++)
		col[i] = columnI(i);

	return col;
}

// compare
/*!
 *   @brief Checks for equality.
 *
 *	 @section Description
 *
 *   @param m1 Constant pointer to a constant matrix of type T.
 *   @param m2 Constant pointer to a constant matrix of type T.
 *
 *   @return Returns true if and only if both matrices are equal.
 */
template<class T>
bool MatrixHelper<T>::areEqual(const boost::numeric::ublas::matrix<T>* const m1, const boost::numeric::ublas::matrix<T>* const m2)
{
	bool res = true;
	assert(m1->size1() == m2->size1() && m1->size2() == m2->size2());

	for(unsigned int i=0; i<m1->size1(); i++)
		for(unsigned int j=0; j<m1->size2(); j++)
			if((*m1)(i,j) != (*m2)(i,j))
				res = false;

	return res;
}

// input and output
/*!
 *   @brief Prints a matrix in the old TN format.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a constant matrix of type T.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::printMatrix(const boost::numeric::ublas::matrix<T>* const m)
{
	for(typename boost::numeric::ublas::matrix<T>::const_iterator1 it1 = m->begin1(); it1 != m->end1(); it1++)
	{
		for(typename boost::numeric::ublas::matrix<T>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); it2++)
			std::cout << *it2 << "\t";
		std::cout << std::endl;
	}
}

/*!
 *   @brief Prints a matrix in the old TN format.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a constant symmetric matrix of type T.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::printMatrix(const boost::numeric::ublas::symmetric_matrix<T>* const m)
{
	for(typename boost::numeric::ublas::symmetric_matrix<T>::const_iterator1 it1 = m->begin1(); it1 != m->end1(); it1++)
	{
		for(typename boost::numeric::ublas::symmetric_matrix<T>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); it2++)
			std::cout << *it2 << "\t";
		std::cout << std::endl;
	}
}

/*!
 *   @brief Prints a matrix in the old TN format.
 *
 *	 @section Description
 *
 *   @param m Constant pointer to a constant upper triangular matrix of type T.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::printMatrix(const boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m)
{
	unsigned int n = m->size1();
	for(typename boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>::const_iterator1 it1 = m->begin1(); it1 != m->end1(); it1++)
	{
		for(typename boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>::const_iterator2 it2 = it1.begin(); it2 != it1.end(); it2++)
			std::cout << *it2 << "\t";
		std::cout << std::endl;
	}
}

/*!
 *   @brief Prints a vector.
 *
 *	 @section Description
 *
 *   @param d Constant pointer to a constant valarray of type T.
 *
 *   @return void
 */
template<class T>
void MatrixHelper<T>::printDiagonalVector(const std::valarray<T>* const d)
{
	for(unsigned int i=0; i<d->size(); i++)
		std::cout << (*d)[i] << " ";
	std::cout << std::endl;
}

}
