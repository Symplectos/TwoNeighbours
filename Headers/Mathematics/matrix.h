#pragma once
/*
 * File:   matrix.h
 * Author: Gilles Bellot
 * Date:   13/10/2016 - Dortmund - Germany
 *
 * Description: wrapper for BOOST uBLAS matrices
 *
 * History:	- 13/10/2016: added matrix reading (boost) and printing (human)
 * 			- 14/10/2016: added support for unit matrix
 * 			- 15/10/2016: added base change (add multiple) for matrices
 * 			- 17/10/2016: added base change (add multiple) for gramians
 * 			- 17/10/2016: added base change (swap) for matrices
 * 			- 02/11/2016: added base change (swap) for gramians
 * 			- 17/11/2016: added base change (add multiple, multiply) for triangular matrices
 * 			- 17/11/2016: added print support for triangular matrices and vectors
 * 			- 25/11/2016: added matrix multiplication and inversion for triangular matrices
 * 			- 29/11/2016: fixed an error in trig matrix multiplication and inversion
 * 			- 05/01/2017: fixed a template error
 * 			- 05/01/2017: deleted readMatrix
 * 			- 17/01/2017: added function to delete zero columns from basis matrix
 * 			- 18/01/2017: deleted makeUnit function -> replaced by boosts own assign function
 * 			- 18/01/2017: deleted mult function -> replaced by boosts own prod function
 * 			- 18/01/2017: deleted read function -> replaced by boosts own read function on std::cin
 * 			- 18/01/2017: deleted printSymmetrixMatrix function
 *			- 18/01/2017: clean up of printMatrix functions
 *			- 18/01/2017: added det function (with lu)
 *			- 18/01/2017: changed print functions to const
 *			- 18/01/2017: changed inverse function to const
 *			- 19/01/2017: added gauss bareiss for determinant computation
 *			- 19/01/2017: added gauss bareiss for rank computation
 *			- 19/01/2017: added orthogonality check
 *			- 20/01/2017: added orthoginality check with inner product
 *			- 20/01/2017: added comparison function
 *			- 21/01/2017: moved rank, det and inverse function to the linear algebra class
 *			- 01/08/2017: added const declarations
 *
 * ToDo: everything :(
 */

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

// helper functions to use BOOST uBLAS matrices
template<class T>
class MatrixHelper
{
private:

public:
	// elementary column operations
	static void colSwap(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j);				// swap columns i and j
	static void colAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, T t);	// column i = column i + t * column j
	static void colMultiply(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const T t);					// col i = t * col i
	static std::valarray<T> getColumn(const boost::numeric::ublas::matrix<T>* const m, const unsigned int i);				// returns the i-th column of m
	static std::valarray<T> getColumn(const boost::numeric::ublas::symmetric_matrix<T>* const m, const unsigned int i);		// returns the i-th column of m

	// elementary row operations
	static void rowAddMultiple(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const unsigned int j, T t);	// row i = row i + T * row j
	static void rowMultiply(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const T t);						// row i *= t
	static void rowAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, const T t);										// row i = row i + t * row j

	// base change for gram matrices
	static void baseSwap(boost::numeric::ublas::symmetric_matrix<T>* const g, unsigned int i, unsigned int j);							// swap base vectors i and j
	static void baseAddMultiple(boost::numeric::ublas::symmetric_matrix<T>* const g, const unsigned int i, const unsigned int j, T t);	// b_i = b_i + q * b_j

	// general utility functions for base or gram matrices
	static void baseDelZeroes(boost::numeric::ublas::matrix<T>* const b, const unsigned int rk);								// delete zero columns of m
	static void setColumn(boost::numeric::ublas::matrix<T>* const b, const std::valarray<T>* const col, const unsigned int i);	// sets column i of m to col

	// check properties
	static bool isOrthogonal(const boost::numeric::ublas::matrix<T>* const m, const boost::numeric::ublas::symmetric_matrix<T>* const g = nullptr);	// returns true iff m^t*m=id or m^t*g*m=g

	// compare
	static bool areEqual(const boost::numeric::ublas::matrix<T>* const m1, const boost::numeric::ublas::matrix<T>* const m2);	// returns true iff m1 == m2

	// input and output
	static void printMatrix(const boost::numeric::ublas::matrix<T>* const m);											// print matrix in tn format
	static void printMatrix(const boost::numeric::ublas::symmetric_matrix<T>* const m);									// print symmetric matrix in tn format
	static void printMatrix(const boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m); 	// print upper triangular matrx in tn format
	static void printDiagonalVector(const std::valarray<T>* const d);													// print diagonal vector
};


// IMPLEMENTATION ///////////////////////////////////////////////////////////////////////

// base change for matrices
template<class T>
void MatrixHelper<T>::colSwap(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j)
{
	assert(i<m->size2() && j<m->size2());

	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnI(*m, i);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnJ(*m, j);

	boost::swap(columnI, columnJ);
}

template<class T>
void MatrixHelper<T>::colAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, const T t)
{
	assert(i < m->size2() && j < m->size2());

	// get the columns in question
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnI(*m, i);
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnJ(*m, j);

	// add multiple
	for(unsigned int k=0; k<columnI.size(); k++)
		columnI(k)+=t*columnJ(k);
}

template<class T>
void MatrixHelper<T>::colMultiply(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const T t)
{
	assert(i < m->size2());

	// get the column in question
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > columnI(*m, i);

	// multiply
	for(unsigned int k=0; k<columnI.size(); k++)
		columnI(k)*=t;
}

// base change for triangular matrices
template<class T>
void MatrixHelper<T>::rowAddMultiple(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const unsigned int j, const T t)
{
	assert(i < m->size1() && j < m->size1());
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

template<class T>
void MatrixHelper<T>::rowMultiply(boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>* const m, const unsigned int i, const T t)
{
	unsigned int n = m->size1();
	assert(i < n);

	for(unsigned int k=i; k<n; k++)
		(*m)(i,k) *= t;
}


// base change for gram matrices
template<class T>
void MatrixHelper<T>::baseSwap(boost::numeric::ublas::symmetric_matrix<T>* const g, unsigned int i, unsigned int j)
{
	unsigned int dim = g->size1();
	assert(i<dim && j<dim);

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

template<class T>
void MatrixHelper<T>::baseAddMultiple(boost::numeric::ublas::symmetric_matrix<T>* const g, const unsigned int i, const unsigned int j, T t)
{	assert(i < g->size1() && j < g->size1());

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
template<class T>
void MatrixHelper<T>::rowAddMultiple(boost::numeric::ublas::matrix<T>* const m, const unsigned int i, const unsigned int j, T t)
{
	assert(i < m->size1() && j < m->size1());

	// get columns
	boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<T> > rowI((*m), i);
	boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<T> > rowJ((*m), j);

	// add multiple
	for(unsigned int k=0; k<rowI.size(); k++)
		rowI[k] += t*rowJ[k];
}

// general utility functions for base and gram matrices
template<class T>
void MatrixHelper<T>::baseDelZeroes(boost::numeric::ublas::matrix<T>* const b, const unsigned int rk)
{
	// delete zero columns - expects the first n-rk columns to be zero!
	assert(rk <= b->size2());

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

template<class T>
void MatrixHelper<T>::setColumn(boost::numeric::ublas::matrix<T>* const b, const std::valarray<T>* const col, const unsigned int i)
{
	// set column i of b to col
	assert(i < b->size2());
	boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<T> > colI(*b, i);
	for(unsigned int j=0; j<b->size1(); j++)
		colI(j) = (*col)[j];
}

// check properties
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
template<class T>
std::valarray<T> MatrixHelper<T>::getColumn(const boost::numeric::ublas::matrix<T>* const m, const unsigned int i)
{
	assert(i<m->size2());

	// get column
	boost::numeric::ublas::matrix_column<const boost::numeric::ublas::matrix<T> > columnI(*m, i);

	// save entries in valarray
	std::valarray<T> col(columnI.size());
	for(unsigned int i=0; i<col.size(); i++)
		col[i] = columnI(i);

	return col;
}

template<class T>
std::valarray<T> MatrixHelper<T>::getColumn(const boost::numeric::ublas::symmetric_matrix<T>* const m, const unsigned int i)
{
	assert(i<m->size2());

	// get column
	boost::numeric::ublas::matrix_column<const boost::numeric::ublas::symmetric_matrix<T> > columnI(*m, i);

	// save entries in valarray
	std::valarray<T> col(columnI.size());
	for(unsigned int i=0; i<col.size(); i++)
		col[i] = columnI(i);

	return col;
}


// compare
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

template<class T>
void MatrixHelper<T>::printDiagonalVector(const std::valarray<T>* const d)
{
	for(unsigned int i=0; i<d->size(); i++)
		std::cout << (*d)[i] << " ";
	std::cout << std::endl;
}

}
