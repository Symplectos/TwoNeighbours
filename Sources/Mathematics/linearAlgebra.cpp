/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	linearAlgebra.cpp
 * @date	14/10/2016 - Eichlinghofen - Germany
 *
 * @brief	Implementation for general linear algebra and lattice algorithms.
 *
 * @section Bibliography
 * - [CH] Cohen, H. --- A Course in Computational Algebraic Number Theory
 * - [HB] Hemkemeier, B. --- Algorithmische Konstruktion von Gittern
  *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/
#include "../../Headers/Mathematics/linearAlgebra.h"
#include "../../Headers/Mathematics/matrix.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <numeric>

#include <gmp.h>
#include <gmpxx.h>

namespace mathematics
{



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Comparison //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Lexicographically compares two pairs consisting of a multi-precision integer and a valarray of the same type.
 *
 *   @param a Reference to a constant pair of a mpz_class integer and a mpz_class valarray.
 *   @param b Reference to a constant pair of a mpz_class integer and a mpz_class valarray.
 *
 *   @return True if and only if the a is lexicographically smaller than b.
 */
bool LALA::comparePairMPZValarray(const std::pair<mpz_class, std::valarray<mpz_class> >& a, const std::pair<mpz_class, std::valarray<mpz_class> >& b)
{
	if(a.first < b.first)
		return true;
	else if(a.first > b.first)
		return false;
	else
	{
		mpz_class na, nb;
		for(unsigned int i=0; i<a.second.size(); i++)
			na += a.second[i]*a.second[i];
		for(unsigned int i=0; i<b.second.size(); i++)
			nb += b.second[i]*b.second[i];
		return na < nb;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Private Scalar Products//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the scalar product of two basis vectors.
 *
 *   @param b Constant pointer to a constant matrix of mpreals whose columns represent a basis for the underlying module.
 *   @param i Constant integer to specify the i-th basis vector.
 *   @param j Constant integer to specify the j-th basis vector.
 *   @param innerProduct Constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return The value of the scalar product of the two specified basis vectors.
 *
 *   In the case where no inner product is specified, the loops are unrolled manually for a significant speed burst.
 *   @todo Unroll the loops in the presence of an inner product as well. Register blocking or prefetching will probably lead to a great increase in speed.
 */
mpfr::mpreal LALA::lllBasisDot(const boost::numeric::ublas::matrix<mpfr::mpreal>* const b, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	assert(i < b->size2() && j < b->size2());

	mpfr::mpreal result = 0;

	if(innerProduct == nullptr)
	{
		unsigned int n = b->size1();
		if(n%4==0)
		{
			mpfr::mpreal dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0;
			for(unsigned int k = 0; k < n; k += 4)
			{
				dot1 += (*b)(k,i) * (*b)(k,j);
				dot2 += (*b)(k + 1,i) * (*b)(k + 1,j);
				dot3 += (*b)(k + 2,i) * (*b)(k + 2,j);
				dot4 += (*b)(k + 3,i) * (*b)(k + 3,j);
			}
			result += dot1 + dot2 + dot3 + dot4;
		}
		else if(n%2==0)
		{
			mpfr::mpreal dot1 = 0, dot2 = 0;
			for(unsigned int k = 0; k < n; k += 2)
			{
				dot1 += (*b)(k,i) * (*b)(k,j);
				dot2 += (*b)(k + 1,i) * (*b)(k + 1,j);
			}
			result += dot1 + dot2;
		}
		else
		{
			mpfr::mpreal dot1 = 0, dot2 = 0;
			unsigned int nn = n/2;
			for(unsigned int k = 0; k < n && nn>0; k += 2, nn--)
			{
				dot1 += (*b)(k,i) * (*b)(k,j);
				dot2 += (*b)(k + 1,i) * (*b)(k + 1,j);
			}
			result += dot1 + dot2 + (*b)(n-1,i) * (*b)(n-1,j);;
		}
	}
	else
	{
		// copy inner product into mpfr matrix
		boost::numeric::ublas::symmetric_matrix<mpfr::mpreal> copyInnerProduct(innerProduct->size1(), innerProduct->size2());
		for(unsigned int a=0; a<innerProduct->size1(); a++)
			for(unsigned int b=0; b<innerProduct->size2(); b++)
				mpfr_set_z(copyInnerProduct(a,b).mpfr_ptr(), (*innerProduct)(a,b).get_mpz_t(), MPFR_RNDN);

		// get columns
		boost::numeric::ublas::vector<mpfr::mpreal> colI(b->size1()), colJ(b->size1());
		for(unsigned int k=0; k<b->size1(); k++)
		{
			colI(k) = (*b)(k,i);
			colJ(k) = (*b)(k,j);
		}

		colJ = boost::numeric::ublas::prod(copyInnerProduct, colJ);
		result = boost::numeric::ublas::inner_prod(colI, colJ);
	}

	return result;
}

/*!
 *   @brief Computes the scalar product of two basis vectors.
 *
 *   @param b Constant pointer to a constant matrix of mpreals whose columns represent a basis for the underlying module.
 *   @param bStar Constant pointer to a constant matrix of mpreals whose columns represent a basis for the underlying module.
 *   @param i Constant integer to specify the i-th basis vector in b.
 *   @param j Constant integer to specify the j-th basis vector in bStar.
 *   @param innerProduct Constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return The value of the scalar product of the two specified basis vectors.
 *
 *   In the case where no inner product is specified, the loops are unrolled manually for a significant speed burst.
 *   @todo Unroll the loops in the presence of an inner product as well. Register blocking or prefetching will probably lead to a great increase in speed.
 */
mpfr::mpreal LALA::lllBasisDot(const boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	assert(i < b->size2() && j < bStar->size2());

	mpfr::mpreal result = 0;

	if(innerProduct == nullptr)
	{
		unsigned int n = b->size1();
		if(n%4==0)
		{
			mpfr::mpreal dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0;
			for(unsigned int k = 0; k < n; k += 4)
			{
				dot1 += (*b)(k,i) * (*bStar)(k,j);
				dot2 += (*b)(k + 1,i) * (*bStar)(k + 1,j);
				dot3 += (*b)(k + 2,i) * (*bStar)(k + 2,j);
				dot4 += (*b)(k + 3,i) * (*bStar)(k + 3,j);
			}
			result += dot1 + dot2 + dot3 + dot4;
		}
		else if(n%2==0)
		{
			mpfr::mpreal dot1 = 0, dot2 = 0;
			for(unsigned int k = 0; k < n; k += 2)
			{
				dot1 += (*b)(k,i) * (*bStar)(k,j);
				dot2 += (*b)(k + 1,i) * (*bStar)(k + 1,j);
			}
			result += dot1 + dot2;
		}
		else
		{
			mpfr::mpreal dot1 = 0, dot2 = 0;
			unsigned int nn = n/2;
			for(unsigned int k = 0; k < n && nn>0; k += 2, nn--)
			{
				dot1 += (*b)(k,i) * (*bStar)(k,j);
				dot2 += (*b)(k + 1,i) * (*bStar)(k + 1,j);
			}
			result += dot1 + dot2 + (*b)(n-1,i) * (*bStar)(n-1,j);;
		}
	}
	else
	{
		// copy inner product into mpfr matrix
		boost::numeric::ublas::symmetric_matrix<mpfr::mpreal> copyInnerProduct(innerProduct->size1(), innerProduct->size1());
		for(unsigned int a=0; a<innerProduct->size1(); a++)
			for(unsigned int b=0; b<innerProduct->size2(); b++)
				mpfr_set_z(copyInnerProduct(a,b).mpfr_ptr(), (*innerProduct)(a,b).get_mpz_t(), MPFR_RNDN);

		// get columns
		boost::numeric::ublas::vector<mpfr::mpreal> colI(b->size1()), colJ(b->size1());
		for(unsigned int k=0; k<b->size1(); k++)
		{
			colI(k) = (*b)(k,i);
			colJ(k) = (*bStar)(k,j);
		}
		colJ = boost::numeric::ublas::prod(copyInnerProduct, colJ);
		result = boost::numeric::ublas::inner_prod(colI, colJ);
	}

	return result;
}

/*!
 *   @brief Computes the scalar product of two basis vectors.
 *
 *   @param b Constant pointer to a constant matrix of multi-precision integers whose columns represent a basis for the underlying module.
 *   @param i Constant integer to specify the i-th basis vector.
 *   @param j Constant integer to specify the j-th basis vector.
 *   @param innerProduct Constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return The value of the scalar product of the two specified basis vectors.
 *
 *   In the case where no inner product is specified, the loops are unrolled manually for a significant speed burst.
 *   @todo Unroll the loops in the presence of an inner product as well. Register blocking or prefetching will probably lead to a great increase in speed.
 */
mpz_class LALA::lllBasisDot(const boost::numeric::ublas::matrix<mpz_class>* const b, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	assert(i < b->size2() && j < b->size2());

	mpz_class result = 0;

	if(innerProduct == nullptr)
	{
		unsigned int n = b->size1();
		if(n%4==0)
		{
			mpz_class dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0;
			for(unsigned int k = 0; k < n; k += 4)
			{
				dot1 += (*b)(k,i) * (*b)(k,j);
				dot2 += (*b)(k + 1,i) * (*b)(k + 1,j);
				dot3 += (*b)(k + 2,i) * (*b)(k + 2,j);
				dot4 += (*b)(k + 3,i) * (*b)(k + 3,j);
			}
			result += dot1 + dot2 + dot3 + dot4;
		}
		else if(n%2==0)
		{
			mpz_class dot1 = 0, dot2 = 0;
			for(unsigned int k = 0; k < n; k += 2)
			{
				dot1 += (*b)(k,i) * (*b)(k,j);
				dot2 += (*b)(k + 1,i) * (*b)(k + 1,j);
			}
			result += dot1 + dot2;
		}
		else
		{
			mpz_class dot1 = 0, dot2 = 0;
			unsigned int nn = n/2;
			for(unsigned int k = 0; k < n && nn>0; k += 2, nn--)
			{
				dot1 += (*b)(k,i) * (*b)(k,j);
				dot2 += (*b)(k + 1,i) * (*b)(k + 1,j);
			}
			result += dot1 + dot2 + (*b)(n-1,i) * (*b)(n-1,j);;
		}
	}
	else
	{
		boost::numeric::ublas::vector<mpz_class> colI(b->size1()), colJ(b->size1());
		for(unsigned int k=0; k<b->size1(); k++)
		{
			colI(k) = (*b)(k,i);
			colJ(k) = (*b)(k,j);
		}
		colJ = boost::numeric::ublas::prod((*innerProduct), colJ);
		result = boost::numeric::ublas::inner_prod(colI, colJ);
	}
	return result;
}

/*!
 *   @brief Computes the scalar product of two basis vectors.
 *
 *   @param b Constant pointer to a constant matrix of multi-precision integers whose columns represent a basis for the underlying module.
 *   @param bStar Constant pointer to a constant matrix of mpreals whose columns represent a basis for the underlying module.
 *   @param i Constant integer to specify the i-th basis vector in b.
 *   @param j Constant integer to specify the j-th basis vector in bStar.
 *   @param innerProduct Constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return The value of the scalar product of the two specified basis vectors.
 *
 *   In the case where no inner product is specified, the loops are unrolled manually for a significant speed burst.
 *   @todo Unroll the loops in the presence of an inner product as well. Register blocking or prefetching will probably lead to a great increase in speed.
 */
mpfr::mpreal LALA::lllBasisDot(const boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, const unsigned int i, const unsigned int j, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	assert(i < b->size2() && j < bStar->size2());

	mpfr::mpreal result = 0;

	if(innerProduct == nullptr)
	{
		unsigned int n = b->size1();
		if(n%4==0)
		{
			mpfr::mpreal dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0;
			for(unsigned int k = 0; k < n; k += 4)
			{
				dot1 += (*b)(k,i).get_mpz_t() * (*bStar)(k,j);
				dot2 += (*b)(k + 1,i).get_mpz_t() * (*bStar)(k + 1,j);
				dot3 += (*b)(k + 2,i).get_mpz_t() * (*bStar)(k + 2,j);
				dot4 += (*b)(k + 3,i).get_mpz_t() * (*bStar)(k + 3,j);
			}
			result += dot1 + dot2 + dot3 + dot4;
		}
		else if(n%2==0)
		{
			mpfr::mpreal dot1 = 0, dot2 = 0;
			for(unsigned int k = 0; k < n; k += 2)
			{
				dot1 += (*b)(k,i).get_mpz_t() * (*bStar)(k,j);
				dot2 += (*b)(k + 1,i).get_mpz_t() * (*bStar)(k + 1,j);
			}
			result += dot1 + dot2;
		}
		else
		{
			mpfr::mpreal dot1 = 0, dot2 = 0;
			unsigned int nn = n/2;
			for(unsigned int k = 0; k < n && nn>0; k += 2, nn--)
			{
				dot1 += (*b)(k,i).get_mpz_t() * (*bStar)(k,j);
				dot2 += (*b)(k + 1,i).get_mpz_t() * (*bStar)(k + 1,j);
			}
			result += dot1 + dot2 + (*b)(n-1,i).get_mpz_t() * (*bStar)(n-1,j);;
		}
	}
	else
	{
		boost::numeric::ublas::vector<mpfr::mpreal> colI(b->size1()), colJ(b->size1());
		for(unsigned int k=0; k<b->size1(); k++)
		{
			mpfr_set_z(colI[k].mpfr_ptr(), (*b)(k,i).get_mpz_t(), MPFR_RNDN);
			colJ(k) = (*bStar)(k,j);
		}
		// copy matrix into mpfr::mpreal
		boost::numeric::ublas::symmetric_matrix<mpfr::mpreal> copyInnerProduct(innerProduct->size1(), innerProduct->size2());
		for(unsigned int a=0; a<innerProduct->size1(); a++)
			for(unsigned int b=0; b<innerProduct->size2(); b++)
				mpfr_set_z(copyInnerProduct(a,b).mpfr_ptr(), (*innerProduct)(a,b).get_mpz_t(), MPFR_RNDN);

		colJ = boost::numeric::ublas::prod(copyInnerProduct, colJ);
		result = boost::numeric::ublas::inner_prod(colI, colJ);
	}

	return result;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Cholesky ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the rational Cholesky decomposition of a symmetric matrix.
 *
 *   @param[in] q Constant pointer to a constant symmetric matrix of multi-precision integers representing a quadratic form.
 *   @param[out] r Constant pointer to an upper triangular matrix of multi-precision rationals.
 *   @param[out] d Constant pointer to a valarray of multi-precision rationals.
 *
 *	 Given a symmetric positive definite matrix q, this algorithm outputs an upper triangular matrix r and a diagonal matrix d such that \f$q = r^t \cdot d \cdot r\f$.
 *
 *   @return void
 *
 *   @todo Assert that q is a quadratic form.
 *
 */
void LALA::choleskyDecomposition(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const q, boost::numeric::ublas::triangular_matrix<mpq_class, boost::numeric::ublas::upper>* const r, std::valarray<mpq_class>* const d)
{
	// given a symmetric positive definite matrix q (representing a quadratic form), this algorithm outputs an upper triangular matrix r and a diagonal matrix d such that
	// q = r^t*d*r

	// step 0: initialize
	unsigned int n = q->size1();

	if(r->size1()!=n)
		r->resize(n,n,false);
	if(d->size() != n)
		d->resize(n);
	assert(r->size1()==n && d->size() == n);

	// step 1: copy q into r
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=i; j<n; j++)
			(*r)(i,j) = (*q)(i,j);

	// step 2: compute r
	for(unsigned int j=0; j<n-1; j++)
		for(unsigned int i=j+1; i<n; i++)
			MatrixHelper<mpq_class>::rowAddMultiple(r, i, j, -(*r)(j,i)/(*r)(j,j));

	// step 3: get diagonal entries from r
	for(unsigned int i=0; i<n; i++)
		(*d)[i] = (*r)(i,i);

	// step 4: divide row i of r by the inverse of d_i
	for(unsigned int i=0; i<n; i++)
		MatrixHelper<mpq_class>::rowMultiply(r, i, 1/(*d)[i]);
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// FP-LLL //////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Tests the LLL condition on given basis vectors.
 *
 *   @param[in, out] b Constant pointer to a matrix of mpreals whose columns specify a basis of the underlying module.
 *   @param[in, out] bStar Constant pointer to a matrix of mpreals whose columns specify an orthogonal basis.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] B Constant pointer to a valarray of mpreals holding the length of the basis vectors in b.
 *   @param[in, out] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in, out] k A constant pointer to an integer. Specifies the depth of the algorithm, or the current column to modify.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector that was already modified.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 This recursive function tests whether a basis fulfils the LLL condition. To do so, a reduction step is done.
 *	 If after that reduction, the LLL condition is fulfilled, the function returns, else a swap step is done and k is reset before the function recalls itself.
 *
 *   @return void
 */
void LALA::lllBasisTestCondition(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, int* const k, const int kmax, const int dim)
{
	// test LLL condition
	LALA::lllBasisRed(b, h, mu, *k, *k-1, dim);
	if((*B)[*k] < (0.75-(*mu)[*k*dim+*k-1]*(*mu)[*k*dim+*k-1])*(*B)[*k-1])
	{
		// swap
		LALA::lllBasisSwap(b, bStar, h, B, mu, *k, kmax, dim);
		(*k-1 > 1) ? (*k=*k-1) : (*k=1);
		LALA::lllBasisTestCondition(b, bStar, h, B, mu, k, kmax, dim);
	}
}

/*!
 *   @brief Length reduction in the LLL algorithm.
 *
 *   @param[in, out] b Constant pointer to a matrix of mpreals whose columns specify a basis of the underlying module.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] l A constant integer specifying the basis vector to reduce with.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 Performs length reduction if necessary, i.e. sets \f$b_k := b_k - \lfloor \mu_{k,l} \rfloor \cdot b_l\f$.
 *
 *   @return void
 */
void LALA::lllBasisRed(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const mu, const int k, const int l, const int dim)
{
	mpfr::mpreal mukl = (*mu)[k*dim+l];
	mpfr::mpreal aux = mpfr::abs(mukl, MPFR_RNDN);
	if(aux > 0.5+eps)	// if |mu_{k,l}| > 1/2
	{
		// get the nearest integer to mu_{k,l}
		mpfr::mpreal q = mpfr::rint(mukl, MPFR_RNDN);

		// change matrix columns for H: H_k = H_k - q * H_l
		if(h != nullptr)
			MatrixHelper<mpfr::mpreal>::colAddMultiple(h,k,l,-q);

		// change base: b_k = b_k - q * b_l;
		MatrixHelper<mpfr::mpreal>::colAddMultiple(b,k,l,-q);

		// change the gram-schmidt coefficients accordingly
		(*mu)[k*dim+l] -= q;
		for(int i=0; i<=l-1; i++)
			(*mu)[k*dim+i] -= q*(*mu)[l*dim+i];
	}
}

/*!
 *   @brief Swap step in the LLL algorithm.
 *
 *   @param[in, out] b Constant pointer to a matrix of mpreals whose columns specify a basis of the underlying module.
 *   @param[in, out] bStar Constant pointer to a matrix of mpreals whose columns specify an orthogonal basis.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] B Constant pointer to a valarray of mpreals holding the length of the basis vectors in b.
 *   @param[in, out] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector yet modified.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 Swaps two basis vectors and updates all relevant information, such as the Gram-Schmidt coefficients.
 *
 *   @return void
 */
void LALA::lllBasisSwap(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, const int k, const int kmax, const int dim)
{
	// swap h_k and h_{k-1}
	if(h != nullptr)
		MatrixHelper<mpfr::mpreal>::colSwap(h,k,k-1);

	// swap b_k and b_{k-1}
	MatrixHelper<mpfr::mpreal>::colSwap(b, k, k-1);

	// swap gram-schmidt coefficients accordingly
	if(k>1)
		for(int j=0; j<=k-2; j++)
			boost::swap((*mu)[k*dim+j], (*mu)[(k-1)*dim+j]);

	// change gram-schmidt coefficients - get new quotient
	mpfr::mpreal u((*mu)[k*dim+k-1]), BB((*B)[k] + u*u*(*B)[k-1]);
	(*mu)[k*dim+k-1] = u*(*B)[k-1] / BB;

	std::valarray<mpfr::mpreal> bb(dim);
	for(int i=0; i<dim; i++)
	{
		bb[i] = (*bStar)(i,k-1);
		(*bStar)(i, k-1) = (*bStar)(i,k) + u*bb[i];
		(*bStar)(i,k) = -(*mu)[k*dim+k-1]*(*bStar)(i,k) + ((*B)[k]/BB)*bb[i];
	}

	(*B)[k] = (*B)[k-1]*(*B)[k] / BB;
	(*B)[k-1] = BB;

	mpfr::mpreal t;

	// recompute gram-schmidt coefficients
	for(int i=k+1; i<=kmax; i++)
	{
		t = (*mu)[i*dim+k];
		(*mu)[i*dim+k] = (*mu)[i*dim+k-1] - u*t;
		(*mu)[i*dim+k-1] = t + (*mu)[k*dim+k-1]*(*mu)[i*dim+k];
	}
}

/*!
 *   @brief Length reduction in the LLL algorithm.
 *
 *   @param[in, out] g Constant pointer to a symmetric matrix of mpreals specifying a gram matrix.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] l A constant integer specifying the basis vector to reduce with.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 Performs length reduction if necessary, i.e. sets \f$b_k := b_k - \lfloor \mu_{k,l} \rfloor \cdot b_l\f$, where \f$b_k, b_l\f$ denote the k-th, resp. the l-th basis vector.
 *
 *   @return void
 */
void LALA::lllGramRed(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const mu, const int k, const int l, const int dim)
{
	mpfr::mpreal mukl = (*mu)[k*dim+l];
	mpfr::mpreal aux = mpfr::abs(mukl, MPFR_RNDN);
	if(aux > 0.5+eps)	// if |mu_{k,l}| > 1/2
	{
		// get the nearest integer to mu_{k,l}
		mpfr::mpreal q = mpfr::rint(mukl, MPFR_RNDN);

		// change matrix columns for H: H_k = H_k - q * H_l
		if(h != nullptr)
			MatrixHelper<mpfr::mpreal>::colAddMultiple(h,k,l,-q);

		// change gramian: b_k = b_k - q * b_l;
		MatrixHelper<mpfr::mpreal>::baseAddMultiple(g,k,l,-q);

		// change the gram-schmidt coefficients accordingly
		(*mu)[k*dim+l] -= q;
		for(int i=0; i<=l-1; i++)
			(*mu)[k*dim+i] -= q*(*mu)[l*dim+i];
	}
}

/*!
 *   @brief Swap step in the LLL algorithm.
 *
 *   @param[in, out] g Constant pointer to a symmetric matrix of mpreals. This is a gram matrix.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] b Constant pointer to a valarray of mpreals holding the length of the basis vectors.
 *   @param[in, out] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector yet modified.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 Swaps two basis vectors and updates all relevant information, such as the Gram-Schmidt coefficients.
 *
 *   @return void
 */
void LALA::lllGramSwap(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const b, std::valarray<mpfr::mpreal>* const mu, const int k, const int kmax, const int dim)
{
	// swap h_k and h_{k-1}
	if(h != nullptr)
		MatrixHelper<mpfr::mpreal>::colSwap(h,k,k-1);

	// swap b_k and b_{k-1}
	MatrixHelper<mpfr::mpreal>::baseSwap(g, k, k-1);

	// swap gram-schmidt coefficients accordingly
	if(k>1)
		for(int j=0; j<=k-2; j++)
			boost::swap((*mu)[k*dim+j], (*mu)[(k-1)*dim+j]);

	// change gram-schmidt coefficients - get new quotient
	mpfr::mpreal u((*mu)[k*dim+k-1]), B((*b)[k] + u*u*(*b)[k-1]);
	(*mu)[k*dim+k-1] = u*(*b)[k-1] / B;
	(*b)[k] = (*b)[k-1]*(*b)[k] / B;
	(*b)[k-1] = B;

	mpfr::mpreal t;
	// recompute gram-schmidt coefficients
	for(int i=k+1; i<=kmax; i++)
	{
		t = (*mu)[i*dim+k];
		(*mu)[i*dim+k] = (*mu)[i*dim+k-1] - u*t;
		(*mu)[i*dim+k-1] = t + (*mu)[k*dim+k-1]*(*mu)[i*dim+k];
	}
}

/*!
 *   @brief Tests the LLL condition on a given gram matrix.
 *
 *   @param[in] g Constant pointer to a symmetric matrix of mpreals defining a gram matrix.
 *   @param[in] h Constant pointer to a matrix of mpreals. This is the basis change matrix.
 *   @param[in] b Constant pointer to a valarray of mpreals holding the length of the basis vectors.
 *   @param[in] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in, out] k A constant pointer to an integer. Specifies the basis vector to be modified.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector yet seen.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 This recursive function tests whether a gram matrix fulfils the LLL condition. To do so, a reduction step is done.
 *	 If after that reduction, the LLL condition is fulfilled, the function returns, else a swap step is done and k is reset before the function recalls itself.
 *
 *   @return void
 */
void LALA::lllGramTestCondition(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, std::valarray<mpfr::mpreal>* const b, std::valarray<mpfr::mpreal>* const mu, int* const k, const int kmax, const int dim)
{
	// test LLL condition
	LALA::lllGramRed(g,h,mu,*k,*k-1, dim);
	if((*b)[*k] < (0.75-(*mu)[*k*dim+*k-1]*(*mu)[*k*dim+*k-1])*(*b)[*k-1])
	{
		// swap
		LALA::lllGramSwap(g, h, b, mu, *k, kmax, dim);
		(*k-1 > 1) ? (*k=*k-1) : (*k=1);
		LALA::lllGramTestCondition(g,h,b,mu,k,kmax, dim);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Integral LLL ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Tests the LLL condition on a given gram matrix.
 *
 *   @param[in] g Constant pointer to a symmetric matrix of multi-precision integers defining a gram matrix.
 *   @param[in] h Constant pointer to a matrix of multi-precision integers. This is the basis change matrix.
 *   @param[in] d Constant pointer to a valarray of multi-precision integers holding the sub-determinants of the gram matrix.
 *   @param[in] lambda Constant pointer to a valarray of multi-precision integers holding the Gram-Schmidt coefficients.
 *	 @param[in, out] k A constant pointer to an integer. Specifies the depth of the algorithm, or the vector to be modified.
 *	 @param[in] kmax A constant integer specifying the maximal vector that has been modified already.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 This recursive function tests whether a gram matrix fulfils the LLL condition. To do so, a reduction step is done.
 *	 If after that reduction, the LLL condition is fulfilled, the function returns, else a swap step is done and k is reset before the function recalls itself.
 *
 *   @return void
 */
void LALA::lllGramTestCondition(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpz_class>* const d, std::valarray<mpz_class>* const lambda, int* const k, const int kmax, const int dim)
{
	// test LLL condition
	LALA::lllGramRed(g,h,d,lambda,*k,*k-1, dim);
	if(4 * (*d)[*k+1] * (*d)[*k-1] < 3 * (*d)[*k]*(*d)[*k] - 4 * (*lambda)[*k*dim+*k-1]*(*lambda)[*k*dim+*k-1])
	{
		// swap
		LALA::lllGramSwap(g, h, d, lambda, *k, kmax, dim);
		(*k-1 > 1) ? (*k=*k-1) : (*k=1);
		LALA::lllGramTestCondition(g,h,d,lambda,k,kmax, dim);
	}
}

/*!
 *   @brief Length reduction in the LLL algorithm.
 *
 *   @param[in, out] g Constant pointer to a symmetric matrix of mpreals specifying a gram matrix.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in] d Constant pointer to a valarray of multi-precision integers holding the sub-determinants of the gram matrix.
 *   @param[in, out] lambda Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] l A constant integer specifying the basis vector to reduce with.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 Performs length reduction if necessary, i.e. sets \f$b_k := b_k - \lfloor \mu_{k,l} \rfloor \cdot b_l\f$, where \f$b_k, b_l\f$ denote the k-th, resp. the l-th basis vector.
 *
 *   @return void
 */
void LALA::lllGramRed(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpz_class>* const d, std::valarray<mpz_class>* const lambda, const int k, const int l, const int dim)
{
	mpz_class lkl = (*lambda)[k*dim+l], dl = (*d)[l+1];
	mpz_class aux(2*lkl);
	mpz_abs(aux.get_mpz_t(), aux.get_mpz_t());
	if(aux > dl)	// if |2*lambda_{k,l}| > d_l
	{
		// get the quotient of the euclidean division of 2 * \lambda_{k,l}+d_l by 2*d_l
		mpz_class q(0), num(2*lkl+dl), den(2*dl);
		mpz_fdiv_q(q.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());

		// change matrix columns for H: H_k = H_k - q * H_l
		if(h != nullptr)
			MatrixHelper<mpz_class>::colAddMultiple(h,k,l,-q);

		// change gramian: b_k = b_k - q * b_l;
		MatrixHelper<mpz_class>::baseAddMultiple(g,k,l,-q);

		// change the gram-schmidt coefficients accordingly
		(*lambda)[k*dim+l] -= q * dl;
		for(int i=0; i<=l-1; i++)
			(*lambda)[k*dim+i] -= q*(*lambda)[l*dim+i];
	}
}

/*!
 *   @brief Swap step in the LLL algorithm.
 *
 *   @param[in, out] g Constant pointer to a symmetric matrix of multi-precision integers. This is a gram matrix.
 *   @param[in, out] h Constant pointer to a matrix of multi-precision integers. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] d Constant pointer to a valarray of multi-precision integers holding the sub-determinants of the gram matrix.
 *   @param[in, out] lambda Constant pointer to a valarray of multi-precision integers holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector yet modified.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *
 *	 Swaps two basis vectors and updates all relevant information, such as the Gram-Schmidt coefficients.
 *
 *   @return void
 */
void LALA::lllGramSwap(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpz_class>* const d, std::valarray<mpz_class>* const lambda, const int k, const int kmax, const int dim)
{
	// swap h_k and h_{k-1}
	if(h != nullptr)
		MatrixHelper<mpz_class>::colSwap(h,k,k-1);

	// swap b_k and b_{k-1}
	MatrixHelper<mpz_class>::baseSwap(g, k, k-1);

	// swap gram-schmidt coefficients accordingly
	if(k>1)
		for(int j=0; j<=k-2; j++)
			boost::swap((*lambda)[k*dim+j], (*lambda)[(k-1)*dim+j]);

	// change gram-schmidt coefficients - get new quotient
	mpz_class l = (*lambda)[k*dim+k-1];
	mpz_class num(((*d)[k-1]*(*d)[k+1] + l*l)), den((*d)[k]);
	mpz_class b,q,t;
	mpz_tdiv_q(b.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());

	// recompute gram-schmidt coefficients
	for(int i=k+1; i<=kmax; i++)
	{
		t = (*lambda)[i*dim+k];
		num = ((*d)[k+1] * (*lambda)[i*dim+k-1]-l*t), den = (*d)[k];
		mpz_tdiv_q(q.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());
		(*lambda)[i*dim+k] = q;

		num = (b*t+l*(*lambda)[i*dim+k]), den = (*d)[k+1];
		mpz_tdiv_q(q.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());
		(*lambda)[i*dim+k-1] = q;
	}
	// set new sub-determinant
	(*d)[k] = b;
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Integral M-LLL /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Swap step in the LLL algorithm.
 *
 *   @param[in, out] b Constant pointer to a symmetric matrix of multi-precision integers specifying a vector system of the underlying module.
 *   @param[in, out] bStar Constant pointer to a matrix of mpreals whose columns specify an orthogonal basis.
 *   @param[in, out] h Constant pointer to a matrix of multi-precision integers. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in, out] B Constant pointer to a valarray of mpreals holding the length of the vectors in b.
 *   @param[in, out] mu mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector yet modified.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *	 @param[in] nVectors A constant integer specifying the dimension of the underlying module.
 *
 *	 Swaps two vectors and updates all relevant information, such as the Gram-Schmidt coefficients.
 *
 *   @return void
 */
void LALA::mlllBasisSwap(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, const int k, const int kmax, const int dim, const int nVectors)
{
	// swap h_k and h_{k-1}
	if(h != nullptr)
		MatrixHelper<mpz_class>::colSwap(h,k,k-1);

	// swap b_k and b_{k-1}
	MatrixHelper<mpz_class>::colSwap(b, k, k-1);

	// swap gram-schmidt coefficients accordingly
	if(k>1)
		for(int j=0; j<=k-2; j++)
			boost::swap((*mu)[k*nVectors+j], (*mu)[(k-1)*nVectors+j]);

	// change gram-schmidt coefficients - get new quotient
	mpfr::mpreal u((*mu)[k*nVectors+k-1]), BB((*B)[k] + u*u*(*B)[k-1]);

	// if BB = 0, i.e. B_k = 0 and u = 0
	if(BB == 0)
	{
		boost::swap((*B)[k], (*B)[k-1]);
		MatrixHelper<mpfr::mpreal>::colSwap(bStar, k, k-1);
		for(int i=k+1; i<=kmax; i++)
			boost::swap((*mu)[i*nVectors+k], (*mu)[i*nVectors+k-1]);
	}
	else if(u!=0 && (*B)[k] == 0)
	{
		(*B)[k-1] = BB;
		MatrixHelper<mpfr::mpreal>::colMultiply(bStar, k-1, u);
		(*mu)[k*nVectors+k-1] = 1/u;
		for(int i=k+1; i<=kmax; i++)
			(*mu)[i*nVectors+k-1] /= u;
	}
	else if((*B)[k] != 0)
	{
		mpfr::mpreal t = (*B)[k-1] / BB;
		(*mu)[k*nVectors+k-1] = u*t;

		std::valarray<mpfr::mpreal> bb(dim);
		for(int i=0; i<dim; i++)
		{
			bb[i] = (*bStar)(i,k-1);
			(*bStar)(i, k-1) = (*bStar)(i,k) + u*bb[i];
			(*bStar)(i,k) = -(*mu)[k*nVectors+k-1]*(*bStar)(i,k) + ((*B)[k]/BB)*bb[i];
		}

		(*B)[k] *= t;
		(*B)[k-1] = BB;

		// recompute gram-schmidt coefficients
		for(int i=k+1; i<=kmax; i++)
		{
			t = (*mu)[i*nVectors+k];
			(*mu)[i*nVectors+k] = (*mu)[i*nVectors+k-1] - u*t;
			(*mu)[i*nVectors+k-1] = t + (*mu)[k*nVectors+k-1]*(*mu)[i*nVectors+k];
		}
	}
}

/*!
 *   @brief Length reduction in the mLLL algorithm.
 *
 *   @param[in, out] b Constant pointer to a matrix of mpreals whose columns specify a vector system of the underlying module.
 *   @param[in, out] h Constant pointer to a matrix of mpreals. This is the basis change matrix from the original vector system to an LLL-reduced basis.
 *   @param[in, out] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in] k A constant pointer to an integer specifying the basis vector to be reduced.
 *	 @param[in] l A constant integer specifying the basis vector to reduce with.
 *	 @param[in] nVectors A constant integer specifying the dimension of the underlying module.
 *
 *	 Performs length reduction if necessary, i.e. sets \f$b_k := b_k - \lfloor \mu_{k,l} \rfloor \cdot b_l\f$.
 *
 *   @return void
 */
void LALA::mlllBasisRed(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpfr::mpreal>* const mu, const int k, const int l, const int nVectors)
{
	mpfr::mpreal mukl = (*mu)[k*nVectors+l];
	mpfr::mpreal aux = mpfr::abs(mukl, MPFR_RNDN);
	if(aux > 0.5+eps)	// if |mu_{k,l}| > 1/2
	{
		// get the nearest integer to mu_{k,l}
		mpfr::mpreal q = mpfr::rint(mukl, MPFR_RNDN);

		// change matrix columns for H: H_k = H_k - q * H_l
		if(h != nullptr)
			MatrixHelper<mpz_class>::colAddMultiple(h,k,l,-q.toLong());

		// change base: b_k = b_k - q * b_l;
		MatrixHelper<mpz_class>::colAddMultiple(b,k,l,-q.toLong());

		// change the gram-schmidt coefficients accordingly
		(*mu)[k*nVectors+l] -= q;
		for(int i=0; i<=l-1; i++)
			(*mu)[k*nVectors+i] -= q*(*mu)[l*nVectors+i];
	}
}

/*!
 *   @brief Tests the LLL condition on given basis vectors.
 *
 *   @param[in] b Constant pointer to a matrix of multi-precision integers whose columns specify a basis of the underlying module.
 *   @param[in] bStar Constant pointer to a matrix of mpreals whose columns specify an orthogonal basis.
 *   @param[in] h Constant pointer to a matrix of multi-precision integers. This is the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in] B Constant pointer to a valarray of mpreals holding the length of the basis vectors in b.
 *   @param[in] mu Constant pointer to a valarray of mpreals holding the Gram-Schmidt coefficients.
 *	 @param[in, out] k A constant pointer to an integer. Specifies the current basis vector to be modified.
 *	 @param[in] kmax A constant integer specifying the "deepest" basis vector yet modified.
 *	 @param[in] dim A constant integer specifying the dimension of the underlying module.
 *	 @param[in] nVectors A constant integer specifying the number of vectors in the basis.
 *
 *	 This recursive function tests whether a basis fulfils the LLL condition. To do so, a reduction step is done.
 *	 If after that reduction, the LLL condition is fulfilled, the function returns, else a swap step is done and k is reset before the function recalls itself.
 *
 *   @return void
 */
void LALA::mlllBasisTestCondition(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const bStar, boost::numeric::ublas::matrix<mpz_class>* const h, std::valarray<mpfr::mpreal>* const B, std::valarray<mpfr::mpreal>* const mu, int* const k, const int kmax, const int dim, const int nVectors)
{
	// test LLL condition
	LALA::mlllBasisRed(b, h, mu, *k, *k-1, nVectors);
	if((*B)[*k] < (0.75-(*mu)[*k*nVectors+*k-1]*(*mu)[*k*nVectors+*k-1])*(*B)[*k-1])
	{
		// swap
		LALA::mlllBasisSwap(b, bStar, h, B, mu, *k, kmax, dim, nVectors);
		(*k-1 > 1) ? (*k=*k-1) : (*k=1);
		LALA::mlllBasisTestCondition(b, bStar, h, B, mu, k, kmax, dim, nVectors);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Short Vectors //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Fincke-Pohst preprocessing for the short vectors algorithm.
 *
 *   @param[in] g A constant pointer to a constant symmetric matrix of multi-precision integers. This is the quadratic form.
 *   @param[out] h A constant pointer to a matrix of mpreals; holds the reduced basis.
 *   @param[out] p A constant pointer to a matrix of mpreals; the permutation matrix.
 *   @param[out] q A constant pointer to a matrix of mpreals; holds the modified quadratic form.

 *
 *	Reordering of the vectors in the computation of the rational Cholesky decomposition of the Gram matrix will decrease the likelyhood that a partial short vector can not be extended.
 *  In addition, the LLL algorithm is used to modify the quadratic form obtained from the Gram matrix of the lattice basis, which will diminish the ranges for the components of the partial short vectors.
 *  <br>Based on [CH] Algorithm 2.7.7
 *
 *   @todo
 *   - Change algorithm to compute h from the LLL algorithm without inversion.
 *
 *   @return void
 */
void LALA::shortVectorsFinckePohst(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h, boost::numeric::ublas::matrix<mpfr::mpreal>* const p, boost::numeric::ublas::matrix<mpfr::mpreal>* const q)
{
	// Fincke-Pohst preprocessing - see [CH] Algorithm 2.7.7

	// step 0: initialize
	unsigned int n = g->size1();

	// step 1: Cholesky decomposition of the matrix g
	boost::numeric::ublas::triangular_matrix<mpq_class, boost::numeric::ublas::upper> r(n,n);		// upper triangular matrix
	std::valarray<mpq_class> d(n);																	// diagonal entries
	LALA::choleskyDecomposition(g, &r, &d);

	// temp variables to hold floating point cholesky, ie U^t * U = g, with U = d^(1/2) * r
	boost::numeric::ublas::triangular_matrix<mpfr::mpreal, boost::numeric::ublas::upper> u(n,n), invu(n,n);
	std::valarray<mpfr::mpreal> invd(n);
	std::valarray<mpfr::mpreal> sqrtd(n);
	mpfr::mpreal tempSqrt;

	// copy r into rFloat
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=i; j<n; j++)
			u(i,j) = r(i,j).__get_mp();

	// compute r^-1
	LALA::inverse(&u, &invu);

	// get square roots and inverse of diagonal elements
	for(unsigned int i=0; i<n; i++)
	{
		tempSqrt = d[i].__get_mp();
		mpfr_sqrt(tempSqrt.mpfr_ptr(), tempSqrt.mpfr_srcptr(), MPFR_RNDN);
		sqrtd[i] = tempSqrt;
		invd[i] = 1/sqrtd[i];
	}

	// compute u = d^(1/2) * r
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=i; j<n; j++)
			u(i,j) *= sqrtd[i];

	// compute u'-1 = r^-1 * d^(1/2)^-1
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=i; j<n; j++)
			invu(i,j) *= invd[j];

	// step 2: LLL reduction of the vectors formed by the rows of invu

	// copy transpose into matrix
	boost::numeric::ublas::matrix<mpfr::mpreal> uLLL(n,n);
	uLLL.assign(boost::numeric::ublas::identity_matrix<mpfr::mpreal> (n));
	uLLL = boost::numeric::ublas::trans(invu);

	// get reduced basis
	boost::numeric::ublas::matrix<mpfr::mpreal> hInv(n,n);
	LALA::lllBasis(&uLLL, &hInv);
	hInv = boost::numeric::ublas::trans(hInv);

	// get h --- todo: change algorithm to compute h from the LLL algorithm
	h->resize(n,n,false);
	LALA::inverse(&hInv, h);

	// step 3: re-order the columns of S (inverse of lll reduced invu) --- todo: change algorithm to not use inversion

	// transpose again
	boost::numeric::ublas::matrix<mpfr::mpreal> S(n,n), invS(n,n);
	invS = boost::numeric::ublas::trans(uLLL);

	// invert --- todo: change algorithm to compute h from the LLL algorithm
	S = boost::numeric::ublas::prod(u,*h);

	// get the norms of the rows of S^-1 and put them into a vector
	std::vector<std::pair<mpfr::mpreal, int> > rowNorms;
	for(unsigned int i=0; i<n; i++)
	{
		mpfr::mpreal norm = 0;
		for(unsigned int j=0; j<n; j++)
			norm += invS(i,j)*invS(i,j);
		rowNorms.push_back(std::pair<mpfr::mpreal, int>(norm,i));
	}

	// sort vector --- to do: write compare functions for those pairs
	std::sort(rowNorms.begin(), rowNorms.end());
	std::reverse(rowNorms.begin(), rowNorms.end());

	// re-order the columns of S in the same order
	boost::numeric::ublas::matrix<mpfr::mpreal> copyS(S);
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=0; j<n; j++)
			S(j, i) = copyS(j, rowNorms[i].second);

	// create permutation matrix
	p->resize(n,n,false);
	for(unsigned int i=0; i<n; i++)
		(*p)(rowNorms[i].second,i) = 1;

	// step 4: get new gram matrix and do cholesky again
	copyS = boost::numeric::ublas::trans(S);
	S = boost::numeric::ublas::prod(copyS,S);

	// init cholesky
	q->resize(n,n, false);
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=i; j<n; j++)
			(*q)(i,j) = S(i,j);

	// cholesky i-loop
	for(unsigned int i=0; i<n; i++)
	{
		for(unsigned int j=i+1; j<n; j++)
		{
			(*q)(j,i) = (*q)(i,j);
			(*q)(i,j) = (*q)(i,j) / (*q)(i,i);
		}

		// main cholesky loop
		for(unsigned int k=i+1; k<n; k++)
			for(unsigned int l=k; l<n; l++)
				(*q)(k,l) -= (*q)(k,i)*(*q)(i,l);
	}
}

/*!
 *   @brief Enumerates all lattice vectors of a given lattice with a length smaller or equal to a given constant.
 *
 *   @param[in] q A constant pointer to a constant matrix of mpreals specifying a quadratic form.
 *   @param[in] c A constant multi-precision integer.
 *   @param[out] X A constant pointer to a vector of pairs of a multi-precision integer (length) and a valarray (short lattice vector); holds all enumerated short vectors.
 *
 * Given a quadratic form q and an integer c, this algorithm enumerates all vectors \f$x\f$ with \f$q(x) \leq c\f$.<br>
 * Based on [CH] Algorithm 2.7.5
 *
 *   @return An Expected with an unsigned long specifying the number of enumerated short vectors, or an Expected with an exception if an error occurred.
 *
 *   @exception Throws a std::runtime_error if a negative bound is encountered.
 */
util::Expected<unsigned long> LALA::shortVectorsPrivate(const boost::numeric::ublas::matrix<mpfr::mpreal>* const q, const mpz_class c, std::vector<std::pair<mpz_class, std::valarray<mpz_class> > >* const X)
{
	// given quadratic form q and an integer c, this algorithm outputs all vectors x with q(x) <= c
	// see [CH] Algorithm 2.7.5
	// step 0: initialize
	unsigned int n = q->size1();
	unsigned int i = n-1;							// main loop variable
	std::valarray<mpfr::mpreal> L(n), T(n), U(n);	// bounds
	std::valarray<mpfr::mpreal> x(n);				// solution
	X->clear();										// this vector holds all the vectors with norm <= c
	T[i] = c.__get_mp(); U[i] = 0;

	bool newVector = false, fromElse = false;		// loop conditions
	for(;;)
	{
		if(!newVector && !fromElse)
		{
			// step 1: compute new bounds

			// Z = sqrt(T_i / q_i,i)
			mpfr::mpreal Z = T[i] / (*q)(i,i);
			Z += LALA::SHVEC_ELLIPSOID_EPSILON;

#ifndef NDEBUG
			if(Z < 0)
				return std::runtime_error("Critical error in LALA::shortVectors: Negative bound!");
#endif
			Z = mpfr::sqrt(Z, MPFR_RNDN);

			// L_i = floor(Z-U_i)
			L[i] = mpfr::floor(Z-U[i]);

			// x_i = ceil(-Z-U_i)-1
			x[i] = mpfr::ceil(-Z-U[i])-1;
		} // end if not from new vector and not from else

		// step 2: main loop
		newVector = false; fromElse = false;
		x[i] = x[i] + 1;

		if(x[i] <= L[i])
		{
			// step 3: change bounds
			if(i != 0)
			{

				T[i-1] = T[i] - (*q)(i,i)*(x[i]+U[i])*(x[i]+U[i]);
				i--;
				U[i] = 0;
				for(unsigned int j=i+1; j<n; j++)
					U[i] += (*q)(i,j)*x[j];
			}
			else
			{
				// step 4: if x==0, terminate the algorithm
				bool nul = true;
				for(unsigned int k=0; k<n; k++)
					if(x[k] != 0)
					{
						nul = false;
						break;
					}
				if(nul)
					return X->size();

				// otherwise output x and q(x)
				mpz_class nx;
				mpfr_get_z(nx.get_mpz_t(), (c.__get_mp()-T[0]+(*q)(0,0)*(x[0]+U[0])*(x[0]+U[0])).mpfr_ptr(), MPFR_RNDN);
				std::valarray<mpz_class> lx(n);
				for(unsigned int k=0; k<n; k++)
					mpfr_get_z(lx[k].get_mpz_t(), x[k].mpfr_srcptr(), MPFR_RNDN);
				X->push_back(std::pair<mpz_class, std::valarray<mpz_class>>(nx,lx));
				newVector = true;
			}
		}
		else
		{
			fromElse = true;
			i++;
		}
	}
	return X->size();
}

// public functions
/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Linear Independence //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Tests a vector system for linear independency over \f$\mathbb{R}\f$.
 *
 *   @param vectorSystem A constant pointer to a constant vector of valarrays of multi-precision integers.
 *
 *   @return True if and only if the given vector system is linearly independent over \f$\mathbb{R}\f$.
 *
 *   Uses Gauss-Bareiss for rank determination.
 */
bool LALA::linearlyIndependentR(const std::vector<std::valarray<mpz_class> >* const vectorSystem)
{
	assert(vectorSystem->size() > 0);

	// copy vectors into a matrix
	unsigned int nVectors = vectorSystem->size(), dim = (*vectorSystem)[0].size();

	boost::numeric::ublas::matrix<mpfr::mpreal> vectorSystemMatrix(dim, nVectors);
	for(unsigned int i=0; i<nVectors; i++)
		for(unsigned int j=0; j<dim; j++)
			vectorSystemMatrix(j,i) = (*vectorSystem)[i][j].get_mpz_t();

	// use gauss-bareiss for rank determination
	return LALA::rank(&vectorSystemMatrix) == vectorSystem->size();
}

/*!
 *   @brief Tests a vector system for linear independency over \f$\mathbb{Z}\f$.
 *
 *   @param vectorSystem A constant pointer to a constant vector of valarrays of multi-precision integers.
 *
 *   @return True if and only if the given vector system is linearly independent over \f$\mathbb{Z}\f$.
 */
bool LALA::linearlyIndependentZ(const std::vector<std::valarray<mpz_class> >* const vectorSystem)
{
	// copy vectors into matrix
	boost::numeric::ublas::matrix<mpz_class> m(vectorSystem[0].size(), vectorSystem->size());
	for(unsigned int i=0; i<m.size2(); i++)
		for(unsigned int j=0; j<m.size1(); j++)
			m(j,i) = (*vectorSystem)[i][j];

	// compute rank of the matrix
	unsigned long rk = LALA::rank(&m);

	if(rk < vectorSystem->size())
		return false;
	else
		return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Element Of /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Tests whether a given vector is in the linear span of a vector system or not.
 *
 *   @param b A constant pointer to a constant vector of valarrays of multi-precision integers. This is the vector system.
 *   @param v A constant pointer to a constant valarray of multi-precision integers. This is the vector to test.
 *
 *   @return True if and only if \f$v \in \operatorname{Lin}\{b_1, \dots, b_n\}\f$, where \f$n\f$ denotes the number of vectors in the vector system b.
 */
bool LALA::inLinearSpanZ(const std::vector<std::valarray<mpz_class> >* const b, const std::valarray<mpz_class>* const v)
{
	// create working copy and add v to the vectorSystem
	std::vector<std::valarray<mpz_class> > vs((*b));
	vs.push_back((*v));
	return !LALA::linearlyIndependentZ(&vs);
}

/*!
 *   @brief Tests whether a given vector is in a given lattice or not.
 *
 *   @param[in] m A constant pointer to a constant matrix of multi-precision rationals. This is the inverse of the gram matrix of the lattice under question.
 *   @param[in] v A constant pointer to a constant valarray of multi-precision integers. This is the vector to test.
 *   @param[in] dim A constant unsigned integer specifying the dimension of the lattice.
 *   @param[out] mv A constant pointer to a vector of multi-precision rationals. Holds \f$m \cdot v\f$, if that result is desired; else NULL may be passed.
 *
 *   @return True if and only if \f$v \in L\f$, where \f$L\f$ denotes the lattice with Gram matrix \f$m^{-1}\f$.
 *
 *   See [HB] 3.4.1 for further details.
 */
bool LALA::elementOfLattice(const boost::numeric::ublas::matrix<mpq_class>* const m, const std::valarray<mpz_class>* const v, const unsigned int dim, boost::numeric::ublas::vector<mpq_class>* const mv)
{
	// see [HB] 3.4.1
	unsigned int n = v->size();

	// compute m*v
	boost::numeric::ublas::vector<mpq_class> ublasV(n);
	boost::numeric::ublas::matrix<mpq_class> boostM((*m));

	for(unsigned int i=0; i<n; i++)
		ublasV(i) = (*v)[i];

	ublasV = boost::numeric::ublas::prod(boostM,ublasV);

	if(mv != nullptr)
		(*mv) = ublasV;

	// the first dim coefficients must be integral
	for(unsigned int i=0; i<dim; i++)
		if(ublasV(i).get_den() != 1)
			return false;

	// all the other coefficients must be 0
	for(unsigned int i=dim; i<n; i++)
		if(ublasV(i) != 0)
			return false;

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Public Scalar Products ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the scalar product of two vectors.
 *
 *   @param x Constant pointer to a constant valarray of multi-precision integers.
 *   @param y Constant pointer to a constant valarray of multi-precision integers.
 *   @param innerProduct Constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return The value of the scalar product of \f$x\f$ and \f$y\f$, i.e. either \f$x^t \cdot y\f$ or \f$x^t \cdot A \cdot y\f$, where \f$A\f$ denoted the given inner product.
 *
 *   In the case where no inner product is specified, the loops are unrolled manually for a significant speed burst.
 *
 *   @todo Unroll the loops in the presence of an inner product as well. Register blocking or prefetching will probably lead to a great increase in speed.
 */
mpz_class LALA::scalarProduct(const std::valarray<mpz_class>* const x, const std::valarray<mpz_class>* const y, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	mpz_class result = 0;

	if(innerProduct == nullptr)
	{
		unsigned int n = x->size();

		if(n%4==0)
		{
			mpz_class dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0;
			for(unsigned int i = 0; i < n; i += 4)
			{
				dot1 += (*x)[i] * (*y)[i];
				dot2 += (*x)[i + 1] * (*y)[i + 1];
				dot3 += (*x)[i + 2] * (*y)[i + 2];
				dot4 += (*x)[i + 3] * (*y)[i + 3];
			}
			result = dot1 + dot2 + dot3 + dot4;
		}
		else if(n%2==0)
		{
			mpz_class dot1 = 0, dot2 = 0;
			for(unsigned int i = 0; i < n; i += 2)
			{
				dot1 += (*x)[i] * (*y)[i];
				dot2 += (*x)[i + 1] * (*y)[i + 1];
			}
			result = dot1 + dot2;
		}
		else
		{
			mpz_class dot1 = 0, dot2 = 0;
			unsigned int nn = n/2;
			for(unsigned int i = 0; i < n && nn>0; i += 2, nn--)
			{
				dot1 += (*x)[i] * (*y)[i];
				dot2 += (*x)[i + 1] * (*y)[i + 1];
			}
			result = dot1 + dot2 + (*x)[n-1] * (*y)[n-1];
		}
	}
	else
	{
		assert(x->size() == y->size() && innerProduct->size2() == y->size() && x->size() == innerProduct->size1());
		unsigned int n = x->size();

		// copy valarrays into boost vectors
		boost::numeric::ublas::vector<mpz_class> copyX(x->size()), copyY(y->size());
		for(unsigned int i=0; i<n; i++)
		{
			copyX(i) = (*x)[i];
			copyY(i) = (*y)[i];
		}

		copyY = boost::numeric::ublas::prod((*innerProduct), copyY);
		result = boost::numeric::ublas::inner_prod(copyX, copyY);
	}

	return result;
}

/*!
 *   @brief Computes the scalar product of two vectors.
 *
 *   @param x Constant pointer to a constant valarray of multi-precision rationals.
 *   @param y Constant pointer to a constant valarray of multi-precision rationals.
 *   @param innerProduct Constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return The value of the scalar product of \f$x\f$ and \f$y\f$, i.e. either \f$x^t \cdot y\f$ or \f$x^t \cdot A \cdot y\f$, where \f$A\f$ denoted the given inner product.
 *
 *   In the case where no inner product is specified, the loops are unrolled manually for a significant speed burst.
 *
 *   @todo Unroll the loops in the presence of an inner product as well. Register blocking or prefetching will probably lead to a great increase in speed.
 */
mpq_class LALA::scalarProduct(const std::valarray<mpq_class>* const x, const std::valarray<mpq_class>* const y, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	mpq_class result = 0;

	if(innerProduct == nullptr)
	{
		unsigned int n = x->size();

		if(n%4==0)
		{
			mpq_class dot1 = 0, dot2 = 0, dot3 = 0, dot4 = 0;
			for(unsigned int i = 0; i < n; i += 4)
			{
				dot1 += (*x)[i] * (*y)[i];
				dot2 += (*x)[i + 1] * (*y)[i + 1];
				dot3 += (*x)[i + 2] * (*y)[i + 2];
				dot4 += (*x)[i + 3] * (*y)[i + 3];
			}
			result = dot1 + dot2 + dot3 + dot4;
		}
		else if(n%2==0)
		{
			mpq_class dot1 = 0, dot2 = 0;
			for(unsigned int i = 0; i < n; i += 2)
			{
				dot1 += (*x)[i] * (*y)[i];
				dot2 += (*x)[i + 1] * (*y)[i + 1];
			}
			result = dot1 + dot2;
		}
		else
		{
			mpq_class dot1 = 0, dot2 = 0;
			unsigned int nn = n/2;
			for(unsigned int i = 0; i < n && nn>0; i += 2, nn--)
			{
				dot1 += (*x)[i] * (*y)[i];
				dot2 += (*x)[i + 1] * (*y)[i + 1];
			}
			result = dot1 + dot2 + (*x)[n-1] * (*y)[n-1];
		}
	}
	else
	{
		assert(x->size() == y->size() && innerProduct->size2() == y->size() && x->size() == innerProduct->size1());
		unsigned int n = x->size();

		// copy valarrays into boost vectors
		boost::numeric::ublas::vector<mpq_class> copyX(x->size()), copyY(y->size());
		for(unsigned int i=0; i<n; i++)
		{
			copyX(i) = (*x)[i];
			copyY(i) = (*y)[i];
		}

		copyY = boost::numeric::ublas::prod((*innerProduct), copyY);
		result = boost::numeric::ublas::inner_prod(copyX, copyY);
	}

	return result;
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Gram-Schmidt ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Performs Gram-Schmidt-Orthogonalization.
 *
 *   @param basis Constant pointer to vector of valarrays: the given basis to orthogonalize.
 *   @param g Constant pointer to a constant symmetric matrix specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return void
 *
 *   Computes an orthogonal basis w.r.t. standard euclidean basis or g, starting with the given basis.
 */
void LALA::gramSchmidt(std::vector<std::valarray<mpq_class> >* const basis, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g)
{
	unsigned int dim = basis->size();

	std::vector<std::valarray<mpq_class> > vHut((*basis));

	boost::numeric::ublas::matrix<mpq_class> r(dim,dim);
	std::vector<std::valarray<mpq_class> > q(dim);
	for(unsigned int i=0; i<dim; i++)
		q[i].resize(dim);

	for(unsigned int i=0; i<dim; i++)
	{
		if(g != nullptr)
			r(i,i) = LALA::scalarProduct(&vHut[i],&vHut[i], g);
		else
			r(i,i) = LALA::scalarProduct(&vHut[i],&vHut[i]);

		for(unsigned int j=0; j<dim; j++)
			q[i][j] = vHut[i][j] / r(i,i);

		for(unsigned int j=i+1; j<dim; j++)
		{
			if(g != nullptr)
				r(i,j) = LALA::scalarProduct(&q[i], &vHut[j], g);
			else
				r(i,j) = LALA::scalarProduct(&q[i], &vHut[j]);

			for(unsigned int k=0; k<dim; k++)
				vHut[j][k] -= r(i,j)*q[i][k]*r(i,i);
		}
	}

	// copy the basis back
	(*basis) = vHut;
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// FP-LLL /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Floating point LLL.
 *
 *   @param[in, out] b Constant pointer to a matrix of mpfreals. The input is the basis to reduce and the output is the LLL-reduced basis.
 *   @param[out] h Constant pointer to a matrix of mpreals. If not a null pointer, this will be the basis change matrix from the original basis to the LLL-reduced basis.
 *
 *   @return A void Expected, or a void Expected with a std::runtime_error exception if the given vector system was not a basis.
 *
 *   After the function returns, b will be lll reduced and h holds the coordinates of the new LLL-reduced basis in terms of the old one.
 *   See [CH] Definition 2.6.1 for further details.
 */
util::Expected<void> LALA::lllBasis(boost::numeric::ublas::matrix<mpfr::mpreal>* const b, boost::numeric::ublas::matrix<mpfr::mpreal>* const h)
{
	// the columns of b encode the base vectors b_1, ..., b_n
	// b will be lll reduced - see [CH] Definition 2.6.1
	// h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one

	// step 1: initialize the algorithm

	// "global" variables
	int dim = b->size1();
	boost::numeric::ublas::matrix<mpfr::mpreal> bStar(dim,dim);	// holds the Gram-Schmidt vectors
	std::valarray<mpfr::mpreal> B(dim);		// b_i = det(<b_r, b_s>)_1<=r,s<=i)
	std::valarray<mpfr::mpreal> mu(dim*dim);// mu_i,j = d_j * u_i,j , where u_i,j is the coefficient from the Gram-Schmidt orthogonalisation

	// set initial values
	int k=1, kmax=0;

	for(int i=0; i<dim; i++)
		bStar(i,0) = (*b)(i,0);

	B[0] = LALA::lllBasisDot(b,0,0);

	if(h != nullptr)
		h->assign(boost::numeric::ublas::identity_matrix<mpfr::mpreal> (b->size1()));

	// step 2: Gram-Schmidt
	while(k<dim)
	{
		if(k > kmax)	// if we haven't worked with this column yet, do so now
		{
			kmax = k;

			// set b*_k = b_k
			for(int i=0; i<dim; i++)
				bStar(i,k) = (*b)(i,k);

			for(int j=0; j<k; j++)
			{
				// mu_{k,j} = <b_k, b*_j> / B_j
				mu[k*dim+j] = LALA::lllBasisDot(b, &bStar, k, j) / B[j];

				// set b*_k = b*_k - mu_{k,j}*b*_j
				for(int i=0; i<dim; i++)
					bStar(i,k) = bStar(i,k) - mu[k*dim+j]*bStar(i,j);
			}

			// set B_k = b*_k * b*_k
			B[k] = LALA::lllBasisDot(&bStar, k, k);

			if(B[k] == 0)
			{
				// linear dependency
				return std::runtime_error("Critical error in LLL: given vector system is linearly dependent!");
			}
		}

		// step 3: test LLL condition
		LALA::lllBasisTestCondition(b, &bStar, h,&B,&mu, &k, kmax, dim);

		// reduce columns
		for(int i=k-2; i>=0; i--)
			LALA::lllBasisRed(b,h,&mu,k,i, dim);
		k++;
	}

	// return success
	return { };
}

/*!
 *   @brief Integral modified LLL.
 *
 *   @param[in, out] b Constant pointer to a matrix of multi-precision integers, whose columns are, as input, the basis to reduce and, as output, the LLL-reduced basis.
 *   @param[out] h Constant pointer to a matrix of multi-precision integers. If not a null pointer, this will be the basis change matrix from the original basis to the LLL-reduced basis.
 *   @param[in] innerProduct A constant pointer to a constant symmetric matrix of multi-precision integers specifying a scalar product. This is NULL by default (standard scalar product).
 *
 *   @return An unsigned int: the rank of the given vector system.
 *
 *	 The columns of b encode a generating system. b will be a LLL-reduced basis. See [CH] Algorithm 2.6.8 for further details.
 *	 If so desired, h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one.
 */
unsigned int LALA::mlllBasis(boost::numeric::ublas::matrix<mpz_class>* const b, boost::numeric::ublas::matrix<mpz_class>* const h, boost::numeric::ublas::symmetric_matrix<mpz_class>* const innerProduct)
{
	// the columns of b encode a generating system
	// b will be a lll reduced basis - see [CH] Algorithm 2.6.8
	// h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one

	// step 1: initialize the algorithm

	// "global" variables
	int dim = b->size1();
	int nVectors = b->size2();
	boost::numeric::ublas::matrix<mpfr::mpreal> bStar(dim,nVectors);	// holds the Gram-Schmidt vectors
	std::valarray<mpfr::mpreal> B(nVectors);							// b_i = det(<b_r, b_s>)_1<=r,s<=i)
	std::valarray<mpfr::mpreal> mu(nVectors*nVectors);					// mu_i,j = d_j * u_i,j , where u_i,j is the coefficient from the Gram-Schmidt orthogonalisation

	// set initial values
	int k=1, kmax=0;

	for(int i=0; i<dim; i++)
		bStar(i,0) = (*b)(i,0).get_mpz_t();

	if(innerProduct == nullptr)
	{
		mpfr::mpreal dot;
		mpfr_set_z(dot.mpfr_ptr(), LALA::lllBasisDot(b,0,0).get_mpz_t() , MPFR_RNDN);
		B[0] = dot;
	}
	else
	{
		mpfr::mpreal dot;
		mpfr_set_z(dot.mpfr_ptr(), LALA::lllBasisDot(b,0,0, innerProduct).get_mpz_t() , MPFR_RNDN);
		B[0] = dot;
	}

	if(h != nullptr)
		h->assign(boost::numeric::ublas::identity_matrix<mpz_class> (b->size1()));

	// step 2: Gram-Schmidt
	while(k<nVectors)
	{
		if(k > kmax)	// if we haven't worked with this column yet, do so now
		{
			kmax = k;

			// set b*_k = b_k
			for(int i=0; i<dim; i++)
				bStar(i,k) = (*b)(i,k).get_mpz_t();

			for(int j=0; j<k; j++)
			{
				// mu_{k,j} = <b_k, b*_j> / B_j
				if(B[j] != 0)
				{
					if(innerProduct == nullptr)
						mu[k*nVectors+j] = LALA::lllBasisDot(b, &bStar, k, j) / B[j];
					else
						mu[k*nVectors+j] = LALA::lllBasisDot(b, &bStar, k, j, innerProduct) / B[j];
				}
				else
					mu[k*nVectors+j] = 0;

				// set b*_k = b*_k - mu_{k,j}*b*_j
				for(int i=0; i<dim; i++)
					bStar(i,k) = bStar(i,k) - mu[k*nVectors+j]*bStar(i,j);
			}

			// set B_k = b*_k * b*_k
			if(innerProduct == nullptr)
				B[k] = LALA::lllBasisDot(&bStar, k, k);
			else
				B[k] = LALA::lllBasisDot(&bStar, k, k, innerProduct);
		}

		// step 3: test LLL condition
		LALA::mlllBasisTestCondition(b, &bStar, h,&B,&mu, &k, kmax, dim, nVectors);

		// reduce columns
		for(int i=k-2; i>=0; i--)
			LALA::mlllBasisRed(b,h,&mu,k,i, nVectors);
		k++;
	}

	// get rank of the vector system
	int rk = 0;
	for(unsigned int i=0; i<b->size2(); i++)
	{
		bool zeroCol = true;
		for(unsigned int j=0; j<b->size1(); j++)
			if((*b)(j,i) != 0)
				zeroCol = false;
		if(!zeroCol)
			rk++;
	}

	// delete zero columns
	MatrixHelper<mpz_class>::baseDelZeroes(b, rk);

	return rk;
}

/*!
 *   @brief Floating point LLL.
 *
 *   @param[in, out] g Constant pointer to a symmetric matrix of mpfreals. The input is a Gram matrix and the output is the LLL-reduced Gram matrix.
 *   @param[out] h Constant pointer to a matrix of mpreals. If not a null pointer, this will be the basis change matrix from the original basis to the LLL-reduced basis.
 *
 *   @return A void Expected, or a void Expected with a std::runtime_error exception if the given vector system was not a basis.
 *
 *  g is the gramian of a basis \f$v_1, \dots, v_n\f$. g will be LLL-reduced. See [CH] Definition 2.6.1 for further details.<br>
 *  If so desired, h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one.
 */
util::Expected<void> LALA::lllGram(boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const g, boost::numeric::ublas::matrix<mpfr::mpreal>* const h)
{
	// g is the gramian of a basis v_1, ..., v_n
	// g will be lll reduced - see [CH] Definition 2.6.1
	// h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one

	// step 1: initialize the algorithm

	// "global" variables
	int dim = g->size1();
	std::valarray<mpfr::mpreal> a(dim*dim);	// auxiliary array to minimize the number of gram-schmidt operations
	std::valarray<mpfr::mpreal> b(dim);		// b_i = det(<b_r, b_s>)_1<=r,s<=i)
	std::valarray<mpfr::mpreal> mu(dim*dim);	// mu_i,j = d_j * u_i,j , where u_i,j is the coefficient from the gram-schmidt orthogonalisation
	mpfr::mpreal u = 0, atemp, btemp;			// u, atemp and btemp are used in the gram-schmidt process

	// set initial values
	b[0]=(*g)(0,0);
	int k=1, kmax=0;

	if(h != nullptr)
		h->assign(boost::numeric::ublas::identity_matrix<mpfr::mpreal> (g->size1()));

	// step 2: Gram-Schmidt
	while(k<dim)
	{
		if(k > kmax)	// if we haven't worked with this column yet, do so now
		{
			kmax = k;
			for(int j=0; j<k; j++)
			{
				atemp = (*g)(k,j);
				for(int i=0; i<j; i++)
					atemp -= mu[j*dim+i]*a[k*dim+i];
				a[k*dim+j] = atemp;
				mu[k*dim+j] = atemp / b[j];
			}
			btemp = (*g)(k,k);
			for(int i=0; i<k; i++)
				btemp -= mu[k*dim+i]*a[k*dim+i];
			if(btemp != 0)
				b[k] = btemp;
			else
			{
				// linear dependency
				return std::runtime_error("Critical error in LLL: given vector system is linearly dependent!");
			}
		}

		// step 3: test LLL condition
		LALA::lllGramTestCondition(g,h,&b,&mu, &k, kmax, dim);

		// reduce columns
		for(int i=k-2; i>=0; i--)
			LALA::lllGramRed(g,h,&mu,k,i, dim);
		k++;
	}

	// return success
	return { };
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Integral LLL ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Integral LLL.
 *
 *   @param[in, out] g Constant pointer to a symmetric matrix of multi-precision integers. The input is a Gram matrix and the output is the LLL-reduced Gram matrix.
 *   @param[out] h Constant pointer to a matrix of multi-precision integers. If not a null pointer, this will be the basis change matrix from the original basis to the LLL-reduced basis.
 *
 *   @return A void Expected, or a void Expected with a std::runtime_error exception if the given vector system was not a basis.
 *
 *	 g is the gramian of a basis \f$v_1, \dots, v_n\f$. g will be LLL-reduced - see [CH] Definition 2.6.1 for further information.<br>
 *	 If so desired, h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one.
 */
util::Expected<void> LALA::lllGram(boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, boost::numeric::ublas::matrix<mpz_class>* const h)
{
	// g is the gramian of a basis v_1, ..., v_n
	// g will be lll reduced - see [CH] Definition 2.6.1
	// h will be the transformation matrix, i.e. h holds the coordinates of the new basis in terms of the old one

	// step 1: initialize the algorithm

	// "global" variables
	int dim = g->size1();
	std::valarray<mpz_class> d(dim+1);				// d_i = det(<b_r, b_s>)_1<=r,s<=i)
	std::valarray<mpz_class> lambda(dim*dim);		// lambda_i,j = d_j * u_i,j , where u_i,j is the coefficient from the gram-schmidt orthogonalisation
	mpz_class u = 0, num, den;						// u is used to compute determinants and gram-schmidt coefficients ; num and den are used in euclidean division

	// set initial values
	d[0]=1; d[1]=(*g)(0,0);
	int k=1, kmax=0;

	if(h != nullptr)
		h->assign(boost::numeric::ublas::identity_matrix<mpz_class> (g->size1()));

	// step 2: Gram-Schmidt
	while(k<dim)
	{
		if(k > kmax)	// if we haven't worked with this column yet, do so now
		{
			kmax = k;
			for(int j=0; j<=k; j++)
			{
				u = (*g)(k,j);

				for(int i=0; i<j; i++)
				{
					num = (d[i+1]*u - lambda[k*dim+i] * lambda[j*dim+i]);
					den = d[i];
					mpz_tdiv_q(u.get_mpz_t(), num.get_mpz_t(), den.get_mpz_t());
				}

				if(j<k)
					lambda[k*dim+j]=u;			// gram-schmidt coefficient
				else
				{
					if(u==0)
					{
						// sub-determinant zero => linear dependency
						return std::runtime_error("Critical error in LLL: given vector system is linearly dependent!");
					}
					d[k+1]=u;				// sub-determinant
				}
			}
		}

		// step 3: test LLL condition
		LALA::lllGramTestCondition(g,h,&d,&lambda, &k, kmax, dim);

		// reduce columns
		for(int i=k-2; i>=0; i--)
			LALA::lllGramRed(g,h,&d,&lambda,k,i, dim);
		k++;
	}

	// return success
	return { };
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Short Vectors //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Enumerates all vector in a given lattice with a shorter length than a given constant.
 *
 *   @param[in] g Constant pointer to a constant symmetric matrix of multi-precision integers specifying a positive definite quadratic form.
 *   @param[in] c Constant positive multi-precision integer.
 *   @param[out] X Constant pointer to a vector of pairs of multi-precision integers and valarrays of the same type, holding all enumerated vectors. The first entry is the length of the vector, the second the vector itself.
 *
 *   @return An Expected with an unsigned long, the number of short vectors enumerated by the algorithm, a std::runtime_error exception an error was encountered during the computation or a std::invalid_argument exception if c is zero or negative.
 *
 *	 This algorithm outputs all non-zero vectors \f$x \in \mathbb{Z}\f$ such that \f$q(x) < c\f$ and the corresponding values of \f$q(x)\f$.
 *
 *   @todo Check for positive definiteness.
 */
util::Expected<unsigned long> LALA::shortVectors(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g, const mpz_class c, std::vector<std::pair<mpz_class, std::valarray<mpz_class>> >* const X)
{
	unsigned int n = g->size1();

	if(c <= 0)
		return std::invalid_argument("Invalid argument in LALA::shortVectors: c <= 0");

	// Fincke-Pohst preprocessing - see [CH] Algorithm 2.7.7
	boost::numeric::ublas::matrix<mpfr::mpreal> q, h, p;
	LALA::shortVectorsFinckePohst(g, &h, &p, &q);

	// short vector algorithm - see [CH] Algorithm 2.7.5
	unsigned long nx = 0;
	try{ nx = LALA::shortVectorsPrivate(&q, c, X).get(); }
	catch(std::runtime_error& e)
	{
		return e;
	}

	// sort the list, norm the list, i.e. first nonzero entry should be positive and re-transform the vectors with h and p
	std::sort(X->begin(), X->end(), LALA::comparePairMPZValarray);
	for(unsigned int i=0; i<X->size(); i++)
	{
		// copy solution vector into a boost vector
		boost::numeric::ublas::vector<mpfr::mpreal> boostX(n);
		for(unsigned int j=0; j<n; j++)
			boostX(j) = (*X)[i].second[j].__get_mp();

		// compute P*x
		boostX = boost::numeric::ublas::prod(p, boostX);

		// compute h*P*x
		boostX = boost::numeric::ublas::prod(h,boostX);

		// copy the result back into X - normalize the vectors, i.e. first nonzero entry should be positive
		bool firstNonZeroCoordinateNegative = false, firstNonZeroCoordinate = true;
		for(unsigned int j=0; j<n; j++)
		{
			mpz_class t;
			mpfr_get_z(t.get_mpz_t(), boostX(j).mpfr_srcptr(), MPFR_RNDN);
			if(t!=0 && firstNonZeroCoordinate)
			{
				firstNonZeroCoordinate = false;
				if(t<0)
					firstNonZeroCoordinateNegative = true;
			}
				if(firstNonZeroCoordinateNegative)
					(*X)[i].second[j] = -t;
				else
					(*X)[i].second[j] = t;
		}
	}

	return nx;
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Basis from Generating Set //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Constructs a lattice basis from a generating set.
 *
 *   @param[in] gen Constant pointer to a constant vector of valarrays with multi-precision integers specifying a generating system for a lattice.
 *   @param[out] basis Constant pointer to a vector of valarrays of multi-precision integers. This will out the basis constructed from the generating set.
 *
 *   @return void
 *
 *	This algorithm computes a lattice basis from a given generating system. <br>
 *	This is a serial algorithm with pruning. See [HB] chapter 3.5 for further information.
 */
void LALA::basisFromGeneratingSystem(const std::vector<std::valarray<mpz_class> >* const gen, std::vector<std::valarray<mpz_class> >* const basis)
{
	// algorithm to compute a lattice basis from a given generating system
	// serial algorithm with pruning
	// see [HB] chapter 3.5 for further information

	// create working copy
	std::vector<std::valarray<mpz_class> > copyGen = (*gen);

	// step1: initialize
	unsigned long dim = copyGen[0].size();

	// test the first dim vectors of the generating set

	// partial basis matrix - with basis extension
	boost::numeric::ublas::matrix<mpq_class> invBasisMatrix(dim,dim);

	for(unsigned int i=0; i<dim; i++)
		for(unsigned int j=0; j<dim; j++)
			invBasisMatrix(j,i) = copyGen[i][j];

	if(LALA::rank(&invBasisMatrix) == dim)
	{
		// found a basis
		for(unsigned int i=0; i<dim; i++)
			basis->push_back(copyGen[i]);

		return;
	}

	// chose v in gen and then delete v from gen
	std::valarray<mpz_class> v = copyGen[copyGen.size()-1];
	copyGen.pop_back();

	// set first basis vector to v
	basis->clear();
	basis->push_back(v);

	invBasisMatrix.clear();
	// extend to vector space basis
	for(unsigned int i=0; i<dim; i++)
		invBasisMatrix(i,0) = v[i];

	LALA::basisExtensionZ(&invBasisMatrix);

	// invert matrix
	LALA::inverse(&invBasisMatrix, &invBasisMatrix);

	// step2: successively add vectors from the generating set to the basis
	while(!copyGen.empty())
	{
		v = (*gen)[gen->size()-1];
		copyGen.pop_back();

		// if v not in <b_1, ..., b_k>_Z
		if(!LALA::elementOfLattice(&invBasisMatrix, &v, basis->size()))
		{
			// then basis = mLLL(b_1, ..., b_k, v)

			// copy basis vectors into matrix
			boost::numeric::ublas::matrix<mpz_class> m((*basis)[0].size(), basis->size()+1);
			for(unsigned int i=0; i<m.size2()-1; i++)
				for(unsigned int j=0; j<m.size1(); j++)
					m(j,i) = (*basis)[i][j];

			// add v
			for(unsigned int i=0; i<m.size1(); i++)
				m(i,basis->size()) = v[i];

			// construct basis using the mLLL algorithm
			unsigned long rk = LALA::mlllBasis(&m);

			// copy basis back into the vector system
			basis->clear();
			for(unsigned long i = 0; i<rk; i++)
			{
				std::valarray<mpz_class> w(dim);
				for(unsigned int j=0; j<dim; j++)
					w[j] = m(j,i);
				basis->push_back(w);
			}

			if(basis->size() == dim)
				return;

			// extend to vector space basis
			for(unsigned int i=0; i<dim; i++)
				for(unsigned int j=0; j<dim; j++)
					if(i < basis->size())
						invBasisMatrix(j,i) = (*basis)[i][j];
					else
						invBasisMatrix(j,i) = 0;

			LALA::basisExtensionZ(&invBasisMatrix);

			// invert matrix
			LALA::inverse(&invBasisMatrix, &invBasisMatrix);
		}
	}
}

// matrix algorithms

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Gauss /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the echelon form of a matrix.
 *
 *   @param[in, out] m A constant pointer to a matrix with multi-precision rationals whose columns encode a basis for the underlying module. This will be changed to column-echelon form.
 *   @param[out] basisExtension Constant pointer to a vector unsigned ints, specifying standard basis vectors useable for a basis extensions. Default: NULL, no basis extension is computed.
 *
 *   @return void
 */
void LALA::echelonForm(boost::numeric::ublas::matrix<mpq_class>* const m, std::vector<unsigned int>* const basisExtension)
{
	unsigned int nCols = m->size2(), nRows = m->size1();
	unsigned int lead = 0;

	std::vector<unsigned int> swapalot;
	if(basisExtension != nullptr)
		for(unsigned int i=0; i<nRows; i++)
			swapalot.push_back(i);

	bool zeroRow = false;
	unsigned int i=0;
	for(unsigned int r=0; r<nRows; r++)
	{
		zeroRow = false;

		if(nCols <= lead)
			break;

		i = r;

		// find non-zero entry
		while((*m)(i,lead) == 0)
		{
			i++;
			if(nRows == i)
			{
				i = r;
				lead++;
				if(nCols == lead)
				{
					zeroRow = true;
					break;
				}
			}
		}

		if(!zeroRow)
		{
			if(i != r)
			{
				// swap rows i and r
				for(unsigned int j=0; j<nCols; j++)
					boost::swap((*m)(i,j), (*m)(r,j));
				if(basisExtension != nullptr)
					boost::swap(swapalot[i], swapalot[r]);
			}

			// divide
			if((*m)(r,lead) != 0)
			{
				mpq_class q = (*m)(r,lead);
				for(unsigned int j=r; j<nCols; j++)
					(*m)(r,j) /= q;
			}

			// do gauss
			for(i=0; i<nRows; i++)
				if(i != r)
					MatrixHelper<mpq_class>::rowAddMultiple(m, i, r, -(*m)(i,lead));
			lead++;
		}
	}

	if(basisExtension != nullptr)
	{
		// find suitable vectors for basis extension

		// search for zero rows
		for(unsigned int i=0; i<nRows; i++)
		{
			zeroRow = true;
			for(unsigned int j=i; j<nCols; j++)
				if((*m)(i,j) != 0)
					zeroRow = false;

			if(zeroRow)
				basisExtension->push_back(i);
		}

		for(unsigned int i=0; i<swapalot.size(); i++)
		{
			if(swapalot[i] != i)
			{
				std::vector<unsigned int>::iterator it = std::find(basisExtension->begin(), basisExtension->end(), i);
				if(it != basisExtension->end())
					(*it) = swapalot[i];
			}
		}
	}
}

/*!
 *   @brief Computes the echelon form of a matrix.
 *
 *   @param[in] m A constant pointer to a constant matrix with multi-precision rationals whose columns encode a basis for the underlying module. This will be changed to column-echelon form.
 *   @param[out] basisExtension Constant pointer to a vector unsigned ints, specifying standard basis vectors useable for a basis extensions. Default: NULL, no basis extension is computed.
 *
 *   @return void
 */
void LALA::echelonFormConst(const boost::numeric::ublas::matrix<mpq_class>* const m, std::vector<unsigned int>* const basisExtension)
{
	// create working copy
	boost::numeric::ublas::matrix<mpq_class> copyM((*m));
	LALA::echelonForm(&copyM, basisExtension);
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Rank ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the rank of a matrix.
 *
 *   @param m A constant pointer to a constant matrix with mpreals.
 *
 *   @return An unsigned long specifying the rank of the given matrix.
 */
unsigned long LALA::rank(const boost::numeric::ublas::matrix<mpfr::mpreal>* const m)
{
	unsigned int nCols = m->size2(), nRows = m->size1();
	unsigned long rk, dim;
	nCols > nRows ? rk = nRows : rk = nCols;
	dim = rk;

	// create working copy
	boost::numeric::ublas::matrix<mpfr::mpreal> copyM((*m));

	// step 1: initialize
	mpfr::mpreal pivot = 1;
	bool found = false;
	for(unsigned int k=0; k<dim; k++)
	{
		found = false;

		// step 2: increase k
		pivot = copyM(k,k);

		// step 3: is p = 0?
		if(pivot==0)
		{
			int nonZeroEntryPos = k;

			// get k-th column
			boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<mpfr::mpreal> > col(copyM,k);

			// find non-zero coordinate
			for(typename boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<mpfr::mpreal> >::iterator it = col.begin()+k; it != col.end(); it++, nonZeroEntryPos++)
			{
				if(*it != 0)
				{
					found = true;
					break;
				}
			}
			if(!found)
				rk--;

			if(found)
			{
				// now for j = k, ...,n exchange the nonZeroPosEntry,j-th entry with the k,j-entry
				for(unsigned int j=k; j<nCols; j++)
					boost::swap(copyM(nonZeroEntryPos,j), copyM(k,j));

				// now set s=-s and p=k,k-th entry
				pivot = copyM(k,k);
			}

		}

		// step 4: main step (change rows) - p is now non-zero
		if(pivot != 0)
			for(unsigned int i=k+1; i<nRows; i++)
			{
				boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<mpfr::mpreal> > rowI(copyM,i);
				mpfr::mpreal t = rowI(k) / pivot;
				MatrixHelper<mpfr::mpreal>::rowAddMultiple(&copyM, i, k, -t);
			}
		}
	return rk;
}

/*!
 *   @brief Computes the rank of a matrix.
 *
 *   @param m A constant pointer to a constant matrix with multi-precision rational entries.
 *
 *   @return An unsigned long specifying the rank of the given matrix.
 */
unsigned long LALA::rank(const boost::numeric::ublas::matrix<mpq_class>* const m)
{
	unsigned int nCols = m->size2(), nRows = m->size1();
	unsigned long rk, dim;
	nCols > nRows ? rk = nRows : rk = nCols;
	dim = rk;

	// create working copy
	boost::numeric::ublas::matrix<mpq_class> copyM((*m));

	// step 1: initialize
	mpq_class pivot = 1;
	bool found = false;
	for(unsigned int k=0; k<dim; k++)
	{
		found = false;

		// step 2: get pivot
		pivot = copyM(k,k);

		// step 3: is p = 0?
		if(pivot==0)
		{
			int nonZeroEntryPos = k;

			// get k-th column
			boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<mpq_class> > col(copyM,k);

			// find non-zero coordinate
			for(typename boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<mpq_class> >::iterator it = col.begin()+k; it != col.end(); it++, nonZeroEntryPos++)
			{
				if(*it != 0)
				{
					found = true;
					break;
				}
			}
			if(!found)
				rk--;

			if(found)
			{
				// now for j = k, ...,n exchange the nonZeroPosEntry,j-th entry with the k,j-entry
				for(unsigned int j=k; j<nCols; j++)
					boost::swap(copyM(nonZeroEntryPos,j), copyM(k,j));

				// now set s=-s and p=k,k-th entry
				pivot = copyM(k,k);
			}
		}

		// step 4: main step (change rows) - p is now non-zero
		if(pivot != 0)
		for(unsigned int i=k+1; i<nRows; i++)
		{
			boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<mpq_class> > rowI(copyM,i);
			mpq_class t = rowI(k) / pivot;
			MatrixHelper<mpq_class>::rowAddMultiple(&copyM, i, k, -t);
		}
	}

	return rk;
}

/*!
 *   @brief Computes the rank of a matrix.
 *
 *   @param m A constant pointer to a constant matrix of multi-precision integers.
 *
 *   @return An unsigned long specifying the rank of the given matrix.
 */
unsigned long LALA::rank(const boost::numeric::ublas::matrix<mpz_class>* const m)
{
	// copy matrix into mpq matrix
	boost::numeric::ublas::matrix<mpq_class> copyM((*m));
	return LALA::rank(&copyM);
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Basis Extensions //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes a basis extensions for a given basis.
 *
 *   @param[in,out] m A constant pointer to a matrix with multi-precision rational entries. A possible basis extension is appended to its columns.
 *
 *   @return void
 *
 *   The matrix m should have the correct dimensions (filled with 0 cols).
 */
void LALA::basisExtensionZ(boost::numeric::ublas::matrix<mpq_class>* const m)
{
	// matrix m should have the correct dimensions (filled with 0 cols)
	assert(m->size1() == m->size2());
	unsigned int dim = m->size1();

	// get non pivot columns
	std::vector<unsigned int> basisExtension;

	LALA::echelonFormConst(m, &basisExtension);

	// extend basis
	std::vector<int> e(dim);
	e.clear();
	for(unsigned int i=0; i<basisExtension.size(); i++)
	{
		e[basisExtension[i]] = 1;
		for(unsigned int j=0; j<dim; j++)
			(*m)(j, dim-1-i) = e[j];
		e[basisExtension[i]] = 0;
	}
}

/*!
 *   @brief Computes the orthogonal complement of a submodule.
 *
 *   @param[in] latticeBasis A constant pointer to a constant vector of valarrays with multi-precision rational entries specifying a lattice basis.
 *   @param[out] orthCompl A constant pointer to a vector of valarrays with multi-precision rational entries holding the orthogonal complement of the given lattice basis.
 *   @param[in] g The scaler product to use. Default: NULL (standard scalar product).
 *
 *   @return void
 */
void LALA::orthogonalComplement(const std::vector<std::valarray<mpq_class> >* const latticeBasis, std::vector<std::valarray<mpq_class> >* const orthCompl, const boost::numeric::ublas::symmetric_matrix<mpz_class>* const g)
{
	unsigned int rk = latticeBasis->size();
	unsigned int dim = latticeBasis->at(0).size();

	std::vector<std::valarray<mpq_class> > complement;

	// extend to a basis of the vector space
	boost::numeric::ublas::matrix<mpq_class> m(dim,dim);
	for(unsigned int i=0; i<latticeBasis->size(); i++)
		for(unsigned int j=0; j<dim; j++)
			m(j,i) = latticeBasis->at(i)[j];

	LALA::basisExtensionZ(&m);

	// copy basis into vector system
	complement.clear();
	for(unsigned int i=0; i<dim; i++)
	{
		std::valarray<mpq_class> v(dim);
		for(unsigned int j=0; j<dim; j++)
			v[j] = m(j,i);
		complement.push_back(v);
	}

	// compute orthogonal complement via gram schmidt
	if(g != nullptr)
		LALA::gramSchmidt(&complement,g);
	else
		LALA::gramSchmidt(&complement);

	// get complement
	(*orthCompl) = std::vector<std::valarray<mpq_class> >(complement.begin()+rk, complement.end());
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Matrix Inversion ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the inverse of a given upper triangular matrix.
 *
 *   @param[in] a A constant pointer to a constant upper triangular matrix of mpreals specifying the matrix to be inverted.
 *   @param[out] b A constant pointer to an upper triangular matrix of mpreals. This will be the inverted matrix when the function returns.
 *
 *   @return void
 *
 *   @section Note
 *
 * For an upper triangular matrix \f$A\f$ with \f$A_{ii} = 1\f$ for all \f$i\f$, we have:
 * \f$A = I + N\f$, with \f$N\f$ nilpotent, thus \f$A^{-1} = (I+N)^{-1} = I + \sum\limits_{i=1}^n N^i\f$, where \f$n\f$ is the nilpotent index of \f$A\f$.
 */
void LALA::inverse(const boost::numeric::ublas::triangular_matrix<mpfr::mpreal, boost::numeric::ublas::upper>* const a, boost::numeric::ublas::triangular_matrix<mpfr::mpreal, boost::numeric::ublas::upper>* const b)
{
	// for an upper triangular matrix A with A_ii = 1 for all i, we have: A = I + N, with N nilpotent, thus
	// A^-1 = (I+N)^-1 = I + sum_{i=1}^n N^i, where n is the nilpotent index

	// initialize
	unsigned int dim = a->size1();
	for(unsigned int i=0; i<dim; i++)
		assert((*a)(i,i) == 1);

	// get the nilpotent part
	boost::numeric::ublas::triangular_matrix<mpfr::mpreal, boost::numeric::ublas::upper> N(dim, dim), Nk(dim, dim);
	for(unsigned int i=0; i<dim; i++)
		for(unsigned int j=i; j<dim; j++)
		{
			if(i!=j)
			{
				N(i,j) = (*a)(i,j);
				Nk(i,j) = (*a)(i,j);
			}
		}

	for(unsigned int k=1; k<dim; k++)
	{
		// compute N^k
		if(k>1)
			Nk = prod(N, Nk);

		for(unsigned int i=0; i<dim; i++)
			for(unsigned int j=i; j<dim; j++)
				if((k%2)==0)
					(*b)(i,j) += Nk(i,j);
				else
					(*b)(i,j) -= Nk(i,j);
	}

	// add unit matrix again
	for(unsigned int i=0; i<dim; i++)
		(*b)(i,i) = 1;
}

/*!
 *   @brief Computes the inverse of a given matrix.
 *
 *   @param[in] a A constant pointer to a constant matrix of mpreals specifying the matrix to be inverted.
 *   @param[out] b A constant pointer to a matrix of mpreals. This will be the inverted matrix when the function returns.
 *
 *   @return True if and only if the input matrix was regular.
 *
 *   This algorithm uses LU factorization to compute the inverse of the given matrix.
 */
bool LALA::inverse(const boost::numeric::ublas::matrix<mpfr::mpreal>* const a, boost::numeric::ublas::matrix<mpfr::mpreal>* const b)
{
	// create a working copy of the input
	boost::numeric::ublas::matrix<mpfr::mpreal> A(*a);

	// create a permutation matrix for the LU-factorization
	boost::numeric::ublas::permutation_matrix<std::size_t> pm(A.size1());

	// perform LU-factorization
	if(boost::numeric::ublas::lu_factorize(A, pm) != 0)
		return false;

	// create identity matrix
	b->assign(boost::numeric::ublas::identity_matrix<mpfr::mpreal> (A.size1()));

	// backsubstitute to get the inverse
	boost::numeric::ublas::lu_substitute(A, pm, *b);

	return true;
}

/*!
 *   @brief Computes the inverse of a given matrix.
 *
 *   @param[in] a A constant pointer to a constant matrix of multi-precision rationals specifying the matrix to be inverted.
 *   @param[out] b A constant pointer to a matrix of multi-precision rationals. This will be the inverted matrix when the function returns.
 *
 *   @return True if and only if the input matrix was regular.
 *
 *   This algorithm uses LU factorization to compute the inverse of the given matrix.
 */
bool LALA::inverse(const boost::numeric::ublas::matrix<mpq_class>* const a, boost::numeric::ublas::matrix<mpq_class>* const b)
{
	// create a working copy of the input
	boost::numeric::ublas::matrix<mpq_class> A(*a);

	// create a permutation matrix for the LU-factorization
	boost::numeric::ublas::permutation_matrix<std::size_t> pm(A.size1());

	// perform LU-factorization
	if(boost::numeric::ublas::lu_factorize(A, pm) != 0)
		return false;

	// create identity matrix
	b->assign(boost::numeric::ublas::identity_matrix<mpq_class> (A.size1()));

	// backsubstitute to get the inverse
	boost::numeric::ublas::lu_substitute(A, pm, *b);

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Determinant ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the determinant of a given matrix.
 *
 *   @param[in] m A constant pointer to a constant matrix of mpreals specifying the matrix whose determinant to calculate.
 *
 *   @return The determinant of the matrix stored in an mpreal.
 *
 *   This algorithm uses LU factorization to compute the determinant of the given matrix.
 */
mpfr::mpreal LALA::det(const boost::numeric::ublas::matrix<mpfr::mpreal>* const m)
{
    // create a working copy of the input
    boost::numeric::ublas::matrix<mpfr::mpreal> copyM((*m));
    boost::numeric::ublas::permutation_matrix<std::size_t> pivots(copyM.size1());

    if(boost::numeric::ublas::lu_factorize(copyM, pivots))
    	return 0;

    mpfr::mpreal det = 1;

    for(std::size_t i = 0; i < pivots.size(); i++)
    {
        if(pivots(i) != i)						// row change -> change sign of det
            det *= -1;

        det *= copyM(i, i);
    }
    return det;
}

/*!
 *   @brief Computes the determinant of a given matrix.
 *
 *   @param[in] m A constant pointer to a constant symmetric matrix of mpreals specifying the matrix whose determinant to calculate.
 *
 *   @return The determinant of the matrix stored in an mpreal.
 *
 *   This algorithm uses LU factorization to compute the determinant of the given matrix.
 */
mpfr::mpreal LALA::det(const boost::numeric::ublas::symmetric_matrix<mpfr::mpreal>* const m)
{
	// create working copy
	boost::numeric::ublas::matrix<mpfr::mpreal> copyM((*m));
    return LALA::det(&copyM);
}

/*!
 *   @brief Computes the determinant of a given matrix.
 *
 *   @param[in] m A constant pointer to a constant matrix of multi-precision integers specifying the matrix whose determinant to calculate.
 *
 *   @return The determinant of the matrix stored in a multi-precision integer.
 *
 *   This algorithm uses Gauss-Bareiss to compute the determinant.
 */
mpz_class LALA::detGaussBareiss(const boost::numeric::ublas::matrix<mpz_class>* const m)
{
	assert(m->size1() == m->size2());
	unsigned int dim = m->size1();

	// crate working copy
	boost::numeric::ublas::matrix<mpz_class> copyM((*m));

	// initialize
	int s = 1;
	mpz_class p, c = 1;
	for(unsigned int k=0; k<dim-1; k++)
	{
		p = copyM(k,k);

		// is p = 0?
		if(p==0)
		{
			int nonZeroEntryPos = k;
			bool found = false;

			// get k-th column
			boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<mpz_class> > col(copyM,k);
			for(typename boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<mpz_class> >::iterator it = col.begin()+k; it != col.end(); it++, nonZeroEntryPos++)
			{
				if(*it != 0)
				{
					found = true;
					break;
				}
			}
			if(!found)
				return 0;

			// now for j = k, ...,n exchange the i,j-th entry with the k,j-entry
			for(unsigned int j=k; j<dim; j++)
				boost::swap(copyM(nonZeroEntryPos,j), copyM(k,j));

			// now set s=-s and p=k,k-th entry
			s = -s;
			p = copyM(k,k);
		}

		// step 4: main step - p is now non-zero
		for(unsigned int i=k+1; i<dim; i++)
		{
			for(unsigned int j=k+1; j<dim; j++)
			{
				mpz_class t = p * copyM(i,j) - copyM(i,k)*copyM(k,j);
				copyM(i,j) = t / c;
			}
		}

		c=p;
	}

	return s * copyM(dim-1,dim-1);
}

/*!
 *   @brief Computes the determinant of a given matrix.
 *
 *   @param[in] m A constant pointer to a constant symmetric matrix of multi-precision integers specifying the matrix whose determinant to calculate.
 *
 *   @return The determinant of the matrix stored in a multi-precision integer.
 *
 *   This algorithm uses Gauss-Bareiss to compute the determinant.
 */
mpz_class LALA::detGaussBareiss(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const m)
{
	assert(m->size1() == m->size2());

	// copy symmetric matrix into matrix and call detGaussBareiss
	boost::numeric::ublas::matrix<mpz_class> copyM(m->size1(), m->size2());
	for(unsigned int i=0; i<m->size1(); i++)
		for(unsigned int j=0; j<m->size2(); j++)
			copyM(i,j) = (*m)(i,j);

	return LALA::detGaussBareiss(&copyM);
}

}
