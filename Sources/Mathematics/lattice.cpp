/************************************************************************************//**
 * @author	Gilles Bellot
 * @file	lattice.cpp
 * @date	09/01/2017 - Dortmund - Germany
 *
 * @brief	Implementation of the algorithms to work with lattices over \f$\mathbb{Z}\f$.
 *
 * @version	1.0.2.0
 * @bug 	No known bugs.
 *
 * @copyright	Gilles Bellot @ TU Dortmund
 ****************************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////////////////

// c++ includes
#include <set>

// bell0bytes mathematics
#include "../../Headers/Mathematics/linearAlgebra.h"
#include "../../Headers/Utilities/serviceLocator.h"
#include "../../Headers/Mathematics/matrix.h"
#include "../../Headers/Mathematics/lattice.h"

// FUNCTIONS ////////////////////////////////////////////////////////////////////////////
namespace mathematics
{

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Constructor ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Constructs a lattice from a given Gram matrix.
 *
 *   @param[in] gram A constant pointer to a constant symmetric matrix of multi-precision integers.
 *
 *   For now this constructor stores the Gram matrix, computes the determinant and the dual of the lattice.<br>
 *   The short vectors up to the highest diagonal entry of the Gram matrix are computed.<br>
 *   The successive minima are computed.<br>
 *   A decomposition is computed.<br>
 *
 *   @exception Throws a std::runtime_error if the computation of the short vectors failed.
 */
Lattice::Lattice(const boost::numeric::ublas::symmetric_matrix<mpz_class>* const gram) : gram((*gram))
{
	// LLL reduce the gramian
	LALA::lllGram(&this->gram);

	// compute determinant
	this->determinant = LALA::detGaussBareiss(&this->gram);

	// compute the dual of the lattice
	this->dual();

	// get highest entry B on the diagonal and compute short vectors x with N(x) <= B
	mpz_class highestDiagonalEntry = 0;
	for(unsigned int i=0; i<this->gram.size1(); i++)
		if(highestDiagonalEntry < this->gram(i,i))
			highestDiagonalEntry = this->gram(i,i);

	try{ nShVecs = LALA::shortVectors(&this->gram, highestDiagonalEntry, &shVecs).get(); }
	catch(std::runtime_error& e)
	{
		throw e;
	}

	// compute decomposition
	this->constructDecomposition();

	// compute sucMin
	this->successiveMinima();

	// compute automorphism group
}

/*!
 *   @brief Copy constructor.
 *
 *   @param[in] L A reference to a constant lattice to copy from.
 *
 */
Lattice::Lattice(const Lattice& L)
{
	this->gram = L.gram;
	this->dualGram = L.dualGram;
	this->nShVecs = L.nShVecs;
	this->shVecs = L.shVecs;
	this->determinant = L.determinant;
	this->sucMin = L.sucMin;
}

/*!
 *   @brief The destructor destroys!
 */
Lattice::~Lattice()
{
	this->gram.clear();
	this->dualGram.clear();
	this->nShVecs = 0;
	this->shVecs.clear();
	this->determinant = 0;
}


// private functions
/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Dual ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the dual of the lattice.
 *
 *   @return void
 */
void Lattice::dual()
{
	// copy gram into mpq matrix
	unsigned int dim = gram.size1();
	dualGram.resize(dim, dim, false);
	boost::numeric::ublas::matrix<mpq_class> gramCopy(dim,dim), dualGramCopy(dim, dim);
	for(unsigned int i=0; i<dim; i++)
		for(unsigned int j=0; j<dim; j++)
			gramCopy(i,j) = gram(i,j);

	LALA::inverse(&gramCopy, &dualGramCopy);

	for(unsigned int i=0; i<dim; i++)
		for(unsigned int j=0; j<dim; j++)
			dualGram(i,j) = this->determinant * dualGramCopy(i,j);

	LALA::lllGram(&dualGram);
}

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Successive Minima ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the successive minima of the lattice.
 *
 *   @todo Clean up.
 *
 *   @return void
 *
 *   This function computes the successive minima of the lattice, see [HV2] for further details.
 */
void Lattice::successiveMinima()
{
	unsigned int dim = this->gram.size1();
	this->sucMin.resize(dim);

	// partial basis matrix - with basis extension
	boost::numeric::ublas::matrix<mpq_class> invBasisMatrix(dim,dim);

	// copy short vectors in new structure in reverse order
	std::vector<std::pair<mpz_class, std::valarray<mpz_class> > > S(nShVecs);
	std::vector<std::valarray<mpz_class> > U;
	for(unsigned int i=0; i<nShVecs; i++)
		S[i] = shVecs[nShVecs-1-i];

	// initialization with minimal v \in S
	unsigned int n = 0;
	std::valarray<mpz_class> v = S[nShVecs-1].second;
	mpz_class nV = 0;
	sucMin[n] = S[nShVecs-1].first;
	U.push_back(v);
	S.pop_back();

	// set L = Zv , extend to vector space basis to test for inclusion
	for(unsigned int i=0; i<dim; i++)
		invBasisMatrix(i,0) = v[i];

	// extend to vector space basis
	LALA::basisExtensionZ(&invBasisMatrix);

	// invert matrix
	LALA::inverse(&invBasisMatrix, &invBasisMatrix);

	while(!S.empty())
	{
		// chose v with minimal norm and set S = S \ {v}
		v = S[S.size()-1].second;
		nV = S[S.size()-1].first;
		S.pop_back();

		// v in L?
		if(!LALA::elementOfLattice(&invBasisMatrix, &v, U.size()))
		{
			// test linear dependency over R
			U.push_back(v);
			if(!LALA::linearlyIndependentR(&U))
				U.pop_back();
			else
			{
					n++;
					sucMin[n] = nV;
					if(n==sucMin.size()-1)
						break;
			}

			// set L = \sum Zv, v \in sucmin, extend to vector space basis to test for inclusion
			invBasisMatrix.clear();
			for(unsigned int j=0; j<U.size(); j++)
				for(unsigned int i=0; i<dim; i++)
					invBasisMatrix(i,j) = U[j][i];

			LALA::basisExtensionZ(&invBasisMatrix);

			// invert matrix
			LALA::inverse(&invBasisMatrix, &invBasisMatrix);
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Decomposition //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Computes the decomposition of the lattice.
 *
 *   See [HB] chapter 3.
 *
 *   @todo Clean up.
 *
 *   @return void
 */
void Lattice::constructDecomposition()
{
	// see [HB] chapter 3

	// decomposition
	std::vector<std::vector<std::valarray<mpz_class> > > dc, tempDC;
	std::vector<std::valarray<mpz_class> > tildeL;

	// step 1: init
	unsigned int dim = this->gram.size1();

#ifndef NDEBUG
	unsigned int verboseLevel = util::ServiceLocator::getProgOpts()->verboseLevel;
#endif

	boost::numeric::ublas::vector<mpq_class> prod;

	// save short vectors in order of decreasing norms
	std::vector<std::valarray<mpz_class> > sv;
	for(unsigned int i=0; i<this->nShVecs; i++)
		sv.push_back(this->shVecs[this->nShVecs-1-i].second);

	// chose v \in shVecs with N(v) minimal
	std::valarray<mpz_class> v = sv[sv.size()-1];
	// delete v from shVecs
	sv.pop_back();

	// set L_0 to the Z-span of v
	std::vector<std::valarray<mpz_class> > l0;
	l0.push_back(v);
	dc.push_back(l0);

#ifndef NDEBUG
	if(verboseLevel > 9)
	{
		std::cout << "first basis vector: ";
		for(unsigned int i=0; i<dim; i++)
			std::cout << v[i] << " ";
		std::cout << std::endl << std::endl;
	}
#endif

	// keep track of non orthogonal projections
	std::set<unsigned int> nonZeroOrthProj;

	// vector space basis (orthogonal complement)
	unsigned int basisDimension = 1;
	std::vector<std::valarray<mpq_class> > latticeBasis, vectorSpaceComplement;
	std::valarray<mpq_class> vq(dim);
	for(unsigned int i=0; i<dim; i++)
		vq[i] = v[i];
	latticeBasis.push_back(vq);

	LALA::orthogonalComplement(&latticeBasis, &vectorSpaceComplement, &this->gram);

	boost::numeric::ublas::matrix<mpq_class> invBasisMatrix(dim,dim);

	for(unsigned int i=0; i<latticeBasis.size(); i++)
		for(unsigned int j=0; j<dim; j++)
			invBasisMatrix(j,i) = latticeBasis[i][j];

	for(unsigned int i=0; i<vectorSpaceComplement.size(); i++)
		for(unsigned int j=0; j<dim; j++)
			invBasisMatrix(j,i+latticeBasis.size()) = vectorSpaceComplement[i][j];

	// invert
	LALA::inverse(&invBasisMatrix, &invBasisMatrix);

	// successively add vectors from shVecs to the L_i
	while(!sv.empty())
	{
		// chose v \in shVecs with N(v) minimal
		v = sv[sv.size()-1];
		// delete v from shVecs
		sv.pop_back();

#ifndef NDEBUG
	if(verboseLevel > 9)
	{
		std::cout << "Testing vector v: ";
		for(unsigned int i=0; i<dim; i++)
			std::cout << v[i] << " ";
		std::cout << "\t ... \t";

		bool test = LALA::elementOfLattice(&invBasisMatrix, &v, basisDimension, &prod);
		if(test)
			std::cout << " in the lattice \t ... \t discard!\n";
		else
			std::cout << " not in the lattice \t ... \t continue!\n\n";
	}
#endif

		if(!LALA::elementOfLattice(&invBasisMatrix, &v, basisDimension, &prod))
		{
			// if v is not in L_0 + ... + L_k, then v is indecomposable

			// compute \pi_i(v), i.e. the orthogonal projection of v into the i-th indecomposable component of the lattice
			// this can be read from invBasisMatrix * v - see [HB] for details
			for(unsigned int i=0; i<basisDimension; i++)
				if(prod[i].get_den() != 1)
					nonZeroOrthProj.insert(i);

#ifndef NDEBUG
	if(verboseLevel > 9)
	{
		std::cout << "Non orthogonal projections: ";
		for(std::set<unsigned int>::iterator it = nonZeroOrthProj.begin(); it != nonZeroOrthProj.end(); it++)
			std::cout << (*it) << " ";
		std::cout << std::endl << std::endl;
	}
#endif

			// compute basis for L = Zv + sum L_i ; i \in nonZeroOrthProj
			basisDimension = 0;
			tildeL.clear();
			for(struct {std::set<unsigned int>::iterator it;  unsigned int i;} s = {nonZeroOrthProj.begin(), 0}; s.it != nonZeroOrthProj.end(); s.it++)
			{
				while(basisDimension < (*s.it+1))
				{
					basisDimension += dc[s.i].size();
					s.i++;
				}

				// add basis of L_i to tildeL
				for(unsigned int j=0; j<dc[s.i-1].size(); j++)
				{
					std::valarray<mpz_class> tVec(dim);
					for(unsigned int k=0; k<dim; k++)
						tVec[k] = dc[s.i-1][j][k];
					tildeL.push_back(tVec);
				}
			}

			basisDimension = tildeL.size();

			// copy tildeL into matrix to do LLL
			boost::numeric::ublas::matrix<mpz_class> basisMatrix(dim,basisDimension+1);
			for(unsigned int i=0; i<basisDimension; i++)
				for(unsigned int j=0; j<dim; j++)
					basisMatrix(j, i) = tildeL[i][j];

			// add v
			for(unsigned int s=0; s<dim; s++)
				basisMatrix(s,basisDimension) = v[s];
			basisDimension++;

			LALA::mlllBasis(&basisMatrix);

#ifndef NDEBUG
	if(verboseLevel > 9)
	{
		std::cout << "tildeL:\n";
		MatrixHelper<mpz_class>::printMatrix(&basisMatrix);
		std::cout << "\n\n";
	}
#endif

			// copy basis back to vector system
			tildeL.clear();
			for(unsigned int i=0; i<basisDimension; i++)
			{
				std::valarray<mpz_class> tildeV(dim);
				for(unsigned int j=0; j<dim; j++)
					tildeV[j] = basisMatrix(j,i);
				tildeL.push_back(tildeV);
			}

			// re-arrange the list of the L_i
			if(nonZeroOrthProj.empty())
			{
				// if the vector was orthogonal to all the others, add new sublattice spanned by this vector
				dc.push_back(tildeL);
			}
			else
			{
				// else rearrange the lattices
				tempDC = dc;
				dc.clear();
				for(unsigned int i=0, j=i; i<tempDC.size(); i++)
				{
					i==0 ? j=i : j+=tempDC[i-1].size();
					if(nonZeroOrthProj.find(j) == nonZeroOrthProj.end())
						dc.push_back(tempDC[i]);
				}
				dc.push_back(tildeL);
			}

#ifndef NDEBUG
	if(verboseLevel > 9)
	{
		std::cout << "decomp:\n";
		for(unsigned int i=0; i<dc.size(); i++)
		{
			std::cout << "\t#" << i+1 << ": ";
			for(unsigned int j=0; j<dc[i].size(); j++)
			{
				for(unsigned int k=0; k<dim; k++)
					std::cout << dc[i][j][k] << " ";
				if(j != dc[i].size()-1)
					std::cout << " ; ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
#endif

			// prepare to test next vector
			nonZeroOrthProj.clear();

			// get inverse of new basis
			latticeBasis.clear();
			vectorSpaceComplement.clear();
			basisDimension = 0;
			for(unsigned int i=0; i<dc.size(); i++)
			{
				for(unsigned int j=0; j<dc[i].size(); j++)
				{
					std::valarray<mpq_class> vq(dim);
					for(unsigned int s=0; s<dim; s++)
						vq[s] = dc[i][j][s];
					latticeBasis.push_back(vq);
				}
				basisDimension += dc[i].size();
			}

			if(basisDimension == dim)
				break;

			// get orthogonal complement
			LALA::orthogonalComplement(&latticeBasis, &vectorSpaceComplement, &this->gram);

			invBasisMatrix.clear();
			for(unsigned int i=0; i<latticeBasis.size(); i++)
				for(unsigned int j=0; j<dim; j++)
					invBasisMatrix(j,i) = latticeBasis[i][j];

			for(unsigned int i=0; i<vectorSpaceComplement.size(); i++)
				for(unsigned int j=0; j<dim; j++)
					invBasisMatrix(j,i+latticeBasis.size()) = vectorSpaceComplement[i][j];

			// invert
			LALA::inverse(&invBasisMatrix, &invBasisMatrix);
		}
	} // end while
	decomp = dc;
}

// public functions

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Element Of /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Tests whether a given vector is an element of the lattice or not.
 *
 *   @param v A constant pointer to a constant valarray of multi-precision integers specifying the vector to test.
 *
 *   @return True if and only if \f$G^{-1} \cdot v\f$ is integral or not.
 */
bool Lattice::contains(const std::valarray<mpz_class>* const v) const
{
	// check whether G^-1 * v is integral or not

	unsigned int dim = v->size();

	// copy v into a boost vector
	boost::numeric::ublas::vector<mpz_class> boostV(dim);
	boost::numeric::ublas::vector<mpq_class> boostRes(dim);
	for(unsigned int i=0; i<dim; i++)
		boostV(i) = (*v)[i];

	boostRes = boost::numeric::ublas::prod(this->dualGram, boostV);

	for(unsigned int i=0; i<dim; i++)
	{
		mpq_class t = boostRes[i] / this->determinant;
		if(t.get_den() != 1)
			return false;
	}

	return true;
}

/*!
 *   @brief Search for a specific vector in the list of suitable short vectors of L.
 *
 *   @param v A constant pointer to a valarray of multi-precision integers specifying the vector to find.
 *   @param indicesOfSuitableShortVectors The integers in this vector specify the indices of the suitable short vectors in shVecs.
 *
 *   @return The index of v in shVecs, or a std::runtime_exception if the vector was not found.
 *
 *  The return value \f$i\f$ is positive, if \f$v\f$> is the \f$i\f$-th short vector of L.
 *  The return value is negative, if it is the negative of that vector; in that case, v is set to -v as well.
 *  <br>Uses quicksort.
 */
util::Expected<int> Lattice::find(boost::numeric::ublas::vector<mpz_class>* const v, const std::vector<unsigned int>* const indicesOfSuitableShortVectors) const
{
		unsigned int dim = v->size(), low=0, high=indicesOfSuitableShortVectors->size()-1, search=0, i=0;
		int sign = 1, cmp;

		for(i=0; i<dim && (*v)(i) == 0; i++);

		if(i < dim && (*v)(i) < 0)
		{
			sign = -1;
			for(i=0; i<dim; i++)
				(*v)(i) *= -1;
		}
		while(low <= high)
		{
			search = low + (high-low)/2;
			cmp = LALA::compare(v, &shVecs[indicesOfSuitableShortVectors->at(search)].second);
			if(cmp == 1)
				// v is in the upper half
				low = search + 1;
			else if(cmp == -1)
				// v is in the lower half
				high = search - 1;
			else
				// found v
				return sign * search;
		}

		// if low > high, v was not found
		return sign * indicesOfSuitableShortVectors->size()+low;
}
/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Print //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/*!
 *   @brief Prints information about the lattice.
 *
 *   @param flags Those flags specify how much information to print.
 *
 *   @return void
 *
 *   @section Flags
 *   - 0: print only the Gram matrix and the determinant
 *   - 1: in addition, print the number of short vectors
 *   - 2: in addition, print the actual short vectors, the successive minima and the decomposition
 *   - 3: in addition, print the dual and its determinant
 */
void Lattice::print(const unsigned int flags) const
{
	unsigned int n = this->gram.size1();

	// print gramian
	std::cout << "Lattice with gramian:\n";
	MatrixHelper<mpz_class>::printMatrix(&this->gram);
	// print determinant
	std::cout << "Determinant: " << this->determinant << std::endl << std::endl;

	if(flags >= 1)
	{
		// print short vectors
		std::cout << "Number of short vectors: " << nShVecs << std::endl << std::endl;

		if(flags >= 2)
		{
			// print the actual short vectors
			std::cout << "\tShort vectors:\n";
			for(unsigned long i=0; i<nShVecs; i++)
			{
				std::cout << "\t\t#" << i+1 << ": ";
				for(unsigned int j=0; j<n; j++)
					std::cout << shVecs[i].second[j] << " ";
				std:: cout << " (" << shVecs[i].first << ")\n";
			}
			std::cout << std::endl;

			// print the successive minima
			std::cout << "\tSuccessive Minima: ";
			for(unsigned int i=0; i<sucMin.size(); i++)
				std::cout << sucMin[i] << " ";
			std::cout << std::endl << std::endl;

			// print decomposition
			std::cout << "\tOrthogonal Decomposition:\n";
			for(unsigned int i=0; i<decomp.size(); i++)
			{
				std::cout << "\t\t#" << i+1 << ": ";
				for(unsigned int j=0; j<decomp[i].size(); j++)
				{
					for(unsigned int k=0; k<decomp[i][j].size(); k++)
						std::cout << decomp[i][j][k] << " ";
					if(j != decomp[i].size()-1)
						std::cout << " ; ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl << std::endl;
		}

		if(flags >= 3)
		{
			// print dual
			std::cout << "Dual gramian:\n";
			MatrixHelper<mpz_class>::printMatrix(&dualGram);
			std::cout << "Dual Gramian denominator: " << this->determinant << std::endl << std::endl;
		}
	}
}

}
