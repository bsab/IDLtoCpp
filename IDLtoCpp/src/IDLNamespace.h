/*
 * IDLMethod.h
 *
 *  Created on: Jul 3, 2013
 *      Author: remocean
 */

#ifndef IDLMETHOD_H_
#define IDLMETHOD_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <cmath>        // std::abs
#include <complex>
#include "ExceptionDef.h"

#include <iostream>
#include <fstream>

namespace bst = boost::numeric::ublas;
using namespace bst;
using namespace std;

namespace idl {
typedef unsigned size_type;

template<class T>
class Arr: public matrix<T> {
public:

	/**
	 * Default dense matrix constructor. Make a dense matrix of size (0,0)
	 */
	BOOST_UBLAS_INLINE
	Arr<T>() :
			matrix<T>() {
		_nRows = _nCols = 0;
	}

	/**
	 * (serve per fare la transpose)
	 */
	BOOST_UBLAS_INLINE
	Arr<T>(matrix<T> &e) :
			matrix<T>(e) {
		_nRows = e.size1();
		_nCols = e.size2();
	}

	/**
	 * Default dense matrix constructor. Make a dense matrix of size (1,dim)
	 * (vettore 1D)
	 */
	BOOST_UBLAS_INLINE
	Arr<T>(const int dim) :
			matrix<T>(1, dim) {
		_nRows = 1;
		_nCols = dim;
		setToZero();
	}

	/**
	 * Default dense matrix constructor. Make a dense matrix of size (1,dim)
	 * (vettore 2D)
	 */
	BOOST_UBLAS_INLINE
	Arr<T>(const int &s1, const int &s2) :
			matrix<T>(s1, s2) {
		_nRows = s1;
		_nCols = s2;
		setToZero();
	}

	// Accessors
	/** Return the number of colums of the matrix (vettore 1D)
	 */
	BOOST_UBLAS_INLINE
	size_type size() const {
		//return matrix<T>::size2();
		return _nCols;
	}

	BOOST_UBLAS_INLINE
	size_type numRows() const {
		return _nRows;
	}

	BOOST_UBLAS_INLINE
	size_type numCols() const {
		return _nCols;
	}

	/** Azzera la matrice
	 */
	BOOST_UBLAS_INLINE
	void setToZero() const {
		this->clear();
	}

	// Resizing
	/** Resize a matrix to new dimensions
	 * \param size the new number of colums
	 */
	BOOST_UBLAS_INLINE
	void resize(size_type s) {
		//size_type s1 = matrix<T>::size1();
		matrix<T>::resize(_nRows, s);
		_nCols = s;
	}

	/** Resize a matrix to new dimensions
	 * \param s1 the new number of rows
	 * \param s2 the new number of colums
	 */
	BOOST_UBLAS_INLINE
	void resize(size_type s1, size_type s2) {
		matrix<T>::resize(s1, s2);
		_nRows = s1;
		_nCols = s2;
	}

	/**
	 * ATTENZIONE: vale solo per i vettori 1D
	 */
	BOOST_UBLAS_INLINE
	void add(const T &v) {
		//size_type s = size();
		_nCols++;
		resize(_nCols);
		operator ()(_nCols - 1) = v;
	}

	/**
	 *
	 */
	BOOST_UBLAS_INLINE
	//Arr<T>& assign(Arr<T> &A) {
	Arr<T>& assign(T a) {
		for (size_type x = 0; x < this->size1(); ++x) {
			for (size_type y = 0; y < this->size2(); ++y) {
				operator ()(x, y) = a;
			}
		}

		return *this;
	}

	/**
	 * Per l'accesso singolo come vettore 1D
	 */
	BOOST_UBLAS_INLINE
	T & operator ()(size_type i) const {
		return this->at_element(0, i);
	}

	/**
	 * Per l'accesso singolo come vettore 2D
	 */
	BOOST_UBLAS_INLINE
	T & operator ()(size_type i, size_type j) const {
		return this->at_element(i, j);
	}

	/**
	 *
	 */
	BOOST_UBLAS_INLINE
	Arr<T>& operator *=(const T &at) {
		matrix_assign_scalar<scalar_multiplies_assign>(*this, at);

		return *this;
	}

private:
	int _nRows;
	int _nCols;
}
;

// ****************** END CLASS ARR *************************

// ********************* TIPI IDL *************************
typedef char BYTE; 	// An 8-bit unsigned integer ranging in value from 0 to 255.
typedef int INT;	// A 16-bit signed integer ranging from -32,768 to +32,767.
typedef unsigned int UINT;// A 16-bit unsigned integer ranging from 0 to 65535.
typedef long int LONG;					// A 32-bit signed integer
typedef unsigned long int ULONG; 		// A 32-bit unsigned integer
typedef long long int LONG64;			// A 64-bit signed integer
typedef unsigned long long int ULONG64; // A 64-bit unsigned integer
typedef float FLOAT; 		// A 32-bit, single-precision, floating-point number
typedef double DOUBLE; 		// A 64-bit, double-precision, floating-point number
typedef std::complex<float> COMPLEX; // A real-imaginary pair of single-precision, floating-point numbers.
typedef std::complex<double> DCOMPLEX; // A real-imaginary pair of double-precision, floating-point numbers.
typedef string STRING;

// ********************* END TIPI IDL *************************

// ********************* METODI IDL *************************

/**
 *
 */
template<class T>
static T sum_elem(const Arr<T> &mat) {
	T sum = 0;
	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			sum += mat(x, y);
		}
	}

	return sum;
}

/**
 *
 * @param dim
 * @return
 */
static Arr<idl::FLOAT> findgen(size_type dim) {
	Arr<idl::FLOAT> a(dim);
	idl::FLOAT value = 0;
	for (size_type i = 0; i < dim; ++i) {
		a(i) = value++;
	}

	return a;
}

/**
 *
 */
template<class T>
static T max(const Arr<T> &mat) {
	T maximum = 0;

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			if (mat(x, y) > maximum)
				maximum = mat(x, y);
		}
	}

	return maximum;
}

/**
 *
 */
template<class T>
static T min(const Arr<T> &mat) {
	T min = numeric_limits<T>::max();

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			if (mat(x, y) < min)
				min = mat(x, y);
		}
	}

	return min;
}

/**
 *
 */
template<class T, class E>
Arr<T> cpow(const Arr<T> &mat, const E exp) {
	Arr<T> tmat(mat);

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			tmat(x, y) = pow(mat(x, y), exp);
		}
	}

	return tmat;
}

/**
 *
 */
template<class T>
Arr<T> csqrt(const Arr<T> &mat) {
	Arr<T> tmat(mat);

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			tmat(x, y) = sqrt(mat(x, y));
		}
	}

	return tmat;
}

/**
 * Esegue la ricerca del minimo, scartando i primi "scrap" elementi
 * e ritorna la prima coppia di indici che contiene il minimo.
 */
template<class T, class T1, class T2>
static T whereMin(const Arr<T> &mat, const int scrap,
		std::pair<T1, T2> &firstIndex) {
	T min = numeric_limits<T>::max();

	int stop = 0;
	//cout << "scrap: " << scrap << endl;
	int lim1 = mat.size1();	//nota: devono essere int perchè il for deve fare check sul segno
	int lim2 = mat.size2();

	for (int x = lim1 - 1; x >= 0; --x) {
		if (x == 0)
			stop = scrap;
		for (int y = lim2 - 1; y >= stop; --y) {
			//cout << "( " << x << "," << y << ")" << endl;
			if (mat(x, y) < min) {
				min = mat(x, y);
				firstIndex.first = x;
				firstIndex.second = y;
			}
		}
	}

	return min;
}

/**
 * Esegue la ricerca del massimo e ritorna
 * l'ultima coppia di indici che lo contiene.
 */
template<class T, class T1, class T2>
static T whereMax(const Arr<T> &mat, std::pair<T1, T2> &lastIndex) {
	T max = numeric_limits<T>::min();

	int lim1 = mat.size1();	//nota: devono essere int perchè il for deve fare check sul segno
	int lim2 = mat.size2();

	for (int x = 0; x < lim1; ++x) {
		for (int y = 0; y < lim2; ++y) {
			//cout << "( " << x << "," << y << ")" << endl;
			if (mat(x, y) > max) {
				max = mat(x, y);
				lastIndex.first = x;
				lastIndex.second = y;
			}
		}
	}

	return max;
}

/**
 * Resizes a vector or array to dimensions given by the parameters d.
 * The supplied dimensions must be integral multiples or factors of the original dimension.
 * The expansion or compression of each dimension is independent of the others,
 * so that each dimension can be expanded or compressed by a different value.
 *
 */
template<class T>
static Arr<T> rebin(const Arr<T> &a, const LONG dim) throw (int) {
	size_type size = a.size();
	if (((size % dim) != 0) && ((dim % size) != 0))
		throw EXCP_REBIN_DIM_NOT_VALID;

	Arr<T> b(dim);
	size_type step = size % dim;

	size_type pos = 0;
	for (size_type i = 0; i < size; ++i) {
		for (size_type k = 0; k < step; ++k) {
			b(pos + k) = a(i);
		}
		pos += step;
	}

	return b;
}

/**
 *
 * @param a
 * @param dim
 * @return
 *
 * for each row
 *   rebin (col[0]...col[size2-1], dim2
 *
 *   TODO da rivedere meglio, perchè in IDL calcola l'interpolazione anche
 *   (per esempio) fra a[0][3] e a[1][0] (a(n,4)
 */
template<class T, class E>
static Arr<T> rebin(const Arr<T> &a, E dim1, E dim2) throw (int) {
	//if (((size1 % dim1) != 0) && ((dim1 % size1) != 0))
	//	throw EXCP_REBIN_DIM_NOT_VALID;

	int size2 = (int) a.size2();
	if (((size2 % (int) dim2) != 0) && (((int) dim2 % size2) != 0))
		throw EXCP_REBIN_DIM_NOT_VALID;

	Arr<T> b(a);
	b.resize(dim1, dim2);

	size_type rpos = 0;
	size_type pos = 0;
	/*
	 for (size_type row = 0; row < size1; row++) {
	 T f = dim2 / size2; //calcolo il fattore di interpolazione
	 T interp = 0; //interpolante
	 pos = 0;
	 for (size_type k = 0; k < size2 - 1; ++k) {
	 interp = (a(row, k + 1) - a(row, k)) / f;
	 //for (int i = 0; i < dim; ++i) {
	 for (int j = 0; j < (int) f; ++j) {
	 if (j == dim2)
	 break;
	 b(row, pos) = a(row, k) + (j * interp);
	 //cout << " (row=" << row << " - pos=" << pos << "): "
	 //		<< b(row, pos) << endl;
	 pos++;
	 }
	 //}
	 }
	 // The last m/n points of the result are obtained by duplicating element n-1
	 //of the input array because their interpolates would lie outside the input array.
	 rpos++;
	 }
	 */
	rpos = a.size1();
	//copio le righe restanti
	for (size_type i = rpos; i < dim1; ++i) {
		for (size_type j = 0; j < dim2; ++j) {
			b(i, j) = b(i - 1, j);
		}
	}

	pos = a.size2();
	//copio le colonne restanti
	for (size_type i = 0; i < dim1; ++i) {
		for (size_type j = pos; j < dim2; ++j) {
			b(i, j) = b(i, j - 1);
		}
	}

	return b;
}

/**
 *
 * @return Floatarray a una dimensione.
 */
static Arr<FLOAT> fltarr(const LONG dim) {
	return Arr<FLOAT>(dim);
}

/**
 *
 * @return Floatarray a due dimensioni.
 */
static Arr<FLOAT> fltarr(const LONG m, const LONG n) {
	return Arr<FLOAT>(m, n);
}

/**
 *
 * @return Longarray a una dimensione.
 */
static Arr<LONG> lonarr(const LONG dim) {
	return Arr<LONG>(dim);
}

/**
 *
 * @return Longarray a due dimensioni.
 */
static Arr<LONG> lonarr(const LONG m, const LONG n) {
	return Arr<LONG>(m, n);
}

/**
 *
 */
template<class T>
static Arr<T> abs(const Arr<T> &mat) {

	Arr<T> tmat(mat);
	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			tmat(x, y) = std::abs(mat(x, y));
		}
	}

	return tmat;
}

/**
 * Per i tipi complessi, il valore assoluto è definito
 * come la radice quadrata della somma dei quadrati della parte reale
 * e della parte immaginaria.
 */
static Arr<FLOAT> c_abs(const Arr<COMPLEX> &mat) {
	Arr<FLOAT> tmat(mat.numRows(), mat.numCols());

	for (size_type x = 0; x < mat.numRows(); ++x) {
		for (size_type y = 0; y < mat.numCols(); ++y) {
			COMPLEX c = mat(x, y);
			tmat(x, y) = abs(c);
		}
	}

	return tmat;
}

/**
 *
 */
template<class T>
static Arr<T> tanh(const Arr<T> &mat) {
	Arr<T> tmat(mat);

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			tmat(x, y) = std::tanh(mat(x, y));
		}
	}

	return tmat;
}

/**
 *
 */
template<class T>
static Arr<T> sqrt(const Arr<T> &mat) {
	Arr<T> tmat(mat);

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			tmat(x, y) = std::sqrt(mat(x, y));
		}
	}

	return tmat;
}

/**
 * returns the number of elements contained in an expression or variable.
 */
template<class T>
static int n_elements(Arr<T> &mat) {
	return mat.size1() * mat.size2();
}

/**
 *
 */
template<class T>
static Arr<T> transpose(const Arr<T> &mat) {
	size_type s1 = mat.size1();
	size_type s2 = mat.size2();
	Arr<T> tmat(s2, s1);

	for (size_type x = 0; x < s1; ++x) {
		for (size_type y = 0; y < s2; ++y) {
			tmat(y, x) = mat(x, y);
		}
	}

	return tmat;
}

/**
 *
 */
template<class T, class E>
Arr<T> minus(const Arr<T> &A, const E &t) {
	size_type s1 = A.size1();
	size_type s2 = A.size2();
	Arr<T> tmat(s1, s2);

	for (size_type x = 0; x < s1; ++x) {
		for (size_type y = 0; y < s2; ++y) {
			tmat(x, y) = A(x, y) - t;
		}
	}

	return tmat;
}

/**
 *
 */
template<class T, class E>
Arr<T> sum(const Arr<T> &A, const E &t) {
	size_type s1 = A.size1();
	size_type s2 = A.size2();
	Arr<T> tmat(s1, s2);

	for (size_type x = 0; x < s1; ++x) {
		for (size_type y = 0; y < s2; ++y) {
			tmat(x, y) = A(x, y) + t;
		}
	}

	return tmat;
}

/**
 *
 */
template<class T>
Arr<T> sum(const Arr<T> &A, const Arr<T> &B) {
	size_type a1 = A.size1();
	size_type a2 = A.size2();
	size_type b1 = B.size1();
	size_type b2 = B.size2();

	if ((b1 < a1) && (b2 < a2)) {
		cerr << "a1: " << a1 << " a2:" << a2 << " b1:" << b1 << " b2:" << b2
				<< endl;
		throw EXCP_SUM_BOUND_NOT_VALID;
	}
	Arr<T> tmat(a1, a2);

	for (size_type x = 0; x < a1; ++x) {
		for (size_type y = 0; y < a2; ++y) {
			tmat(x, y) = A(x, y) + B(x, y);
		}
	}

	return tmat;
}

/**
 *
 */
template<class T>
Arr<T> minus(const Arr<T> &A, const Arr<T> &B) {
	size_type a1 = A.size1();
	size_type a2 = A.size2();
	size_type b1 = B.size1();
	size_type b2 = B.size2();

	if ((b1 < a1) && (b2 < a2)) {
		cerr << "a1: " << a1 << " a2:" << a2 << " b1:" << b1 << " b2:" << b2
				<< endl;
		throw EXCP_MINUS_BOUND_NOT_VALID;
	}
	Arr<T> tmat(a1, a2);

	for (size_type x = 0; x < a1; ++x) {
		for (size_type y = 0; y < a2; ++y) {
			tmat(x, y) = A(x, y) - B(x, y);
		}
	}

	return tmat;
}

/**
 *
 */
//BOOST_UBLAS_INLINE
template<class T, class E>
Arr<T> dot(const Arr<T> &A, const E &at) {
	Arr<T> tmat(A);

	matrix_assign_scalar<scalar_multiplies_assign>(tmat, at);

	return tmat;
}

/**
 *
 */
//BOOST_UBLAS_INLINE
template<class T>
Arr<T> dot(const Arr<T> &A, const Arr<T> &B) {
	//matrix<T> res = bst::prod(A, B);
	//return res;

	size_type a1 = A.size1();
	size_type a2 = A.size2();
	size_type b1 = B.size1();
	size_type b2 = B.size2();

	if ((b1 < a1) && (b2 < a2)) {
		cerr << "a1: " << a1 << " a2:" << a2 << " b1:" << b1 << " b2:" << b2
				<< endl;
		throw EXCP_MINUS_BOUND_NOT_VALID;
	}
	Arr<T> tmat(a1, a2);

	for (size_type x = 0; x < a1; ++x) {
		for (size_type y = 0; y < a2; ++y) {
			tmat(x, y) = A(x, y) * B(x, y);
		}
	}

	return tmat;
}

/**
 * Torna una sottomatrice specificata dagli indici di range
 * @param A
 * @param r1
 * @param r2
 * @param c1
 * @param c2
 * @return
 */
template<class T>
Arr<T> cut(const Arr<T> &A, const size_type &r1, const size_type &r2,
		const size_type &c1, const size_type &c2) {

	matrix<T> m = A;
	matrix<T> m2 = matrix_range<matrix<T> >(m, range(r1, r2), range(c1, c2));

	//copio il contenuto del range in una matrice
	Arr<T> tmat(r2 - r1, c2 - c1);
	for (size_type i = 0; i < m2.size1(); ++i)
		for (size_type j = 0; j < m2.size2(); ++j)
			tmat(i, j) = m2(i, j);

	return tmat;
}

template<class T>
Arr<T> cut(const Arr<T> &A, const size_type &c1, const size_type &c2) {

	//calcolo il subrange
	//cout << "CUT(" << c1 << " , " << c2 << ")" << endl;

	matrix<T> m2;
	matrix<T> m = A;
	//if (A.size1() > 1)
	//	m2 = matrix_range<matrix<T> >(m, range(c1, c2), range(0, m.size2()));
	//else
	m2 = matrix_range<matrix<T> >(m, range(0, m.size1()), range(c1, c2));

	size_type s1 = m2.size1();
	size_type s2 = m2.size2();
	//se faccio il range su una matrice allora il risultato sarà
	//la matrice trasformata in un vettore
	if (s1 > 1) {
		Arr<T> tmp(m2.size1() * m2.size2());
		size_type i = 0;
		for (size_type x = 0; x < s1; ++x) {
			for (size_type y = 0; y < s2; ++y) {
				tmp(i) = m2(x, y);
				i++;
			}
		}
		return tmp;
	} else
		return Arr<T>(m2);
}

/**
 *
 * @param c1
 * @param c2
 * @return
 */
static bool cgt(const idl::COMPLEX & c1, const idl::COMPLEX & c2) {
	return (c1.real() > c2.real() && c1.imag() > c2.imag());
}

/**
 *
 * @param c1
 * @param c2
 * @return
 */
static bool clt(const idl::COMPLEX & c1, const idl::COMPLEX & c2) {
	return (c1.real() < c2.real() && c1.imag() < c2.imag());
}

/**
 *
 */
static COMPLEX cmax(const Arr<COMPLEX> &mat) {
	COMPLEX maximum(0, 0);

	for (size_type x = 0; x < mat.size1(); ++x) {
		for (size_type y = 0; y < mat.size2(); ++y) {
			COMPLEX elem = mat(x, y);
			if (cgt(elem, maximum))
				maximum = mat(x, y);
		}
	}

	return maximum;
}

/**
 * Ritorna le dimensioni di matrice 3D
 * @param
 * @return
 */
template<class T>
static bst::vector<idl::LONG> getMat3Ddimensions(
		const bst::vector<idl::Arr<T> > & mat3D) {
	bst::vector<idl::LONG> size(3);

	size(0) = mat3D[0].size1(); 	//dimensione X
	size(1) = mat3D[0].size2(); 	//dimensione Y
	size(2) = mat3D.size(); 		//dimensione Z

	return size;
}

/**
 * Riduce il numero di slice (ovvero l'asse Z) di una matrice 3D
 * @param mat3D
 * @param from
 * @param to
 * @return
 */
template<class T>
static bst::vector<idl::Arr<T> > resizeSliceMat3D(
		const bst::vector<idl::Arr<T> > & mat3D, size_type from, size_type to) {
	bst::vector<idl::Arr<T> > tmat(to - from);

	size_type j = 0;
	for (size_type i = from; i < to; i++) {
		tmat(j) = mat3D(i);
		j++;
	}

	return tmat;
}
// ********************* END METODI IDL *************************
}

#endif /* IDLMETHOD_H_ */
