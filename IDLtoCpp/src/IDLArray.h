/*
 * IDLArray.h
 *
 *  Created on: Jul 12, 2013
 *      Author: remocean
 */

#ifndef IDLARRAY_H_
#define IDLARRAY_H_

#include <boost/multi_array.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>

using namespace std;

namespace idl2 {

typedef unsigned size_type;

enum _DIM {
	_1D, _2D, _3D, _UNDEF
};

template<class T>
class Array: private boost::multi_array<T, 3> { //eredito pvt cosi nn espongo i metodi del multiarray
public:

	/**
	 *
	 */
	Array<T>() :
			boost::multi_array<T, 3>() {
		dimensions = _UNDEF;
		_size1 = -1;
		_size2 = -1;
		_size3 = -1;
	}

	/**
	 *
	 * @param x
	 */
	Array<T>(size_type x) :
			boost::multi_array<T, 3>(boost::extents[x][0][0]) {
		dimensions = _1D;
		_size1 = x;
		_size2 = -1;
		_size3 = -1;
		setToZero();
	}

	/**
	 *
	 * @param x
	 * @param y
	 */
	Array<T>(size_type x, size_type y) :
			boost::multi_array<T, 3>(boost::extents[x][y][0]) {
		dimensions = _2D;
		_size1 = x;
		_size2 = y;
		_size3 = -1;
		setToZero();
	}

	/**
	 *
	 * @param x
	 * @param y
	 * @param z
	 */
	Array<T>(size_type x, size_type y, size_type z) :
			boost::multi_array<T, 3>(boost::extents[x][y][z]) {
		dimensions = _3D;
		_size1 = x;
		_size2 = y;
		_size3 = z;
		setToZero();
	}

	/**
	 *
	 * @param x
	 */
	inline
	void resize(size_type x) {
		boost::multi_array<T, 3>::resize(boost::extents[x][0][0]);

		dimensions = _1D;
		_size1 = x;
		_size2 = -1;
		_size3 = -1;
	}

	/**
	 *
	 * @param x
	 * @param y
	 */
	inline
	void resize(size_type x, size_type y) {
		boost::multi_array<T, 3>::resize(boost::extents[x][y][0]);

		dimensions = _2D;
		_size1 = x;
		_size2 = y;
		_size3 = -1;
	}

	/**
	 *
	 * @param x
	 * @param y
	 * @param z
	 */
	inline
	void resize(size_type x, size_type y, size_type z) {
		boost::multi_array<T, 3>::resize(boost::extents[x][y][z]);

		dimensions = _3D;
		_size1 = x;
		_size2 = y;
		_size3 = z;
	}

	/**
	 * Per l'accesso singolo come vettore 1D
	 */
	inline T & operator ()(size_type i) {
		return boost::multi_array<T, 3>::operator [](i, 0, 0);
	}

	/**
	 * Per l'accesso singolo come vettore 2D
	 */
	inline T & operator ()(size_type i, size_type j) {
		return boost::multi_array<T, 3>::operator [](i, j, 0);
	}

	/**
	 * Per l'accesso singolo come vettore 3D
	 */
	inline T & operator ()(size_type i, size_type j, size_type k) {
		return boost::multi_array<T, 3>::operator [](i, j, k);
	}

	/**
	 *
	 * @return
	 */
	_DIM getDimensionality() const {
		return dimensions;
	}

	/**
	 *
	 * @return
	 */
	size_type size1() const {
		return _size1;
	}

	/**
	 *
	 * @return
	 */
	size_type size2() const {
		return _size2;
	}

	/**
	 *
	 * @return
	 */
	size_type size3() const {
		return _size3;
	}

	/**
	 *
	 */
	void setToZero() {
		for (size_type i = 0; i != _size1; ++i) {
			for (size_type j = 0; j != _size2; ++j) {
				for (size_type k = 0; k != _size3; ++k) {
					this->operator ()(i, j, k) = 0;
				}
			}
		}
	}

	/**
	 *
	 * @param mat
	 * @return
	 */
	T max(const Array<T> &mat) {
		T maximum = 0;

		for (size_type i = 0; i != _size1; ++i) {
			for (size_type j = 0; j != _size2; ++j) {
				for (size_type k = 0; k != _size3; ++k) {
					if (mat(i, j, k) > maximum)
						maximum = mat(i, j, k);
				}
			}
		}
		return maximum;
	}

	/**
	 *
	 */
	T min(const Array<T> &mat) {
		T min = numeric_limits<T>::max();

		for (size_type i = 0; i != _size1; ++i) {
			for (size_type j = 0; j != _size2; ++j) {
				for (size_type k = 0; k != _size3; ++k) {
					if (mat(i, j, k) < min)
						min = mat(i, j, k);
				}
			}
		}
		return min;
	}

private:
	size_type _size1, _size2, _size3;
	_DIM dimensions;
};

namespace type {
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

} //end namespace type

namespace methods {

/**
 *
 * @param mat
 * @return
 */
template<class T>
static T sum_elem(const Array<T> &mat) {
	T sum = 0;
	for (size_type i = 0; i != mat.size1(); ++i) {
		for (size_type j = 0; j != mat.size2(); ++j) {
			for (size_type k = 0; k != mat.size3(); ++k) {
				sum += mat(i, j, k);
			}
		}
	}
	return sum;
}

/**
 *
 * @param dim
 * @return
 */
static Array<FLOAT> findgen(int dim) {
	Array<FLOAT> a(dim);
	FLOAT value = 0;

	for (size_type i = 0; i < dim; ++i) {
		a(i) = value++;
	}

	return a;
}

/**
 *
 * @return Floatarray a una dimensione.
 */
static Array<FLOAT> fltarr(size_type dim) {
	return Array<FLOAT>(dim);
}

/**
 *
 * @return Floatarray a due dimensioni.
 */
static Array<FLOAT> fltarr(size_type m, size_type n) {
	return Array<FLOAT>(m, n);
}

/**
 *
 * @return Floatarray a 3 dimensioni.
 */
static Array<FLOAT> fltarr(size_type m, size_type n, size_type k) {
	return Array<FLOAT>(m, n, k);
}

/**
 *
 * @return Longarray a una dimensione.
 */
static Array<LONG> lonarr(size_type dim) {
	return Array<LONG>(dim);
}

/**
 *
 * @return Longarray a due dimensioni.
 */
static Array<LONG> lonarr(size_type x, size_type y) {
	return Array<LONG>(x, y);
}

/**
 *
 * @return Longarray a 3 dimensioni.
 */
static Array<LONG> lonarr(size_type x, size_type y, size_type z) {
	return Array<LONG>(x, y, z);
}

/**
 *
 */
template<class T>
static Array<T> abs(const Array<T> &mat) {

	for (size_type i = 0; i != mat.size1(); ++i) {
		for (size_type j = 0; j != mat.size2(); ++j) {
			for (size_type k = 0; k != mat.size3(); ++k) {
				mat(i, j, k) = std::abs(mat(i, j, k));
			}
		}
	}
	return mat;
}

/**
 *
 */
template<class T>
static Array<T> transpose(const Array<T> &mat) {
	size_type s1 = mat.size1();
	size_type s2 = mat.size2();
	size_type s3 = mat.size3();

	Array<T> tmat(s1, s2, s3);
	for (size_type i = 0; i != s1; ++i) {
		for (size_type j = 0; j != s2; ++j) {
			for (size_type k = 0; k != s3; ++k) {
				tmat(j, i, k) = mat(i, j, k); //fisso l'asse z
			}
		}
	}
	return tmat;
}

/**
 *
 */
template<class T, class E>
Array<T> minus(const Array<T> &mat, const E &t) {
	size_type s1 = mat.size1();
	size_type s2 = mat.size2();
	size_type s3 = mat.size3();

	Array<T> tmat(s1, s2, s3);
	for (size_type i = 0; i != s1; ++i) {
		for (size_type j = 0; j != s2; ++j) {
			for (size_type k = 0; k != s3; ++k) {
				tmat(i, j, k) = mat(i, j, k) - t;
			}
		}
	}
	return tmat;
}

/**
 *
 */
template<class T, class E>
Array<T> mult(const Array<T> &mat, const E &t) {
	size_type s1 = mat.size1();
	size_type s2 = mat.size2();
	size_type s3 = mat.size3();

	Array<T> tmat(s1, s2, s3);
	for (size_type i = 0; i != s1; ++i) {
		for (size_type j = 0; j != s2; ++j) {
			for (size_type k = 0; k != s3; ++k) {
				tmat(i, j, k) = mat(i, j, k) * t;
			}
		}
	}
	return tmat;
}

} //end namespace methods

} //end namespace idl

#endif /* IDLARRAY_H_ */
