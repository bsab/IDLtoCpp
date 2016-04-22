/*
 * Main.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: remocean
 */

#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "IDLNamespace.h"

using namespace std;
namespace fs = boost::filesystem;
using namespace boost::numeric::ublas;

int main() {
	cout << "START" << endl;

	try {
		 idl::Arr<idl::FLOAT> a(4);
		 a(0) = 10;
		 a(1) = 5;
		 a(2) = -3;
		 a(3) = -8;
		 cout << "a:" << a << endl;
		 a.setToZero();
		 cout << "a.clear:" << a << endl;

		 idl::Arr<idl::FLOAT> a2D(3,2);
		 //riga 0
		 a2D(0,0) = 0;
		 a2D(0,1) = 1;
		 //riga 1
		 a2D(1,0) = 2;
		 a2D(1,1) = 3;
		 //riga 2
		 a2D(2,0) = 4;
		 a2D(2,1) = 5;
		 cout << "a2D:" << a2D << endl;
		 idl::Arr<idl::FLOAT> b2D = idl::rebin(a2D, 6, 4);

		 idl::Arr<idl::FLOAT> z = idl::brange(a2D, 0, 2, 1, 2);
		 cout << "z:" << z << endl;

	} catch (int exp) {
		std::cerr << "IDL Exception : " << ExceptionMsg(exp) << std::endl;
	}

	cout << endl << "END" << endl;

	return 0;
}
