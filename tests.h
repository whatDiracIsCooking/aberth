/* Functions to generate coeffs for test polynomials*/
#ifndef _TESTS_H_
#define _TESTS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <vector>
#include <time.h>
#include <complex>

/*
 * Use some test cases from 
 * https://numpi.dm.unipi.it/mpsolve-2.2/mpsolve.pdf
 */

namespace aberth
{

  //----------------------------------------------------------------------------//
  //! easy
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> easyCoeffs(const unsigned int deg);

  //----------------------------------------------------------------------------//
  //! truncated exponential
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> expCoeffs(const unsigned int deg);

  //----------------------------------------------------------------------------//
  //! kam1
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> kam1coeffs(const double c);

  //----------------------------------------------------------------------------//
  //! roots of unity
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> unityRootsCoeffs(const unsigned int deg);

}


#endif  // #ifndef _TESTS_H_
