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

#include "tests.h"



namespace aberth
{

  using namespace std::complex_literals;

  //----------------------------------------------------------------------------//
  //! easy
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> easyCoeffs(const unsigned int deg)
  {

    // assert constraint on degree of polynomial
    assert(deg <= 100 and "Degree of polynomial for easyCoeffs must be <= 100");

    std::vector<std::complex<double>> coeffs(deg + 1);
    for(unsigned int i = 0; i <= deg; i++) {
      coeffs[i] = static_cast<std::complex<double>>(i + 1);
    }
    return coeffs;
  }

  //----------------------------------------------------------------------------//
  //! truncated exponential
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> expCoeffs(const unsigned int deg)
  {

    // assert constraint on degree of polynomial
    assert(deg <= 10 and "Degree of polynomial for expCoeffs must be <= 10");

    std::vector<std::complex<double>> coeffs(deg + 1);
    std::complex<double> k = {1, 0};
    for(unsigned int i = 0; i <= deg; i++) {
      if (2 <= i) k *= static_cast<std::complex<double>>(i);
      coeffs[i] = static_cast<std::complex<double>>(1) / k;
    }
    return coeffs;
  }

  //----------------------------------------------------------------------------//
  //! kam1
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> kam1coeffs(const double c)
  {

    // assert constraints on c
    assert(1e-20 < c and c < 1e-6 and "c must be in [1e-20, 1e-6]");

    unsigned int deg = 7;
    std::vector<std::complex<double>> coeffs(deg + 1);
    for(unsigned int i = 0; i <= deg; i++) {
      switch (i) {
        case 0: 
          coeffs[i] = static_cast<std::complex<double>>(9 * pow(c, 4)); break;
        case 1:
          coeffs[i] = static_cast<std::complex<double>>(6 * pow(c, 2)); break;
        case 2:
          coeffs[i] = static_cast<std::complex<double>>(1); break;
        case 7:
          coeffs[i] = static_cast<std::complex<double>>(c) * 1i; break;
        default:
          coeffs[i] = static_cast<std::complex<double>>(0); break;
      }
    }
    return coeffs;
  }


  //----------------------------------------------------------------------------//
  //! roots of unity
  //----------------------------------------------------------------------------//
  std::vector<std::complex<double>> unityRootsCoeffs(const unsigned int deg)
  {

    std::vector<std::complex<double>> coeffs(deg + 1);
    for(unsigned int i = 0; i <= deg; i++) {
      if (!i) coeffs[i] = static_cast<std::complex<double>>(-1);
      else if (i == deg) coeffs[i] = static_cast<std::complex<double>>(1);
      else coeffs[i] = static_cast<std::complex<double>>(0);
    }
    return coeffs;

  }

}
