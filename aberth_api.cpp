// includes, system
#include <vector>
#include <string>
#include <cassert>
#include <complex>
#include <memory>
#include <cfloat>
#include <algorithm>

#include "aberth.h"

namespace aberth
{

  //----------------------------------------------------------------------------//
  //! non-default constructor
  //----------------------------------------------------------------------------//
  template <typename T>
  ZeroFinder<T>::ZeroFinder(const std::vector<std::complex<T>> coeffs,
      const T tol,
      const unsigned int maxIters,
      const std::string initMode) :
    Polynomial<T>::Polynomial(coeffs), // call base class constructor
    tol_(tol),
    maxIters_(maxIters),
    initMode_(initMode)
  {

    // tolerance must be above smallest positive double
    assert(DBL_MIN <= tol and "tol must be at least DBL_MIN");

    // only two methods of initializing zeros 
    assert((initMode == "rand" or initMode == "symm") and
        "initMode must be either 'rand' or 'symm'");

    // initialize data
    initData();

    // initialize pointer to the polynomial's derivative
    deriv_ = std::make_unique<Polynomial<T>>(this->compDerivCoeffs());

    // instantiate reverse polynomial and its derivative
    std::vector<std::complex<T>> revCoeffs(coeffs.size());
    for(unsigned int i = 0; i < coeffs.size(); i++) {
      revCoeffs[i] = coeffs[coeffs.size() - 1 - i];
    }
    rev_ = std::make_unique<Polynomial<T>>(revCoeffs);
    revDeriv_ = std::make_unique<Polynomial<T>>(rev_->compDerivCoeffs());
  }
  
  //----------------------------------------------------------------------------//
  //! use Aberth's method to approximate zeros
  //! @param  verbose  whether or not extra info will be printed
  //! returns true iff all zeros have converged to within tol_
  //----------------------------------------------------------------------------//
  template <typename T>
  bool ZeroFinder<T>::compZeros(const bool verbose)
  {

    // run until converged or maxIters is reached
    for(; iters_ < maxIters_ and !allConv_; iters_++) {

      // compute diffSum then update zeros
      compInvDiffSum();
      newtonStep();

      // determine if all zeros have been converged
      allConv_ = true;
      for(const bool convI : conv_) {if (!convI) {allConv_ = false; break;};}
    }

    // sort zeros by their magnitude
    sort(zeros_.begin(), zeros_.end(), [](const auto v, const auto w)
        {return abs(v) < abs(w);});

    // print error message if not all zeros have converged
    if (iters_ == maxIters_ and !allConv_)
      printf("Failed to converge all zeros after maximum (%u) iterations\n",
             maxIters_);

    // if requested, print verbose output
    if (verbose) {

      // obtain maximum error and index at which it occured
      T err = 0, maxErr = 0;
      unsigned int errInd = 0, maxInd = 0;
      for (auto zerosI : zeros_) {
        err = static_cast<T>(abs(this->eval(zerosI)));
        if (maxErr < err) {maxErr = err; maxInd = errInd;}
        errInd++;
      }
      printf("Iterations performed = %u\n", iters_);
      printf("Index of max error   = %u\n", maxInd);
      printf("Max error            = %.15f\n", maxErr);
    }

    return allConv_;
  }

  //----------------------------------------------------------------------------//
  //! return zeros of polynomial
  //----------------------------------------------------------------------------//
  template <typename T>
  std::vector<std::complex<T>> ZeroFinder<T>::getZeros() const {return zeros_;}

  //----------------------------------------------------------------------------//
  // specialization
  //----------------------------------------------------------------------------//
  template class Polynomial<float>;
  template class Polynomial<double>;

  template class ZeroFinder<float>;
  template class ZeroFinder<double>;

}
