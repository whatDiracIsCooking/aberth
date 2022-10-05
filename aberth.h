/* Structs to organize variables */
#ifndef _ABERTH_H_
#define _ABERTH_H_

// includes, system
#include <vector>
#include <string>
#include <complex>
#include <memory>

namespace aberth
{

  //----------------------------------------------------------------------------//
  //! Simple polynomial class
  //----------------------------------------------------------------------------//
  template <typename T>
  class Polynomial
  {

    //**************************************
    // data members
    //**************************************
    
    protected:

      //! coefficients
      std::vector<std::complex<T>> coeffs_;

      //! degree of polynomial
      unsigned int deg_;

    //**************************************
    // member functions
    //**************************************

    public:

      //! constructor
      Polynomial(const std::vector<std::complex<T>> coeffs);

      //! return coeffs
      std::vector<std::complex<T>> getCoeffs() const ;

      //! use Horner's method to evaluate polynomial at x
      std::complex<T> eval(const std::complex<T> x) const;

      //! compute coeffs of this polynomial's derivative
      std::vector<std::complex<T>> compDerivCoeffs() const;

  };

  //----------------------------------------------------------------------------//
  //! Derived from Polynomial class
  //----------------------------------------------------------------------------//
  template <typename T>
  class ZeroFinder : public Polynomial<T>
  {

    //**************************************//
    // data members
    //**************************************//

    protected:

      //! (not necessarily converged) zeros of polynomial
      std::vector<std::complex<T>> zeros_;

      //! tolerance for calculation
      T tol_;

      //! maximum number of iterations
      unsigned int maxIters_;

      //! method for initializing zeros
      std::string initMode_;

      //! number of iterations performed
      unsigned int iters_;

      //! whether or not all zeros are converged to within tol
      bool allConv_;

      //! whether or not individual zeros are converged to within tol
      std::vector<bool> conv_;

      //! summed reciprocal differences of current zeros
      std::vector<std::complex<T>> invDiffSum_;

      //! instance of this polynomial's derivative
      std::unique_ptr<Polynomial<T>> deriv_;

      //! instance of this polynomial's reverse
      std::unique_ptr<Polynomial<T>> rev_;

      //! instance of the derivative of this polynomial's reverse
      std::unique_ptr<Polynomial<T>> revDeriv_;


    //**************************************//
    // member functions
    //**************************************//

    public:

      //! non-default constructor
      ZeroFinder(const std::vector<std::complex<T>> coeffs,
          const T tol = 1e-9,
          const unsigned int maxIters = 200,
          const std::string initMode = "rand");

      //! default destructor suffices
      ~ZeroFinder() = default;

      //! copy constructor
      ZeroFinder(const ZeroFinder<T> &other);

      //! move constructor
      ZeroFinder(ZeroFinder<T> &&other);

      //! use copy and swap
      void swap(ZeroFinder<T> &first, ZeroFinder<T> &second);

      //! assignment
      ZeroFinder<T>& operator=(ZeroFinder<T> other);

      //! use Aberth's method to approximate zeros
      bool compZeros(const bool verbose = false);

      //! return (not necessarily converged) zeros of polynomial
      std::vector<std::complex<T>> getZeros() const;

    protected:

      //! initialize data
      void initData();

      //! compute initial guesses for zeros
      void initZeros();

      //! get upper bound on magnitude of zeros
      T compBound() const;

      //! compute invDiffSum elts
      void compInvDiffSum();

      //! update zeros with a Newton step
      void newtonStep();

  };
}

#endif  // #ifndef _ABERTH_H_
