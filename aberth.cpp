// includes, system
#include <vector>
#include <cassert>
#include <complex>

// includes, project
#include "aberth.h"

namespace aberth
{

  //----------------------------------------------------------------------------//
  //! constructor with initializer list
  //----------------------------------------------------------------------------//
  template <typename T>
  Polynomial<T>::Polynomial(const std::vector<std::complex<T>> coeffs) :
    coeffs_(coeffs),
    deg_(coeffs.size() - 1) {}

  //----------------------------------------------------------------------------//
  //! return coeffs
  //----------------------------------------------------------------------------//
  template <typename T>
  std::vector<std::complex<T>> Polynomial<T>::getCoeffs() const {return coeffs_;}

  //----------------------------------------------------------------------------//
  //! use Horner's method to evaluate polynomial at x
  //----------------------------------------------------------------------------//
  template <typename T>
  std::complex<T> Polynomial<T>::eval(std::complex<T> x) const
  {

    // initialize return value as leading coeff
    std::complex<T> val = *coeffs_.crbegin();
    for(auto coeffsI = ++coeffs_.crbegin(); coeffsI < coeffs_.crend(); coeffsI++) {
      val = *coeffsI + val*x;
    }

    return val;
  }

  //----------------------------------------------------------------------------//
  //! compute coeffs of this polynomial's derivative
  //----------------------------------------------------------------------------//
  template <typename T>
  std::vector<std::complex<T>> Polynomial<T>::compDerivCoeffs() const
  {

    // return value and copy of this polynomial's coffs
    std::vector<std::complex<T>> derivCoeffs(this->deg_),
                                 coeffs = this->coeffs_;

    // if this polynomial is constant, we return the zero polynomial
    if (!this->deg_) derivCoeffs[0] = static_cast<std::complex<T>>(0);

    // else we compute coeffs
    else {
      for (unsigned int i = 1; i < coeffs.size(); i++) {
        derivCoeffs[i - 1] = coeffs[i] * static_cast<std::complex<T>>(i);
      }
    }

    return derivCoeffs;
  }

  //----------------------------------------------------------------------------//
  //! copy constructor
  //----------------------------------------------------------------------------//
  template <typename T>
  ZeroFinder<T>::ZeroFinder(const ZeroFinder<T> &other) :
    Polynomial<T>::Polynomial(other.getCoeffs()), // call base class constructor
    zeros_(other.zeros_),
    tol_(other.tol_),
    maxIters_(other.maxIters_),
    initMode_(other.initMode_),
    iters_(other.iters_),
    allConv_(other.allConv_),
    conv_(other.conv_),
    invDiffSum_(other.invDiffSum_),
    deriv_(std::make_unique<Polynomial<T>>(other.deriv_->getCoeffs())),
    rev_(std::make_unique<Polynomial<T>>(other.rev_->getCoeffs())),
    revDeriv_(std::make_unique<Polynomial<T>>(other.revDeriv_->getCoeffs()))
  {}

  //----------------------------------------------------------------------------//
  //! move constructor
  //----------------------------------------------------------------------------//
  template <typename T>
  ZeroFinder<T>::ZeroFinder(ZeroFinder<T> &&other) :
    Polynomial<T>::Polynomial(other.getCoeffs()), // call base class constructor
    zeros_(std::move(other.zeros_)),
    tol_(std::move(other.tol_)),
    maxIters_(std::move(other.maxIters_)),
    initMode_(std::move(other.initMode_)),
    iters_(std::move(other.iters_)),
    allConv_(std::move(other.allConv_)),
    conv_(std::move(other.conv_)),
    invDiffSum_(std::move(other.invDiffSum_)),
    deriv_(std::move(other.deriv_)),
    rev_(std::move(other.rev_)),
    revDeriv_(std::move(other.revDeriv_))
  {}

  //----------------------------------------------------------------------------//
  //! swap
  //! @param  first   instance of ZeroFinder<T> to be swapped
  //! @param  second  instance of ZeroFinder<T> to be swapped
  //----------------------------------------------------------------------------//
  template <typename T>
  void ZeroFinder<T>::swap(ZeroFinder<T> &first, ZeroFinder<T> &second)
  {
    std::swap(first.zeros_, second.zeros_);
    std::swap(first.tol_, second.tol_);
    std::swap(first.maxIters_, second.maxIters_);
    std::swap(first.initMode_, second.initMode_);
    std::swap(first.iters_, second.iters_);
    std::swap(first.allConv_, second.allConv_);
    std::swap(first.conv_, second.conv_);
    std::swap(first.invDiffSum_, second.invDiffSum_);
    std::swap(first.deriv_, second.deriv_);
    std::swap(first.rev_, second.rev_);
    std::swap(first.revDeriv_, second.revDeriv_);
  }

  //----------------------------------------------------------------------------//
  //! assignment
  //! @param  other  (possibly distinct) instance of ZeroFinder<T>
  //----------------------------------------------------------------------------//
  template <typename T>
  ZeroFinder<T>& ZeroFinder<T>::operator=(ZeroFinder<T> other)
  {
    if (this != &other) swap(*this, other);
    return *this;
  }

  //----------------------------------------------------------------------------//
  // compute initial guesses for zeros
  // also initialize other values
  //----------------------------------------------------------------------------//
  template <typename T>
  void ZeroFinder<T>::initData()
  {

    // compute initial guesses for zeros
    initZeros();

    // resize and fill vector members
    conv_.resize(this->deg_);
    std::fill(conv_.begin(), conv_.end(), false);

    invDiffSum_.resize(this->deg_);
    std::fill(invDiffSum_.begin(), invDiffSum_.end(), 0);

    // other data
    iters_ = 0;
    allConv_ = false;
  }

  //----------------------------------------------------------------------------//
  // compute initial guesses for zeros
  //----------------------------------------------------------------------------//
  template <typename T>
  void ZeroFinder<T>::initZeros()
  {

    // real and imag component of guesses and upper bound on all zeros' magnitude
    T a,
      b,
      rMax = ZeroFinder<T>::compBound();

    // initial zero for symm
    std::complex<T> z;

    // compute a guess for each root
    zeros_.resize(this->deg_);
    for(unsigned int i = 0; i < this->deg_; i++) {

      // generate components of guess
      if (initMode_ == "rand" or (initMode_ == "symm" and !i)) {
        a = rMax * (static_cast<T>(std::rand()) / static_cast<T>(RAND_MAX));
        b = sqrt(rMax*rMax - a*a) * (static_cast<T>(std::rand()) /
            static_cast<T>(RAND_MAX));
        z = {a, b};
      }

      //**************************************//
      // symmetrically distribute initial guesses
      // by randomly generating a guess, then repeatedly multiplying
      //**************************************//
      if (initMode_ == "symm" and i) z = (z * z) / abs(z);

      //**************************************//
      // randomly distributed initial guesses will repeatedly generate components
      // in the first "if" block of this loop
      //**************************************//
      zeros_[i] = z;
    }
  }

  //----------------------------------------------------------------------------//
  //! obtain cauchy bound on roots' magnitude
  //! https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots
  //----------------------------------------------------------------------------//
  template <typename T>
  T ZeroFinder<T>::compBound() const
  {

    // current ratio, largest ratio, and magnitude of leading coeff
    T curr,
      best = 0,
      lead = abs(this->coeffs_.back());

    // obtain largest ratio
    for (auto coeffsI : this->coeffs_) {
      curr = abs(coeffsI) / lead;
      if (best < curr) best = curr;
    }

    return best + static_cast<T>(1);
  }

  //----------------------------------------------------------------------------//
  //! compute invDiffSum elts
  //----------------------------------------------------------------------------//
  template <typename T>
  void ZeroFinder<T>::compInvDiffSum()
  {

    // compute for each zero
    for(unsigned int i = 0; i < this->deg_; i++) {
      invDiffSum_[i] = {0, 0};
      for(unsigned int j = 0; j < this->deg_; j++) {
        if (i != j) invDiffSum_[i] += static_cast<std::complex<T>>(1) /
                                    (zeros_[i] - zeros_[j]);
      }
    }
  }

  //----------------------------------------------------------------------------//
  //! update zeros
  //! use of reverse polynomial comes from
  //! Bini, Numerical Algorithms 13 (1996) 179-200
  //----------------------------------------------------------------------------//
  template <typename T>
  void ZeroFinder<T>::newtonStep()
  {

    // ratio of polynomial evaluated at z to its derivative evaluated at z
    // correction to be applied
    std::complex<T> ratio,
                    corr;

    // consider each zero
    for(unsigned int i = 0; i < this->deg_; i++) {

      // if magnitude of zero is less than 1, then will likely avoid overflow
      if (abs(zeros_[i]) < 1) ratio = this->eval(zeros_[i]) / deriv_->eval(zeros_[i]);

      // otherwise, we have to use reverse polynomial
      else {
        std::complex<T> gamma = static_cast<std::complex<T>>(1) / zeros_[i];
        ratio = static_cast<std::complex<T>>(1) /
               (static_cast<std::complex<T>>(this->deg_) * gamma - 
               (gamma * gamma * revDeriv_->eval(gamma) / rev_->eval(gamma)));
      }

      // compute and apply new correction
      corr = ratio / (static_cast<std::complex<T>>(1) -
            (ratio * invDiffSum_[i]));
      zeros_[i] -= corr;

      // determine if ith zero has converged after applying correction
      if (!conv_[i]) {
        if (static_cast<T>(abs(this->eval(zeros_[i]))) < tol_) conv_[i] = true;
      }
    }
  }

  //----------------------------------------------------------------------------//
  // specialization
  //----------------------------------------------------------------------------//
  template class Polynomial<float>;
  template class Polynomial<double>;

  template class ZeroFinder<float>;
  template class ZeroFinder<double>;

}
