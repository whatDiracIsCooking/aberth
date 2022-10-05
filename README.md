# The Aberth method in C++

### Summary
This implements a root finding algorithm for univariate polynomials of finite degree.
One provides the coefficients of their polynomial `p(x)` as a `std::vector<std::complex<T>>` where `T` is either `float` or `double`.
Roots are computed and stored in a `std::vector<std::complex<T>>`.

The method works by generating an initial guess for each root, then [each guess is updated](https://en.wikipedia.org/wiki/Aberth_method)
until they are all converged (or the maximum iterations are reached).

Convergence is obtained when `|p(x)| < tolerance` for each root `x`.
Roots are sorted by non-decreasing magnitude.

### Example 1
The zeros of `p(x) = x^2 - x - 6` are approximated as follows:

```
// declare vector for zeros and initialize coeffs
std::vector<std::complex<double>> zeros, coeffs = {{-6, 0}, {-1, 0}, {1, 0}};

// instantiate zero finder
aberth::ZeroFinder<double> *finder = new aberth::ZeroFinder<double>(coeffs);

// compute then copy zeros
finder->compZeros();
zeros = finder->getZeros();
```

### Inital guess generation
There are two methods for inital guess generation, `rand` and `symm`.
In either case, [Cauchy's bound](https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots) `R`
on the magnitude of all zeros is computed.
Then, for the default method `rand`, a guess with magnitude bounded by `R`
is produced for each zero. Alternatively, for `symm`, a guess `z` is produced with `|z| <= R`.
Subsequent guesses are produced by setting `z <- z * z / |z|`, creating a regular n-gon centered at the origin.


### Example 2
The complex roots of `p(x) = x^4 - 1` are computed below at single precision, with tolerance set to `1e-12`,
maximum iterations set to `100`, and initial guesses generated with `symm`. Verbose messages are printed as well.
```
std::vector<std::complex<float>> zeros, coeffs = {{-1, 0}, {0, 0}, {0,0}, {0,0}, {1, 0}};
aberth::ZeroFinder<float> *finder = new aberth::ZeroFinder<float>(coeffs, 1e-12, 100, "symm");
finder->compZeros(true);
zeros = finder->getZeros();
```

Default values for the tolerance, maximum iterations, and method of generating initial guess are detailed in `aberth.h`.

### Recommended use
The implementation here isn't at all sophisticated. It's best used for **friendly** polynomials:
* Roots are well-separated or exactly degenerate 
* Tolerance is at least `1e-12`
* Neither overflow nor underflow occurs when evaluating the polynomial at points whose magnitude is at most Cauchy's bound
  * Degree is less than `100`
  * Coefficients are well-within numeric limits defined in the `<limits>` header
