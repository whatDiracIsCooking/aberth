# This is a C++ implementation of the Aberth method
##### It's a root finding algorithm for polynomials of finite degree

At a minimum, one must provide coefficients as a `std::vector<std::complex<T>>` where `T` is either `float` or `double`.
An initial guess is generated for each root, then updates for each guess are computed until each root is converged
or the maximum iterations are reached.
