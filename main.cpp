#include <vector>
#include <cassert>
#include <complex>

// includes, project
#include "aberth.h"
#include "tests.h"

int main()
{

  namespace ab = aberth;

  // vars for storing coefficients and zeros of polynomials
  std::vector<std::complex<double>> coeffs,
                                    zeros;

  // seed rand
  std::srand(8008335);

  // make sure tests pass
  bool passed = true;

  printf("Testing \"Easy\" polynomials...\n");
  for(unsigned int i = 1; i <= 10; i ++) {

    // generate coeffs
    coeffs = ab::easyCoeffs(10 * i);

    // initialize class
    ab::ZeroFinder<double> *finder = new ab::ZeroFinder<double>(coeffs);
    
    // compute zeros
    passed = finder->compZeros();

    // copy zeros (for demonstration)
    std::vector<std::complex<double>> zeros = finder->getZeros();

    // die if test failed
    if (!passed) {printf("Failed \"Easy\" polynomial with degree %u\n", i); exit(1);}
  }

  // the rest of the tests are like the one above

  printf("Testing \"Exp\" polynomials...\n");
  for(unsigned int i = 1; i <= 10; i++) {
    coeffs = ab::expCoeffs(i);
    ab::ZeroFinder<double> *finder = new ab::ZeroFinder<double>(coeffs);
    passed = finder->compZeros();
    zeros = finder->getZeros();
    if (!passed) {printf("Failed \"Exp\" polynomial with degree %u\n", i); exit(1);}
  }

  printf("Testing \"kam1\" polynomials...\n");
  for(unsigned int i = 1; i <= 10; i++) {
    double c = 1e-6 / (i*10);
    coeffs = ab::kam1coeffs(c);
    ab::ZeroFinder<double> *finder = new ab::ZeroFinder<double>(coeffs);
    passed = finder->compZeros();
    zeros = finder->getZeros();
    if (!passed) {printf("Failed \"kam1\" polynomial with c = %f\n", c); exit(1);}
  }

  printf("Testing \"roots of unity\" polynomials...\n");
  for(unsigned int i = 1; i <= 10; i++) {
    coeffs = ab::unityRootsCoeffs(i);
    ab::ZeroFinder<double> *finder = new ab::ZeroFinder<double>(coeffs);
    passed = finder->compZeros();
    zeros = finder->getZeros();
    if (!passed) {printf("Failed \"roots of unity\" polynomial with degree %u\n", i); exit(1);}
  }

  if (passed) printf("All tests passed! :)\n");
//
}
