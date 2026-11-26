// This diagnostic adapter includes an external GPLv3 CLASS-independent
// thermal spectral-function example, gw3g/spectral v1.1-spf. Not vendored.
// Build against that separately obtained source tree.
#define main spectral_example_default_main
#include "examples/dileptons.cpp"
#undef main
#include <iomanip>

int main(int argc, char** argv) {
  if (argc == 2) tol = atof(argv[1]);
  if (argc > 2 || !(tol>0.0)) return 2;
  Rho_V vector_trace;
  double M, p;
  while (std::cin >> M >> p) {
    if (!(M >= 4. && p > 0.)) return 2;
    k = p;
    k0 = hypot(M,p);
    vector_trace();
    std::cout << std::setprecision(16) << M << " " << p << " "
              << vector_trace.lo << " " << vector_trace.nlo << std::endl;
  }
  return std::cin.eof() ? 0 : 2;
}
