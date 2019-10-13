#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <vector>
#include "SquareMatrix.h"
#include "Schur.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 11 ;
    }
  }
}

int main() {
  //srand(time(NULL));
  //int size = 4;
  //SquareMatrix sm(size);
  //load_random_values(sm);

  int side;
  std::cout << "Enter side: ";
  std::cin >> side;
  SquareMatrix sm(side);
  std::cin >> sm;

  std::vector<double> poly = sm.get_characteristic_polynomial();

  std::cout << "Characteristic Polynomial: ";
  std::cout << "     [";
  for(auto i = poly.rbegin(); i != poly.rend(); i++)
    std::cout << std::setw(11) << *i;
  std::cout << " ]" << std::endl;

  std::cout << "Determinant: " << sm.get_determinant() << std::endl;
  std::cout << "Trace: " << sm.get_trace() << std::endl;

  std::cout << "Using factor():" << std::endl;
  std::cout << "Matrix:" << std::endl << sm;
  sm.factor();
  sm.printCoefficients();
  sm.printEigenvalues();
  std::cout << "Eigenvectors: " << std::endl << sm.get_eigenvectors() << std::endl;

  std::cout << "Using Schur:" << std::endl;
  std::cout << "Matrix:" << std::endl << sm;
  Schur ss(sm);
  ss.printEigenvalues();
  std::cout << "Eigenvectors: " << std::endl << *(ss.get_eigenvector());

  return 0;
}
  // 3 1 -4 2 -2 4 3 -1 1 2 -2 3 1 -1 4 2

  // 5 4 2 1 0 1 -1 -1 -1 -1 3 0 1 1 -1 2

  // 2 4 -6 0 4 6 -3 -4 0 0 4 0 0 4 -6 2

  // 5 2 -6 -1 0 1 3 1 -4

  // 1 -3 3 3 -5 3 6 -6 4

  // -2 -4 2 -2 1 2 4 2 5

  // 4 8 -2 -3 -6 1 9 12 -5

  // 0 -6 -4 5 -11 -6 -6 9 4

  // 3 -1 0 -1 2 -1 0 -1 3

  // 5 8 16 4 1 8 -4 -4 -11
