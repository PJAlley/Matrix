#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 11 + 1;
    }
  }
}

int main() {
  srand(time(NULL));
  std::cout << "Testing Matrix..." << std::endl;

  SquareMatrix sm(5);
  load_random_values(sm);
  std::cout << "sm:" << std::endl << sm;

  SquareMatrix r = sm.get_inverse();
  std::cout << "r:" << std::endl << r;

  SquareMatrix t = sm * r;
  std::cout << "t:" << std::endl << t;
  //std::cout << "inv r:" << std::endl << r.get_inverse();
  return 0;
}