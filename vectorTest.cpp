#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <complex>
#include "Matrix.h"
#include "ColumnVector.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 21 - 10;
    }
  }
}

void load_random_values(ColumnVector& col) {
  for (int i = 0; i < col.size(); ++i) {
    col[i] = rand() % 21 - 10;
  }
}

int main() {
  srand(time(NULL));
  std::cout << "Testing Matrix and Vectors..." << std::endl;

  Matrix m(3, 3);
  std::cin >> m;

  std::cout << "Matrix: " << std::endl << m;

  ColumnVector c(3);
  std::cin >> c;

  std::cout << "Vector: " << std::endl << c;

  ColumnVector res = m * c;
  std::cout << res;
  return 0;
}
