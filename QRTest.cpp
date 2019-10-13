#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "QRMatrix.h"
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
  std::cout << "Testing QR..." << std::endl;
  int size = 4;
  SquareMatrix sm(size);
  load_random_values(sm);
  std::cout << "sm:" << std::endl << sm;

  ColumnVector b(size);
  load_random_values(b);
  std::cout << "b: " << b.size() << std::endl << b;

  QRMatrix qr(sm);
  std::cout << "QT: " << std::endl << *(qr.get_qt());
  std::cout << "R : " << std::endl << *(qr.get_r());

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  ColumnVector x = qr.solve(b);
  std::cout << "x:" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  std::cout << "sm:" << std::endl << sm;
  Matrix m(size, 3);
  load_random_values(m);
  std::cout << "m:" << std::endl << m;
  Matrix res = qr.solve(m);
  std::cout << "res:" << std::endl << res;

  return 0;
}