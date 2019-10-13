#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "SVD.h"
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
  std::cout << "Testing SVD..." << std::endl;
  Matrix m(4, 5);
  load_random_values(m);
  std::cout << "m:" << std::endl << m;

  SVD svd(m);
  svd.compute();

  Matrix r = svd.range();
  std::cout << "Ramge:" << std::endl << r;

  Matrix n = svd.nullspace();
  std::cout << "Nullspace:" << std::endl << n;

  std::cout << "U(sv):" << std::endl << *(svd.get_u());
  std::cout << "W(sv):" << std::endl << svd.get_w();
  std::cout << "Vt(sv):" << std::endl << svd.get_vt();

  ColumnVector b(4);
  load_random_values(b);
  std::cout << "b: " << b.size() << std::endl << b;

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  ColumnVector x = svd.solve(b);
  std::cout << "x:" << std::endl << x;

  Matrix s(4, 3);
  load_random_values(s);
  std::cout << "s:" << std::endl << s;
  std::cout << std::endl << "Solving for matrix..." << std::endl;
  Matrix y = svd.solve(s);
  std::cout << "y:" << std::endl << y;

  return 0;
}