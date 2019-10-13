#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "ColumnVector.h"
#include "SVD.h"
#include "QRMatrix.h"
#include "LUMatrix.h"

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
  std::cout << "Testing decompositions..." << std::endl;

  int size = 4;
  SquareMatrix sm(size);
  load_random_values(sm);
  std::cout << "sm:" << std::endl << sm;

  ColumnVector b(size);
  load_random_values(b);
  std::cout << "b: " << b.size() << std::endl << b;

  Matrix m(size, 3);
  load_random_values(m);
  std::cout << "m:" << std::endl << m;

  std::cout << std::endl << "Solving matrix..." << std::endl;
  std::cout << std::endl << "Solving for column vector..." << std::endl;
  ColumnVector x = sm.solve(b);
  std::cout << "x:" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  Matrix res = sm.solve(m);
  std::cout << "res:" << std::endl << res;

  std::cout << std::endl << "Testing LU..." << std::endl;
  LUMatrix lu(sm);
  std::cout << "Upper: " << std::endl << *(lu.get_upper());
  std::cout << "Lower: " << std::endl << *(lu.get_lower());

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  x = lu.solve(b);
  std::cout << "x:" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  res = lu.solve(m);
  std::cout << "res:" << std::endl << res;

  std::cout << std::endl << "Testing SVD..." << std::endl;
  SVD svd(sm);
  std::cout << "U(sv):" << std::endl << *(svd.get_u());
  std::cout << "W(sv):" << std::endl << svd.get_w();
  std::cout << "Vt(sv):" << std::endl << svd.get_vt();

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  x = svd.solve(b);
  std::cout << "x:" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  res = svd.solve(m);
  std::cout << "res:" << std::endl << res;

  std::cout << std::endl << "Testing QR..." << std::endl;
  QRMatrix qr(sm);
  qr.compute();
  std::cout << "QT: " << std::endl << *(qr.get_qt());
  std::cout << "R : " << std::endl << *(qr.get_r());

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  x = qr.solve(b);
  std::cout << "x:" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  res = qr.solve(m);
  std::cout << "res:" << std::endl << res;
  
  return 0;
}