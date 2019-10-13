#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "LUPMatrix.h"
#include "LUMatrix.h"
#include "ColumnVector.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 10;
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
  std::cout << "Testing LUP..." << std::endl;
  //int size = 3;
  //SquareMatrix sm(size);
  //load_random_values(sm);
  //std::cout << "sm:" << std::endl << sm;

  int size;
  std::cout << "Enter side: ";
  std::cin >> size;
  SquareMatrix sm(size);
  std::cin >> sm;
  std::cout << "sm:" << std::endl << sm;

  ColumnVector b(size);
  load_random_values(b);
  std::cout << "b: " << b.size() << std::endl << b;

  std::cout << "Determinant (Square Matrix):  " << sm.get_determinant() << std::endl;
  LUPMatrix lup(sm);
  std::cout << "Determinant (LUP Matrix):  " << lup.determinant() << std::endl;
  std::cout << "Upper: " << std::endl << *(lup.get_upper());
  std::cout << "Lower: " << std::endl << *(lup.get_lower());
  std::cout << "Permutation Matrix: " << std::endl << *(lup.get_permutation_matrix());

  const SquareMatrix lur = *(lup.get_lower()) * *(lup.get_upper());
  std::cout << "L * U : " << std::endl << lur;

  SquareMatrix pa = *(lup.get_permutation_matrix()) * sm;
  std::cout << "P * A : " << std::endl << pa;

  std::cout << "Testing LU..." << std::endl;
  LUMatrix lu(sm);
  std::cout << "Determinant (LU Matrix):  " << lu.determinant() << std::endl;
  std::cout << "Upper: " << std::endl << *(lu.get_upper());
  std::cout << "Lower: " << std::endl << *(lu.get_lower());
  const SquareMatrix lurp = *(lu.get_lower()) * *(lu.get_upper());
  std::cout << "L * U : " << std::endl << lurp;

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  ColumnVector x = lup.solve(b);
  std::cout << "x (LUP):" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  std::cout << "sm:" << std::endl << sm;
  Matrix m(size, 3);
  load_random_values(m);
  std::cout << "m:" << std::endl << m;
  Matrix res = lup.solve(m);
  std::cout << "res (LUP):" << std::endl << res;

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  ColumnVector c = lu.solve(b);
  std::cout << "x (LU):" << std::endl << c;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  Matrix resm = lu.solve(m);
  std::cout << "res (LU):" << std::endl << resm;

  std::cout << std::endl << "Inverse: " << std::endl;
  IdentityMatrix i(size);
  Matrix lupi = lup.solve(i);
  std::cout << "Inverse (LUP):" << std::endl << lupi;
  Matrix lui = lu.solve(i);
  std::cout << std::endl << "Inverse (LU):" << std::endl << lui;

  std::cout << std::endl << "Inverse: " << std::endl << sm.get_inverse();

  return 0;
}

//2 7 6 9 5 1 4 3 8
//5 4 2 1 0 1 -1 -1 -1 -1 3 0 1 1 -1 2
//2 4 -6 0 4 6 -3 -4 0 0 4 0 0 4 -6 2