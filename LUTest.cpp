#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "LUMatrix.h"
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
  std::cout << "Testing LU..." << std::endl;
  //int size = 4;
  //SquareMatrix sm(size);
  //load_random_values(sm);
  //std::cout << "sm:" << std::endl << sm;

  int size;
  std::cout << "Enter side: ";
  std::cin >> size;
  SquareMatrix sm(size);
  std::cin >> sm;
  std::cout << "sm:" << std::endl << sm;

  std::vector<double> poly = sm.get_characteristic_polynomial();

  std::cout << "Characteristic Polynomial: ";
  std::cout << "     [";
  for(auto i = poly.rbegin(); i != poly.rend(); i++)
    std::cout << std::setw(11) << *i;
  std::cout << " ]" << std::endl;

  ColumnVector b(size);
  load_random_values(b);
  std::cout << "b: " << b.size() << std::endl << b;

  LUMatrix lu(sm);
  std::cout << "Determinant (LU Decomposition): " << lu.determinant() << std::endl;
  std::cout << "Determinant (Square Matrix):    " << sm.get_determinant() << std::endl;
  std::cout << "Upper: " << std::endl << *(lu.get_upper());
  std::cout << "Lower: " << std::endl << *(lu.get_lower());
  const SquareMatrix lur = *(lu.get_lower()) * *(lu.get_upper());
  std::cout << "L * U : " << std::endl << lur;

  std::cout << std::endl << "Solving for column vector..." << std::endl;
  ColumnVector x = lu.solve(b);
  std::cout << "x:" << std::endl << x;

  std::cout << std::endl << "Solving for matrix..." << std::endl;
  std::cout << "sm:" << std::endl << sm;
  Matrix m(size, 3);
  load_random_values(m);
  std::cout << "m:" << std::endl << m;
  Matrix res = lu.solve(m);
  std::cout << "res:" << std::endl << res;

  return 0;
}