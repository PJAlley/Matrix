#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "Eigen.h"
#include "Hessen.h"
#include "QRMatrix.h"
#include "Schur.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 21 ;
    }
  }
}

int main() {
  srand(time(NULL));
  //int size = 4;
  //SquareMatrix sm(size);
  //load_random_values(sm);

  int side;
  std::cout << "Enter side: ";
  std::cin >> side;
  SquareMatrix sm(side);
  std::cin >> sm;

  std::cout << "Matrix:" << std::endl << sm;

  Hessen sh(sm);
  std::cout << "Hessen(sm):" << std::endl;
  sh.printEigenvalues();
  std::cout << "Eigenvector(sm):" << std::endl << *(sh.get_eigenvector());
  SquareMatrix diag = sh.get_diagonal_matrix();
  std::cout << "Diagonal(sm):" << std::endl << diag;
  SquareMatrix inv = sh.get_eigenvector()->transpose();
  std::cout << "Inverse(sm):" << std::endl << inv;
  SquareMatrix answ = *(sh.get_eigenvector()) * diag * inv;
  std::cout << "Matrix(sb):" << std::endl << answ;
  std::cout << std::endl;


  Eigen se(sm);
  std::cout << "Eigen(sm):" << std::endl;
  se.printEigenvalues();
  std::cout << "Eigenvector(sm):" << std::endl << *(se.get_eigenvector());
  SquareMatrix diag2 = se.get_diagonal_matrix();
  std::cout << "Diagonal(sm):" << std::endl << diag2;
  SquareMatrix inv2 = se.get_eigenvector()->get_inverse();
  std::cout << "Inverse(sb):" << std::endl << inv2;
  SquareMatrix ans = *(se.get_eigenvector()) * diag2 * inv2;
  std::cout << "Matrix(sm):" << std::endl << ans;

  std::cout << std::endl << "Matrix:" << std::endl << sm;

  Schur ss(sm);
  std::cout << std::endl << "Schur:" << std::endl;

  ss.printEigenvalues();
  SquareMatrix h = ss.get_diagonal_matrix();
  std::cout << "Eigenvector(sm):" << std::endl << *(ss.get_eigenvector());
  SquareMatrix v = ss.get_eigenvector()->get_inverse();
  SquareMatrix q = *(ss.get_eigenvector()) * h * v;
  std::cout << "Eigenvector(sm):" << std::endl << *(ss.get_eigenvector());
  std::cout << "Diagonal(sm):" << std::endl << h;
  std::cout << "Inverse(sb):" << std::endl << v;
  std::cout << "Q(sm):" << std::endl << q;
  //3 1 -4 2 -2 4 3 -1 1 2 -2 3 1 -1 4 2

  //5 4 2 1 0 1 -1 -1 -1 -1 3 0 1 1 -1 2

  //2 4 -6 0 4 6 -3 -4 0 0 4 0 0 4 -6 2

  //3 -1 5 2 -3 -9 5 2 1

  //-1 4 -2 -3 4 0 -3 1 3

  return 0;
}