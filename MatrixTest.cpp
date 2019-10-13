#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <complex>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "RowVector.h"
#include "ColumnVector.h"
#include "Eigen.h"
#include "Hessen.h"
#include "SVD.h"
#include "QRMatrix.h"
#include "LUMatrix.h"
#include "LUPMatrix.h"
#include "Schur.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 21 - 10;
    }
  }
}

void load_random_values(RowVector& row) {
  for (int i = 0; i < row.size(); ++i) {
    row[i] = rand() % 21 - 10;
  }
}

void load_random_values(ColumnVector& col) {
  for (int i = 0; i < col.size(); ++i) {
    col[i] = rand() % 21 - 10;
  }
}

int main() {
  srand(time(NULL));
  std::cout << "Testing Matrix..." << std::endl;
  int size = 6;
  SquareMatrix sm(size);
  load_random_values(sm);
  std::cout << "sm:" << std::endl << sm;

  std::cout << "inverse(sm):" << std::endl << sm.get_inverse();
  SquareMatrix sb(size);
  load_random_values(sb);
  std::cout << "sb:" << std::endl << sb;
  std::cout << "d:" << std::endl << sm * sb;

  double dp = dotProduct(sm, sb);
  std::cout << "dp:" << dp << std::endl;

  double det = sb.get_determinant();
  std::cout << "Determinant: " << det << std::endl;
  std::cout << "sb:" << std::endl << sb;

  sb.factor();
  sb.printCoefficients();
  sb.printEigenvalues();

  // LU Decomposition
  std::cout << std::endl << "LU Decomposition:" << std::endl;
  std::cout << "sm:" << std::endl << sm;
  LUMatrix lu(sm);
  std::cout << "Determinant: " << lu.determinant() << std::endl;
  std::cout << "Upper: " << std::endl << *(lu.get_upper());
  std::cout << "Lower: " << std::endl << *(lu.get_lower());
  const SquareMatrix lur = *(lu.get_lower()) * *(lu.get_upper());
  std::cout << "L * U : " << std::endl << lur;

  // LUP Decomposition
  std::cout << std::endl << "LUP Decomposition:" << std::endl;
  std::cout << "sm:" << std::endl << sm;
  LUPMatrix lup(sm);
  std::cout << "Determinant: " << lup.determinant() << std::endl;
  std::cout << "Upper: " << std::endl << *(lup.get_upper());
  std::cout << "Lower: " << std::endl << *(lup.get_lower());
  std::cout << "Permutation Matrix: " << std::endl << *(lup.get_permutation_matrix());
  const SquareMatrix lupr = *(lup.get_lower()) * *(lup.get_upper());
  std::cout << "L * U : " << std::endl << lupr;
  std::cout << "P * A : " << std::endl << *(lup.get_permutation_matrix()) * sm;

  // QR Decomposition
  std::cout << std::endl << "QR Decomposition:" << std::endl;
  std::cout << "sm:" << std::endl << sm;
  QRMatrix qr(sm);
  std::cout << "QT: " << std::endl << *(qr.get_qt());
  std::cout << "R : " << std::endl << *(qr.get_r());
  const SquareMatrix qtr = *(qr.get_qt()) * *(qr.get_r());
  std::cout << "Qt * R : " << std::endl << qtr;
  const SquareMatrix qtq = *(qr.get_qt()) * *(qr.get_q());
  std::cout << "Qt * Q : " << std::endl << qtq;


  Eigen hess(sm);
  std::cout << std::endl << "Eigen Decomposition:" << std::endl;
  std::cout << "sm:" << std::endl << sm;
  try {
    std::cout << "Evector(sb):" << std::endl << *(hess.get_eigenvector());
    hess.printEigenvalues();
  }
  catch(std::exception& e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  // Schur Decomposition
  Schur schur(sm);
  std::cout << std::endl << "Schur Decomposition:" << std::endl;
  std::cout << "sm:" << std::endl << sm;
  try {
    std::cout << "Evector(sb):" << std::endl << *(schur.get_eigenvector());
    schur.printEigenvalues();
    SquareMatrix d = schur.get_diagonal_matrix();
    std::cout << "Eigenvalue Matrix(sb):" << std::endl << d;
    SquareMatrix q = *(schur.get_eigenvector()) * d * schur.get_eigenvector()->get_inverse();
    std::cout << "V * D * V' : " << std::endl << q;
  }
  catch(std::exception& e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  // SVD
  std::cout << std::endl << "Single Value Decomposition:" << std::endl;
  std::cout << "sm:" << std::endl << sm;
  SVD svd(sm);
  try {
    std::cout << "U(sb):" << std::endl << *(svd.get_u());
    std::cout << "W(sb):" << std::endl << svd.get_w();
    std::cout << "V(sb):" << std::endl << *(svd.get_v());
  }
  catch(std::exception& e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  Matrix sv(5, 6);
  load_random_values(sv);
  std::cout << "sv:" << std::endl << sv;
  SVD svdm(sv);
  try {
    std::cout << "U(sv):" << std::endl << *(svdm.get_u());
    std::cout << "W(sv):" << std::endl << svdm.get_w();
    std::cout << "V(sv):" << std::endl << *(svdm.get_v());
    //svdm.printW();
    Matrix ds = (*(svdm.get_u()) * svdm.get_w()) * svdm.get_vt();
    std::cout << "ds:" << std::endl << ds;
  }
  catch(std::exception& e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  std::cout << std::endl << "Row Vector:" << std::endl;
  RowVector rv(sv, 2);
  std::cout << "rv:" << std::endl << rv;

  RowVector rv2(sv[3], sv.cols());
  std::cout << "rv2:" << std::endl << rv2;

  RowVector rv3 = rv2;
  rv3[1] = 15;
  std::cout << "rv3:" << std::endl << rv3;
  std::cout << "rv2:" << std::endl << rv2;

  std::cout << std::endl << "Column Vector:" << std::endl;
  ColumnVector cv(sv, 2);
  std::cout << "cv:" << std::endl << cv;

  ColumnVector colvec(size);
  load_random_values(colvec);
  std::cout << "sm:" << std::endl << sm;
  std::cout << "colvec:" << std::endl << colvec;

  ColumnVector res = sm * colvec;
  std::cout << "res:" << std::endl << res;

  double resr = dotProduct(rv3, res);
  std::cout << "Dot Product: " << resr << std::endl;

  Matrix cross = outerProduct(colvec, res);
  std::cout << "cross:" << std::endl << cross;

  Matrix scross = outerProduct(colvec, rv3);
  std::cout << "scross:" << std::endl << scross;

  return 0;
}
