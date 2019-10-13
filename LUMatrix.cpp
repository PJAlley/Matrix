#include <iostream>
#include <cmath>
//#include "Matrix.h"
#include "SquareMatrix.h"
#include "LUMatrix.h"
#include "RowVector.h"
#include "ColumnVector.h"

static const double TINY = std::pow(10, -10);

// LU Matrix
LUMatrix::LUMatrix(SquareMatrix& s) : m(&s), lower(new IdentityMatrix(s.sides())), 
  upper(new SquareMatrix(s.sides())), side(s.sides()) { compute(); }

LUMatrix::~LUMatrix() {
  delete upper;
  delete lower;
}

// Computes the LU Decomposition: LU = A.
void LUMatrix::compute() {
  for (auto i = 0; i < side; ++i) {
    for (auto j = i; j < side; ++j) {
      double sum = 0;
      for (auto k = 0; k < i; k++) {
        sum += lower->at(i, k) * upper->at(k, j);
      }
      upper->at(i, j) = m->at(i, j) - sum;
      if (i == j && upper->at(i, j) == 0) {
        upper->at(i, j) = TINY;
      }
    }

    for (auto j = i + 1; j < side; ++j) {
      double sum = 0;
      for (auto k = 0; k < i; ++k) {
        sum += lower->at(j, k) * upper->at(k, i);
      }
      lower->at(j, i) = (m->at(j, i) - sum) / upper->at(i, i);
    }
  }
}

double LUMatrix::determinant() {
  double determinant = side % 2 ? -1 : 1;
  for (auto i = 0; i < side; ++i) {
    determinant *= upper->at(i, i);
  }
  return determinant;
}

/* Solve Ax = b for x, given LU and b.
  First, solve Ly = b for y.
  Then solve Ux = y for x.

  Equation: 
  A • x = (L • U) • x = L • (U • x) = b
  L • y = b
  U • x = y
*/
ColumnVector LUMatrix::solve(const ColumnVector& b) {
  if (b.size() != side) throw std::logic_error("Vector size and matrix size do not match");

  // Solve Ly = b for y.
  ColumnVector y(side);
  y[0] = b[0] / lower->at(0, 0);
  for (int i = 1; i < side; i++) {
    double sum = 0;
    for (int j = 0; j < i; j++) {
      sum += lower->at(i, j) * y[j];
    }
    y[i] = (b[i] - sum) / lower->at(i, i);
  }

  // Now solve Ux = y for x.
  ColumnVector x(side);
  x[side - 1] = y[side - 1] / upper->at(side - 1, side - 1);
  
  for (int i = side - 2; i >= 0; i--) {
    double sum = 0;
    for (int j = i + 1; j < side; j++) {
      sum += upper->at(i, j) * x[j];
    }
    x[i] = (y[i] - sum) / upper->at(i, i);
  }
  return x;
}


//Solve sets of linear equations in a matrix using column vectors.
Matrix LUMatrix::solve(const Matrix & b) {
  if (b.rows() != side) throw std::logic_error("Matrix sizes do not match");

  Matrix res(b.rows(), b.cols());
  for (int i = 0; i < b.cols(); i++) {
    ColumnVector cvec(b, i);
    ColumnVector cres = solve(cvec);

    for (int j = 0; j < b.rows(); j++) {
      res(j, i) = cres[j];
    }
  }
  return res;
}
