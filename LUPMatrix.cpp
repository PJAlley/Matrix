#include <iostream>
#include <cmath>
#include "SquareMatrix.h"
#include "LUPMatrix.h"
#include "RowVector.h"
#include "ColumnVector.h"

static const double TINY = std::pow(10, -10);

// LU Matrix
LUPMatrix::LUPMatrix(SquareMatrix& s) : m(new SquareMatrix(s)), 
  lower(new IdentityMatrix(s.sides())), 
  upper(new SquareMatrix(s.sides())), 
  p(new IdentityMatrix(s.sides())), 
  side(s.sides()), sign(1) { compute(); }

LUPMatrix::~LUPMatrix() {
  delete upper;
  delete lower;
  delete p;
}

/*
  Computes the LUP Decomposition.
  The result should be such that PA = LU.
*/
void LUPMatrix::compute() {
  for (int k = 0; k < side; k++) {
    int pivot_row = k;
    double pivot_value = std::fabs(m->at(pivot_row, pivot_row));
    for (int i = k + 1; i < side; i++) {
      double value = std::fabs(m->at(i, k));
      if (value > pivot_value) {
        pivot_row = i;
        pivot_value = value;
      }
    }

    if (k == pivot_row) {
      if (pivot_value == 0) {
        std::cout << "Singular." << std::endl;
        return;
      }
    }
    else {  //Swap rows.
      m->swapRows(k, pivot_row);
      p->swapRows(k, pivot_row);
      sign = -sign;
    }
    for (int i = k + 1; i < side; i++) { //iterate down rows
      double mult = m->at(i, k) / m->at(k, k);
      m->at(i, k) = mult;
      for (int j = k + 1; j < side; j++) {
        // subtract off lower triangle factor times pivot row
        m->at(i, j) = m->at(i, j) - (mult * m->at(k, j));
      }
    }
    for(int j = k; j < side; j++) {
      upper->at(k, j) = m->at(k, j);
    }
    for(int j = 0; j < k; j++) {
      lower->at(k, j) = m->at(k, j);
    }
  }
}

double LUPMatrix::determinant() {
  double determinant = sign;
  for (auto i = 0; i < side; ++i) {
    determinant *= upper->at(i, i);
  }
  return determinant;
}

/* Solve Ax = b for x, given L, U, P and b.
  First, interchange b with the permutation matrix.
  Then, solve Ly = b for y.
  Finally, solve Ux = y for x.

  Equation: 
  A • x = (L • U) • x = L • (U • x) = b
  L • y = b
  U • x = y
*/
ColumnVector LUPMatrix::solve(const ColumnVector& b) {
  if (b.size() != side) throw std::logic_error("Vector size and matrix size do not match");

  // Change values based on permutation matrix.
  ColumnVector c = *p * b;

  // Solve Ly = b for y.
  ColumnVector y(side);
  //y[0] = c[0];
  for (int i = 0; i < side; i++) {
    double sum = 0;
    for (int j = 0; j < i; j++) {
      sum += lower->at(i, j) * y[j];
    }
    y[i] = (c[i] - sum);
  }

  // Now solve Ux = y for x.
  ColumnVector x(side);
  //x[side - 1] = y[side - 1] / upper->at(side - 1, side - 1);
  for (int i = side - 1; i >= 0; i--) {
    double sum = 0;
    for (int j = i + 1; j < side; j++) {
      sum += upper->at(i, j) * x[j];
    }
    x[i] = (y[i] - sum) / upper->at(i, i);
  }
  return x;
}

//Solve sets of linear equations in a matrix using column vectors.
Matrix LUPMatrix::solve(const Matrix & b) {
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

