#include <iostream>
#include <cmath>
#include <vector>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "ColumnVector.h"
#include "QRMatrix.h"

// Returns a with the sign of b.
static inline double sign(double a, double b) {
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

QRMatrix::QRMatrix(SquareMatrix& s) : q(new IdentityMatrix(s.sides())), r(new SquareMatrix(s)), 
  qt(nullptr), side(s.sides())  { compute(); }

QRMatrix::~QRMatrix() {
  delete q;
  delete qt;
  delete r;
}

void QRMatrix::compute() {
  int side = r->sides();
  std::vector<double> c(side);
  std::vector<double> d(side);

  for (auto k = 0; k < side - 1; ++k) {
    double scale = 0;
    for (auto i = k; i < side; ++i) {
      scale = std::fmax(scale, std::abs(r->at(i, k)));
    }
    if (scale == 0) {
      c[k] = d[k] = 0;
    }
    else {
      for (auto i = k; i < side; ++i) {
        r->at(i, k) /= scale;
      }
      double isum = 0;
      for (auto i = k; i < side; ++i) {
        isum += std::pow(r->at(i, k), 2);
      }
      double sigma = sign(std::sqrt(isum), r->at(k, k));
      r->at(k, k) += sigma;
      c[k] = sigma * r->at(k, k);
      d[k] = -1 * scale * sigma;

      for (auto j = k + 1; j < side; ++j) {
        double jsum = 0;
        for (auto i = k; i < side; ++i) {
          jsum += r->at(i, k) * r->at(i, j);
        }
        double tau = jsum / c[k];
        for (auto i = k; i < side; ++i) {
          r->at(i, j) -= tau * r->at(i, k);
        }
      }
    }
  }
  d[side - 1] = r->at(side - 1, side - 1);
  
  for (auto k = 0; k < side - 1; ++k) {
    if (c[k] != 0) {
      for (auto j = 0; j < side; ++j) {
        double sum = 0;
        for (auto i = k; i < side; ++i)
          sum += r->at(i, k) * q->at(i, j);
        sum /= c[k];
        for (auto i = k; i < side; ++i)
          q->at(i, j) -= sum * r->at(i, k);
      }
    }
  }
  qt = new SquareMatrix(q->transpose());
  for (auto i = 0; i < side; ++i) {
    r->at(i, i) = d[i];
    for (auto j = 0; j < i; j++) r->at(i, j) = 0;
  }
}

/* Solve Ax = b for x, given QR and b.
  QR should be computed first using compute().
  First, form Qt • b.
  Then solve R • x = Qt • b for x.

  Equation: 
  Qt • b = y
  R • x = y
*/
ColumnVector QRMatrix::solve(const ColumnVector& b) {
  if (b.size() != side) throw std::logic_error("Vector size and matrix row size do not match");

  // Calculate Qt • b and place in y.
  ColumnVector y(side);
  for (int i = 0; i < side; i++) {
    double sum = 0;
    for (int j = 0; j < side; j++)
      sum += qt->at(i, j) * b[j];
    y[i] = sum;
  }

  // Solve for x: R • x = y = Qt • b.
  ColumnVector x(side);
  for (int i = side - 1; i >= 0; i--) {
    double sum = b[i];
    for (int j = i + 1; j < side; j++)
      sum -= r->at(i, j) * y[j];
    x[i] = sum / r->at(i, i);
  }
  return x;
}

//Solve sets of linear equations in a matrix using column vectors.
Matrix QRMatrix::solve(const Matrix & b) {
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
