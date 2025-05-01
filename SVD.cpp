#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "ColumnVector.h"
#include "SVD.h"

using std::pow;
using std::fabs;
using std::sqrt;
using std::min;
using std::fmax;

static const double epsilon = std::numeric_limits<double>::epsilon();

static const int max_iterations = 30;

// Returns a with the sign of b.
static inline double sign(double a, double b) {
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

// Computes [1 + (a^2 + b^2)]^1/2
static double pythag(const double a, const double b) {
  double abs_a = fabs(a);
  double abs_b = fabs(b);
  double result;

  if (abs_a > abs_b) {
    double c = abs_b / abs_a;
    result = abs_a * sqrt(1 + c * c);
  }
  else if (abs_b > 0) {
    double c = abs_a / abs_b;
    result = abs_b * sqrt(1 + c * c);
  }
  else result = 0;
  return result;
}

SVD::SVD(Matrix& m) : u(new Matrix(m)), v(new SquareMatrix(m.cols())), w(m.cols()), row(m.rows()), col(m.cols()) { compute(); }

SVD::~SVD() {
  delete u;
  delete v;
}

void SVD::compute() {
  decompose();
  reorder();
  
  for(int i = 0; i < w.size(); i++)
    if (fabs(w[i]) < epsilon) w[i] = 0;
}

void SVD::printW() const {
  std::cout << "W:"<< std::endl;
  std::cout << "     [";
  for(const auto i: w)
    std::cout << std::setw(11) << i;
  std::cout << " ]" << std::endl;
}

const SquareMatrix SVD::get_w() {
  SquareMatrix wm(w.size());
  for (int i = 0; i < wm.sides(); i++) {
    wm(i, i) = w[i];
  }
  return wm;
}

int SVD::rank() {
  int rank = 0;
  double thresh = 0.5 * sqrt(row + col + 1) * w[0] * epsilon;
  for (int i = 0; i < col; ++i) {
    if (w[i] > thresh) rank++;
  }
  return rank;
}

int SVD::nullity() {
  int nullity = 0;
  double thresh = 0.5 * sqrt(row + col + 1) * w[0] * epsilon;
  for (int i = 0; i < col; ++i) {
    if (w[i] <= thresh) nullity++;
  }
  return nullity;
}

Matrix SVD::range() {
  Matrix range(row, rank());
  double thresh = 0.5 * sqrt(row + col + 1) * w[0] * epsilon;

  int nr = 0;
  for (int j = 0; j < col; j++) {
    if (w[j] > thresh) {
      for (int i = 0; i < row; i++) 
        range[i][nr] = u->at(i, j);
      nr++;
    }
  }
  return range;
}

Matrix SVD::nullspace() {
  Matrix nullspace(col, nullity());
  double thresh = 0.5 * sqrt(row + col + 1) * w[0] * epsilon;

  int nn = 0;
  for (int j = 0; j < col; j++) {
    if (w[j] <= thresh) {
      for (int i = 0; i < col; i++) 
        nullspace[i][nn] = v->at(i, j);
      nn++;
    }
  }
  return nullspace;
}

// Solve Ax = b for x, given U*W*Vt and b.
ColumnVector SVD::solve(const ColumnVector& b) {
  if (b.size() != row) throw std::logic_error("Vector size and matrix row size do not match");

  ColumnVector y(col);
  double thresh = 0.5 * sqrt(row + col + 1) * w[0] * epsilon;
  
  // Calculate Ut * b
  for (int i = 0; i < col; i++) {
    double s = 0;
    if (w[i] > thresh) {
      for (int j = 0; j < row; j++)
        s += u->at(j, i) * b[j];
        s /= w[i];
    }
    y[i] = s;
  }

  // Multiply by V.
  ColumnVector x(col);
  for (int i = 0; i < col; i++) {
    double s = 0;
    for (int j = 0; j < col; j++)
      s += v->at(i, j) * y[i];
    x[i] = s;
  }
  return x;
}

//Solve sets of linear equations in a matrix using column vectors.
Matrix SVD::solve(const Matrix & b) {
  if (b.rows() != row) throw std::logic_error("Matrix sizes do not match");

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

void SVD::decompose() {
  double norm = 0, scale = 0, g = 0;
  double c, f, h, s, x, y, z;
  int k, l, p;
  std::vector<double> rv(col);

  // Householder reduction.
  for (int i = 0; i < col; i++) {
    l = i + 2;
    rv[i] = scale * g;
    g = s = scale = 0;
    if (i < row) {
      for (k = i; k < row; k++) scale += fabs(u->at(k, i));
      if (scale) {
        for (k = i; k < row; k++) {
          u->at(k, i) /= scale;
          s += u->at(k, i) * u->at(k, i);;
        }
      
        f = u->at(i, i);
        g = -sign(sqrt(s), f);
        h = f * g - s;
        u->at(i, i) = f - g;
        for (int j = l - 1; j < col; j++) {
          for (s = 0, k = i; k < row; k++) s += u->at(k, i) * u->at(k, j);
          f = s / h;
          for (k = i; k < row; k++) u->at(k, j) += f * u->at(k, i);
        }
        for (k = i; k < row; k++) u->at(k, i) *= scale;
      }
    }
    w[i] = scale * g;

    // Right-hand reduction.
    g = s = scale = 0;
    if (i + 1 <= row && i + 1 != col) {
      for (k = l - 1; k < col; k++) scale += fabs(u->at(i, k));
      if (scale) {
        for (k = l - 1; k < col; k++){
          u->at(i, k) /= scale;
          s += u->at(i, k) * u->at(i, k);
        }
        f = u->at(i, l - 1);
        g = -sign(sqrt(s), f);
        h = f * g - s;
        u->at(i, l - 1) = f - g;
        for (k = l - 1; k < col; k++) rv[k] = u->at(i, k) / h;
        for (int j = l - 1; j < row; j++) {
          for (s = 0, k = l - 1; k < col; k++) s += u->at(j, k) * u->at(i, k);
          for (k = l - 1; k < col; k++) u->at(j, k) += s * rv[k];
        }
        for (k = l - 1; k < col; k++) u->at(i, k) *= scale;
      }
    }
    norm = fmax(norm, (fabs(w[i]) + fabs(rv[i])));
  }

  // Accumulation of right-hand transformations.
  for (int i = col - 1; i >= 0; i--) {
    if (i < col - 1) {
      if (g) {
        for (int j = l; j < col; j++) {
          v->at(j, i) = (u->at(i, j) / u->at(i, l)) / g;
        }
        for (int j = l; j < col; j++) {
          for (s = 0, k = l; k < col; k++) s += u->at(i, k) * v->at(k, j); 
          for (k = l; k < col; k++) v->at(k, j) += s * v->at(k, i);
        }
      }
      for (int j = l; j < col; j++) v->at(i, j) = v->at(j, i) = 0;
    }
    v->at(i, i) = 1;
    g = rv[i];
    l = i;
  }

  // Accumulation of left-hand transformations.
  for (int i = min(row, col) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    for (int j = l; j < col; j++) u->at(i, j) = 0;
    if (g) {
      g = 1 / g;
      for (int j = l; j < col; j++) {
        for (s = 0, k = l; k < row; k++) s += u->at(k, i) * u->at(k, j);
        f = (s / u->at(i, i)) * g;
        for (k = i; k < row; k++) u->at(k, j) += f * u->at(k, i);
      }
      for (int j = i; j < row; j++) u->at(j, i) *= g;
    }
    else {
      for (int j = i; j < row; j++) u->at(j, i) = 0;
    }
    u->at(i, i) += 1;
  }

  // Diagonalization of the bidiagonal form:
  // Loop over singular values, and over allowed iterations.
  for (k = col - 1; k >= 0; k--) {
    for (int iter = 0; iter < max_iterations; iter++) {
      bool flag = true;
      // Test for splitting.
      for (l = k; l >= 0; l--) {
        p = l - 1;
        if (!l || fabs(rv[l]) <= epsilon * norm) {
          flag = false;
          break;
        }
        if (fabs(w[p]) <= epsilon * norm) break;
      }
      if (flag) {
        c = 0;  // Cancellation of rv[l], if l > 0
        s = 1;
        for (int i = l; i < k + 1; i++) {
          f = s * rv[i];
          rv[i] = c * rv[i];
          if (fabs(f) <= epsilon * norm) break;
          g = w[i];
          h = pythag(f, g);
          w[i] = h;
          h = 1 / h;
          c = g * h;
          s = -f * h;
          for (int j = 0; j < row; j++) {
            y = u->at(j, p);
            z = u->at(j, i);
            u->at(j, p) = y * c + z * s;
            u->at(j, i) = z * c - y * s;
          }
        }
      }
      z = w[k];
      // Convergence: Singular value is made nonnegative.
      if (l == k) {
        if (z < 0) {
          w[k] = -z;
          for (int j = 0; j < col; j++) v->at(j, k) = -v->at(j, k);
        }
        break;
      }
      if (iter == max_iterations) throw std::runtime_error(std::string("Too many iterations for column ") + (char)k + std::string("."));
      // Shift from bottom 2-by-2 minor.
      x = w[l];
      p = k - 1;
      y = w[p];
      g = rv[p];
      h = rv[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
      g = pythag(f, 1);
      f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
      c = 1, s = 1;
      // Next QR transformation.
      for (int j = l; j <= p; j++) {
        int i = j + 1;
        g = rv[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = pythag(f, h);
        rv[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (int jj = 0; jj < col; jj++) {
          x = v->at(jj, j);
          z = v->at(jj, i);
          v->at(jj, j) = x * c + z * s;
          v->at(jj, i) = z * c - x * s;          
        }
        z = pythag(f, h);
        w[j] = z;
        if (z) {
          z = 1 / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (int jj = 0; jj < row; jj++) {
          y = u->at(jj, j);
          z = u->at(jj, i);
          u->at(jj, j) = y * c + z * s;
          u->at(jj, i) = z * c - y * s;
        }
      }
      rv[l] = 0;
      rv[k] = f;
      w[k] = x;
    }
  }
}

void SVD::reorder() {
  std::vector<double> sr(row);
  std::vector<double> sc(col);
  int inc;

  // Sort. Shell's Sort.
  for (inc = 1; inc < col; inc = 3 * inc + 1);
  do {
    inc /= 3;
    for (int i = inc; i < col; i++) {
      double sw = w[i];
      for (int k = 0; k < row; k++) sr[k] = u->at(k, i);
      for (int k = 0; k < col; k++) sc[k] = v->at(k, i);
      int j = i;
      while (w[j - inc] < sw) {
        w[j] = w[j - inc];
        for (int k = 0; k < row; k++) u->at(k, j) = u->at(k, j - inc);
        for (int k = 0; k < col; k++) v->at(k, j) = v->at(k, j - inc);
        j -= inc;
        if (j < inc) break;
      }
      w[j] = sw;
      for (int k = 0; k < row; k++) u->at(k, j) = sr[k];
      for (int k = 0; k < col; k++) v->at(k, j) = sc[k];
    }
  } while (inc > 1);
  // Flip signs.
  for (int k = 0; k < col; k++) {
    int s = 0;
    for (int i = 0; i < row; i++){
      if (u->at(i, k) < 0) s++;
    }
    for (int j = 0; j < col; j++) {
      if (v->at(j, k) < 0) s++;
    }
    if (s > (row + col) / 2) {
      for (int i = 0; i < row; i++) u->at(i, k) = -u->at(i, k);
      for (int j = 0; j < col; j++) v->at(j, k) = -v->at(j, k);
    }
  }
}