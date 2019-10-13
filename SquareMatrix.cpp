#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <complex>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "ColumnVector.h"
#include "MatrixRow.h"
#include "MatrixColumn.h"

static const double NaN = std::nan("1");

static const double TINY = std::numeric_limits<double>::epsilon();

SquareMatrix::SquareMatrix() : Matrix(), side(0), trace(NaN), inverse(nullptr), coeffs(1), eigenvalues(0) { }

SquareMatrix::SquareMatrix(int s) : Matrix(s, s), side(s), trace(NaN), inverse(nullptr), coeffs(s + 1), eigenvalues(s) { }

SquareMatrix::SquareMatrix(const SquareMatrix & s) : Matrix(s.side, s.side) {
  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = s.matrix[i][j];
    }
  }
  side = s.side;
  trace = s.trace;
  inverse = nullptr;
  if(s.inverse) {
    inverse = new SquareMatrix(*s.inverse);
  }
  coeffs = s.coeffs;
  eigenvalues = s.eigenvalues;
}

SquareMatrix::SquareMatrix(const Matrix & s) : Matrix(s.rows(), s.cols()), coeffs(s.rows() + 1), eigenvalues(s.rows()) {
  if (row != col) throw std::runtime_error("Rows and columns must be equal to convert");
  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = s(i, j);
    }
  }
  side = row;
  trace = NaN;
  inverse = nullptr;
}

SquareMatrix::SquareMatrix(SquareMatrix && s) {
  matrix = s.matrix;
  side = s.side;
  row = s.row;
  col = s.col;
  trace = s.trace;
  inverse = s.inverse;
  coeffs = s.coeffs;
  eigenvalues = s.eigenvalues;

  s.matrix = nullptr;
  s.inverse = nullptr;
  s.side = s.row = s.col = 0;
  s.trace = NaN;
  s.coeffs.clear();
  s.eigenvalues.clear();
}

SquareMatrix::SquareMatrix(Matrix && s) : coeffs(s.rows() + 1), eigenvalues(s.rows()) {
  if (s.rows() != s.cols()) throw std::runtime_error("Rows and columns must be equal to convert");
  matrix = std::move(s.matrix);
  side = s.rows();
  row = std::move(s.rows());
  col = std::move(s.cols());
  trace = NaN;
  inverse = nullptr;

  s.matrix = nullptr;
  s.row = s.col = 0;
}

SquareMatrix& SquareMatrix::operator=(const SquareMatrix & s) {
  if (&s == this) return *this;

  clear();
  side = s.side;
  trace = s.trace;
  row = s.row;
  col = s.col;
  matrix = new double * [row];
  for(int i = 0; i < row; ++i) {
    matrix[i] = new double[col];
  }

  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = s.matrix[i][j];
    }
  }
  inverse = nullptr;
  if(s.inverse) {
    inverse = new SquareMatrix(*s.inverse);
  }
  coeffs = s.coeffs;
  eigenvalues = s.eigenvalues;
  return *this;
}

SquareMatrix& SquareMatrix::operator=(SquareMatrix && s) {
  if (&s == this) return *this;

  clear();
  side = s.side;
  trace = s.trace;
  row = s.row;
  col = s.col;
  matrix = s.matrix;
  inverse = s.inverse;
  coeffs = s.coeffs;
  eigenvalues = s.eigenvalues;
  s.side = s.row = s.col = 0;
  s.trace = NaN;
  s.matrix = nullptr;
  s.inverse = nullptr;
  s.coeffs.clear();
  s.eigenvalues.clear();
  return *this;
}

SquareMatrix& SquareMatrix::operator=(const Matrix & s) {
  if (s.rows() != s.cols()) throw std::runtime_error("Rows and columns must be equal to convert");
  if (&s == this) return *this;

  clear();
  side = s.rows();
  trace = NaN;
  row = side;
  col = side;
  inverse = nullptr;
  coeffs.resize(side + 1);
  eigenvalues.resize(side);
  matrix = new double * [row];
  for(int i = 0; i < row; ++i) {
    matrix[i] = new double[col];
  }

  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = s(i, j);
    }
  }
  return *this;
}

SquareMatrix::~SquareMatrix() {
  if (inverse) delete inverse;
  trace = NaN;
  side = 0;
}


double SquareMatrix::get_determinant() {
  char_poly();
  return coeffs[0];
}

double SquareMatrix::get_trace() {
  char_poly();
  return -coeffs[coeffs.size() - 2];
}

const SquareMatrix& SquareMatrix::get_inverse() {
  if(!get_determinant()) throw std::logic_error("Cannot compute inverse");

  if(inverse != nullptr) delete inverse;
  
  inverse = new SquareMatrix(*this);
  
  std::vector<int> indxc(side);
  std::vector<int> indxr(side);
  std::vector<int> ipiv(side);

  int icol = 0;
  int irow = 0;
  for (auto i = 0; i < side; i++) {
    double big = 0;
    // Look for the pivot element.
    for (auto j = 0; j < side; j++) {
      if (ipiv[j] != 1) {
        for (auto k = 0; k < side; k++) {
          if (ipiv[k] == 0) {
            if (std::abs(inverse->at(j, k)) >= big) {
              big = std::abs(inverse->at(j, k));
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);
    // Found the pivot element. Put the row with the pivot in its diagonal column.
    if (irow != icol) {
      inverse->swapRows(irow, icol);
    }
    indxr[i] = irow;
    indxc[i] = icol;
 
     //Divide the pivot row by the pivot element.
    double pivinv = 1.0 / inverse->at(icol, icol);
    inverse->at(icol, icol) = 1.0;
    for (auto l = 0; l < side; l++) inverse->at(icol, l) *= pivinv;

    // Now, reduce the rows, except for the pivot one.
    for (auto ll = 0; ll < side; ll++) {
      if (ll != icol) {
        double dum = inverse->at(ll, icol);
        inverse->at(ll, icol) = 0.0;
        for (auto l = 0; l < side; l++) inverse->at(ll, l) -= inverse->at(icol, l) * dum;
      }
    }
  }

  // Unscramble.
  for (int l = side - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l]) {
      inverse->swapColumns(indxr[l], indxc[l]);
    }
  }
  return *inverse;
}

SquareMatrix operator/(SquareMatrix& a, SquareMatrix& b) {
  SquareMatrix c = b.get_inverse();
  return a * c;
}

// Balances the matrix. Replaces the original matrix.
void SquareMatrix::balance() {
  const double radix = std::numeric_limits<double>::radix;
  double sqrRadix = radix * radix;
  bool done = false;
  
  while (!done) {
    done = true;
    for (auto i = 0; i < side; ++i) {
      double r = 0;
      double c = 0;
      for (auto j = 0; j < side; ++j) {
        if (j != i) {
          r += std::abs(matrix[i][j]);
          c += std::abs(matrix[j][i]);
        }
      }
      if (c && r) {
        double g = r / radix;
        double f = 1;
        double s = c + r;
        while (c < g) {
          f *= radix;
          c *= sqrRadix;
        }

        g = r * radix;
        while (c > g) {
          f /= radix;
          c /= sqrRadix;
        }
        if ((c + r) / f < 0.95 * s) {
          done = false;
          g = 1 / f;
          //scale[i] *= f;
          for (auto j = 0; j < side; ++j) {
            matrix[i][j] *= g;
            matrix[j][i] *= f;
          }
        }
      }
    }
  }
}

void SquareMatrix::printCoefficients() const {
  std::cout << "Characteristic Polynomial:"<< std::endl;
  std::cout << "     [";
  for(auto i = coeffs.rbegin(); i != coeffs.rend(); i++)
    std::cout << std::setw(11) << *i;
  std::cout << " ]" << std::endl;
}

// Gets the characteristic polynomial of a matrix.
void SquareMatrix::char_poly() {
  SquareMatrix b(side);
  SquareMatrix c(*this);

  for (int i = 0; i < coeffs.size(); i++)
    coeffs[i] = 0;

  auto trace = [](const SquareMatrix& c) -> double {
    double t = 0;
    for (auto i = 0; i < c.side; ++i)
      t += c[i][i];
    return t;
  };

  coeffs[0] = 1;

  for (int l = 0; l < side; l++) {
    if (l) {
      for (int i = 0; i < side; i++) {
        for (int j = 0; j < side; j++) {
          double s = 0;
          for (int k = 0; k < side; k++) {
            s += b[i][k] * matrix[k][j];
          }
          c[i][j] = s;
        }
      }
    }
    double c_trace = trace(c) / (l + 1);
    coeffs[l + 1] = -1 * c_trace * coeffs[0];
    if (l < side - 1) {
      for (int i = 0; i < side; i++) {
        for (int j = 0; j < side; j++) {
          if (j == i)
            b[i][j] = c[i][j] - c_trace;
          else
            b[i][j] = c[i][j];
        }
      }
    }
  }
  for (int i = 0; i < coeffs.size(); i++) {
    if (std::fabs(coeffs[i]) < TINY)
      coeffs[i] = 0;
  }
  std::reverse(coeffs.begin(), coeffs.end());
}

const std::vector<double> SquareMatrix::get_characteristic_polynomial() {
  if (coeffs.back() == 0)
    char_poly();
  //std::vector<double> res = coeffs; 
  return coeffs;
}

void SquareMatrix::factor() {
  char_poly();

  auto f = [](std::vector<double> c, std::complex<double> x) -> std::complex<double> {
    int i = c.size() - 1;
    std::complex<double> p = c[i];
    while (i > 0) 
      p = p * x + c[--i];
    return p;
  };

  auto eq = [](std::vector<std::complex<double>> a, std::vector<std::complex<double>> b) -> bool {
    bool done = true;
    for (int i = 0; i < a.size(); ++i) {
      std::complex<double> delta = a[i] - b[i];
      if (std::fabs(delta) > TINY) {
        done = false;
        break;
      }
    }
    return done;
  };

  auto comp = [](std::complex<double> a, std::complex<double> b) -> bool {
    if (real(a) == real(b))
      return imag(a) > imag(b);
    return real(a) < real(b);
  };

  // Initialize to 0.4 + 0.9i.
  std::complex<double> seed = std::complex<double>(0.4, 0.9);

  for (int i = 0; i < side; i++)
    eigenvalues[i] = std::pow(seed, i);

  // Iterate using the Durand–Kerner method.
  while (true) {
    std::vector<std::complex<double>> temp = eigenvalues;
    for (int i = 0; i < eigenvalues.size(); ++i) {
      std::complex<double> res = 1;
      for (int j = 0; j < eigenvalues.size(); ++j) {
        if (i != j) {
          res *= (eigenvalues[i] - eigenvalues[j]);
        }
      }
      eigenvalues[i] = eigenvalues[i] - f(coeffs, eigenvalues[i]) / res;
    }
    if (eq(eigenvalues, temp)) {
      break;
    }
  }

  const double epsilon = 0.00001;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    if (std::abs(real(eigenvalues[i])) < epsilon) {
      eigenvalues[i] = std::complex<double>(0, imag(eigenvalues[i]));
    }
    if (std::abs(imag(eigenvalues[i])) < epsilon) {
      eigenvalues[i] = std::complex<double>(real(eigenvalues[i]), 0);
    }
  }

  std::sort(eigenvalues.begin(), eigenvalues.end(), comp);
}

void SquareMatrix::printEigenvalues() const {
  std::cout << "Eigenvalues:"<< std::endl;
  std::cout << std::setw(6) << "R [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << real(i);
  std::cout << " ]" << std::endl;

  std::cout << std::setw(6) << "I [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << imag(i);
  std::cout << " ]" << std::endl;
}

// Solve A • x = b for x, given A and b.
// A • x = b <=> x = inverse(A) • b.
ColumnVector SquareMatrix::solve(const ColumnVector& b) {
  if (b.size() != side) throw std::logic_error("Vector size and matrix row size do not match");
  SquareMatrix inv = get_inverse();
  return inv * b;
}

//Solve sets of linear equations in a matrix using column vectors.
Matrix SquareMatrix::solve(const Matrix & b) {
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

Matrix SquareMatrix::get_eigenvectors() {
  using std::fabs;
  factor();

  Matrix res(side, side);
  for (int j = 0; j < eigenvalues.size(); j++) {
    
    double ev;
    if (!imag(eigenvalues[j]))
      ev = real(eigenvalues[j]);
    else {
      double r = real(eigenvalues[j]);
      double i = imag(eigenvalues[j]);
      int sign = (r < 0) ? (i < 0 ? 1 : -1) : (i < 0 ? -1 : 1);

      ev = sign * std::abs(eigenvalues[j]);
    }
    
    SquareMatrix t(*this);
    IdentityMatrix i(side);
    t -= i * ev;
    
    t.reducedRowEchelonForm();

    int index = side;
    while (index > 0) {
      index--;
      MatrixRow mrow(t, index);
      MatrixColumn mcol(t, index);
      
      if (mrow.absSum()) {
        if (std::fabs(eigenvalues[j] - eigenvalues[j + 1]) < TINY)
          j += 2;
        break;
      }

      if (index + 1 < side)
        j++;
      
      MatrixColumn b(res, j);
      b[index] = 1;
      for (int i = 0; i < side; i++) {
        if (mcol[i] && i != index)
          b[i] = -mcol[i];
      }
      double min = fabs(b[0]);
      for (int i = 1; i < side; i++) {
        double val = fabs(b[i]);
        if (val < min && val != 0)
          min = val;
      }
      if (min != 0) {
        double scale = 1 / min;
        for (int k = 0; k < side; k++) {
          b[k] *= scale;
        }
      }
    } 
  }

  return res;
}

// Identity Matrix
IdentityMatrix::IdentityMatrix() : SquareMatrix() { }

IdentityMatrix::IdentityMatrix(int s) : SquareMatrix(s) {
  for (int i = 0; i < side; ++i) {
    matrix[i][i] = 1;
  }
}

IdentityMatrix::~IdentityMatrix() { }

