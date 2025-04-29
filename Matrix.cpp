#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <limits>
#include <algorithm>
#include "Matrix.h"
#include "ColumnVector.h"
#include "MatrixRow.h"

static const double NaN = std::nan("1");

static const double TINY = std::pow(10, -10);

inline void swap_elements(double& a, double& b) {
  double tmp = a;
  a = b;
  b = tmp;
}

Matrix::Matrix() : row(0), col(0), matrix(nullptr) {}

Matrix::Matrix(int r, int c) : row(r), col(c) {
  matrix = new double * [row];
  for(int i = 0; i < row; ++i) {
    matrix[i] = new double[col];
  }

  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = 0;
    }
  }
}

Matrix::~Matrix() { clear(); }

Matrix::Matrix(const Matrix& m) {
  row = m.row;
  col = m.col;
  matrix = new double * [row];
  for(int i = 0; i < row; ++i) {
    matrix[i] = new double[col];
  }

  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = m.matrix[i][j];
    }
  }
}

Matrix::Matrix(Matrix&& m) : row(m.row), col(m.col), matrix(m.matrix) {
  m.row = 0;
  m.col = 0;
  matrix = nullptr;
}

Matrix& Matrix::operator=(const Matrix& m) {
  if (&m == this) return *this;

  clear();
  row = m.row;
  col = m.col;
  matrix = new double * [row];
  for(int i = 0; i < row; ++i) {
    matrix[i] = new double[col];
  }

  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      matrix[i][j] = m.matrix[i][j];
    }
  }
  return *this;
}

Matrix& Matrix::operator=(Matrix&& m) {
  if (&m == this) return *this;

  clear();
  row = m.row;
  col = m.col;
  matrix = m.matrix;

  m.row = m.col = 0;
  m.matrix = nullptr;
  return *this;
}

Matrix& Matrix::swap(Matrix& m) {
  using std::swap;
  swap(this->row, m.row);
  swap(this->col, m.col);
  swap(this->matrix, m.matrix);
  return *this;
}

void Matrix::swapRows(int a, int b) {
  if (a >= row) throw std::range_error("Row subscript out of range: a");
  if (b >= row) throw std::range_error("Row subscript out of range: b");

  double * temp = matrix[a];
  matrix[a] = matrix[b];
  matrix[b] = temp;
}

void Matrix::swapColumns(int a, int b) {
  if (a >= col) throw std::range_error("Column subscript out of range: a");
  if (b >= col) throw std::range_error("Column subscript out of range: b");

  for (auto i = 0; i < row; ++i) {
    swap_elements(matrix[i][a], matrix[i][b]);
  }
}

Matrix Matrix::transpose() {
  Matrix m(col, row);

  for(int i = 0; i < row; ++i) {
    for(int j = 0; j < col; ++j) {
      m(j, i) = matrix[i][j];
    }
  }
  return m;
}

const double& Matrix::operator()(int r, int c) const {
  check(r, c);
  return matrix[r][c];
}

double& Matrix::operator()(int r, int c) {
  return const_cast<double&>(static_cast<const Matrix&>(*this)(r, c));
}

const double& Matrix::at(int r, int c) const {
  return this->operator()(r, c);
}

double& Matrix::at(int r, int c) {
  return this->operator()(r, c);
}

double * Matrix::operator[](int r) {
  check(r);
  return matrix[r];
}

const double * Matrix::operator[](int r) const {
  check(r);
  return matrix[r];
}

void Matrix::reducedRowEchelonForm() {
  int lead = 0;
  for (int r = 0; r < row; r++) {
    if (col <= lead)
      return;
    int i = r;

    while (matrix[i][lead] == 0) {
      i++;
      if (row == i) {
        i = r;
        lead++;
        if(col == lead) return;
      }
    }
    swapRows(i, r);

    MatrixRow rm(*this, r);
    if (rm[lead] != 0) {
      double div = rm[lead];
      rm /= div;
    }
    for(int k = 0; k < row; k++){
      if (k != r) {
        double mult = matrix[k][lead];
        MatrixRow km(*this, k);
        km -= rm * mult;
      }
    }
    lead++;
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < col; j++) {
        if (std::fabs(matrix[i][j]) < TINY)
          matrix[i][j] = 0;    
      }
    }
  }
}

ColumnVector Matrix::gauss(ColumnVector& b) {
  if (b.size() != row) throw std::range_error("Size of vector does not match rows of matrix");
  
  for (int i = 0; i < row; i++) {
    double max = std::fabs(matrix[i][i]);
    int pivot = i;
    for (int j = i + 1; j < row; j++) {
      if (double tmp = std::fabs(matrix[j][i]) > max) {
        max = tmp;
        pivot = j;
      }
    }
    if (max < TINY)
      //throw std::runtime_error("Matrix is (nearly) singular");
      continue;
    
    swapRows(i, pivot);

    for (int j = i + 1; j < row; j++) {
      double div = matrix[j][i] / matrix[i][i];
      b[j] -= div * b[i];
      for (int c = i; c < row; c++) {
        matrix[j][c] -= div * matrix[i][c];
      }
    }
  }

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (std::fabs(matrix[i][j]) < TINY)
        matrix[i][j] = 0;    
    }
  }
  std::cout << " Gaussian Elimination: " << std::endl << *this;

  /*
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++)
      m(i, j) = matrix[i][j];

  for (int k = 0; k < row; k++)
    m(k, col) = b[k];
  
  m.reducedRowEchelonForm();
  std::cout << "RREF: " << std::endl << m << std::endl;
  
  for (int i = 0; i < row; i++) {
    double max = std::fabs(m(i, i));
    int max_row = i;
    for (int j = i + 1; j < row; j++) {
      if (double tmp = std::fabs(m(j, i)) > max) {
        max = tmp;
        max_row = j;
      }
    }
    if (max == 0)
      continue;
    swapRows(i, max_row);
    for (int j = i + 1; j < row; j++) {
      double div = m(j, i) / m(i, i);
      for (int c = i + 1; c < row; c++) {
        m(j, c) -= div * m(i, c);
      }
      m(j, i) = 0;
      b[j] -= div * b[i];
    }
  }
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (std::fabs(m(i, j)) < TINY)
       m(i, j) = 0;    
    }
  }
  */

  ColumnVector x(row);
  for (int i = row - 1; i >= 0; i--) {
    double sum = 0;
    for (int j = i + 1; j < row; j++) {
      sum += matrix[i][j] * x[j];
    }
    x[i] = (b[i] - sum) / matrix[i][i];
  }
  
  return x;
}

Matrix& Matrix::operator+=(const Matrix & m) {
  if (row != m.row || col != m.col) throw std::logic_error("Rows and columns do not match");
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix[i][j] += m.matrix[i][j];
    }
  }
  return *this;
}

Matrix operator+(const Matrix & a, const Matrix & b) {
  Matrix c(a);
  return c += b;
}

Matrix& Matrix::operator+=(double k) {
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix[i][j] += k;
    }
  }
  return *this;
}

Matrix operator+(const Matrix & a, double b) {
  Matrix c(a);
  return c += b;
}

Matrix operator+(double a, const Matrix & b) {
  Matrix c(b);
  return c += a;
}

Matrix& Matrix::operator-=(const Matrix & m) {
  if (row != m.row || col != m.col) throw std::logic_error("Rows and columns do not match");
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix[i][j] -= m.matrix[i][j];
    }
  }
  return *this;
}

Matrix operator-(const Matrix & a, const Matrix & b) {
  Matrix c(a);
  return c -= b;
}

Matrix operator-(const Matrix & a, double b) {
  Matrix c(a);
  return c -= b;
}

Matrix& Matrix::operator-=(double k) {
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix[i][j] -= k;
    }
  }
  return *this;
}

Matrix& Matrix::operator*=(double k) {
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix[i][j] *= k;
    }
  }
  return *this;
}

Matrix operator*(const Matrix & a, double b) {
  Matrix c(a);
  return c *= b;
}

Matrix operator*(double a, const Matrix & b) {
  Matrix c(b);
  return c *= a;
}


Matrix operator*(Matrix & a, Matrix & b) {
  if (a.cols() != b.rows()) throw std::logic_error("Columns and rows do not match");

  Matrix c(a.rows(), b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      for (int k = 0; k < b.rows(); ++k) {
        c(i, j) += a(i, k) * b(k, j);
      }
      if (std::fabs(c(i, j)) < TINY) c(i, j) = 0;
    }
  }

  return c;
}

Matrix operator*(const Matrix & a, const Matrix & b) {
  if (a.cols() != b.rows()) throw std::logic_error("Columns and rows do not match");
  
  Matrix c(a.rows(), b.cols());
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      for (int k = 0; k < b.rows(); ++k) {
        c(i, j) += a(i, k) * b(k, j);
      }
      if (std::fabs(c(i, j)) < TINY) c(i, j) = 0;
    }
  }
  return c;
}

double dotProduct(const Matrix & a, const Matrix & b) {
  if (a.rows() != b.rows() or a.cols() != b.cols())
     throw std::logic_error("Rows and columns do not match");

  double dotProduct = 0;
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      dotProduct += a(i, j) * b(i, j);
    }
  }
  return dotProduct;
}

Matrix& Matrix::operator/=(double k) {
  if (!k) throw std::runtime_error("Attempted to divide by zero");
  for (int i = 0; i < row; ++i) {
    for (int j = 0; j < col; ++j) {
      matrix[i][j] /= k;
    }
  }
  return *this;
}

Matrix operator/(const Matrix & a, double b) {
  Matrix c(a);
  return c /= b;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
  for(int i = 0; i < m.row; ++i) {
    os << std::setw(4) << i << " [";
    for(int j = 0; j < m.col; ++j) {
      os << std::setw(11) << m.matrix[i][j];
    }
    os << std::setw(2) << "]" << std::endl;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Matrix& m) {
  std::cout << "Enter " << m.row * m.col << " values: " << std::endl;
  for(int i = 0; i < m.row; ++i) {
    for(int j = 0; j < m.col; ++j) {
      is >> m.matrix[i][j];
    }
  }
  return is;
}

void Matrix::check(int r, int c) const {
  if(r < 0 || r >= row) throw std::range_error("Row subscript out of range");
  if(c < 0 || c >= col) throw std::range_error("Column subscript out of range");
  return;
}

void Matrix::check(int r) const {
  if(r < 0 || r >= row) throw std::range_error("Row subscript out of range");
  return;
}

void Matrix::clear() {
  for(int i = 0; i < row; ++i) {
    delete [] matrix[i];
  }
  delete [] matrix;
  row = 0;
  col = 0;
}
