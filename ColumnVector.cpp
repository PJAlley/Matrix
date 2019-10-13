#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include "Matrix.h"
#include "RowVector.h"
#include "ColumnVector.h"

ColumnVector::ColumnVector() : length(0) {}

ColumnVector::ColumnVector(int c) : col(c), length(c) {}

ColumnVector::ColumnVector(const ColumnVector& c) {
  length = c.size();
  col = c.col;
}

// Get the r-th column of the matrix and load it in the vector.
// Length is the number of rows in the matrix.
ColumnVector::ColumnVector(const Matrix& m, int c) : length(m.rows()) {
  for (int i = 0; i < length; i++)
    col.push_back(m[i][c]);
}

void ColumnVector::set_size(int c) {
  if (c <= 0) throw std::range_error("Length cannot be negative.");
  col.resize(c);
  length = c;
}

const double& ColumnVector::operator[](int c) const {
  if (c >= col.size())  throw std::range_error("Subscript out of range");
  return col[c];
}

double& ColumnVector::operator[](int c) {
  return const_cast<double&>(static_cast<const ColumnVector&>(*this)[c]);
}

std::ostream& operator<<(std::ostream& os, const ColumnVector& c) {
  std::cout << "     [";
  for(const auto i: c.col)
    std::cout << std::setw(11) << i;
  std::cout << " ]" << std::endl;
  return os;
}

std::istream& operator>>(std::istream& is, ColumnVector& c) {
  std::cout << "Enter " << c.length << " values: " << std::endl;
  for (int i = 0; i < c.length; i++) {
    is >> c.col[i];
  }
  return is;
}

// Multiply a matrix by a column vector. Returns a column vector.
ColumnVector operator*(const Matrix& a, const ColumnVector& b) {
  if (a.cols() != b.size()) throw std::logic_error("Row size and column size do not match");

  ColumnVector c(a.rows());

  for (auto i = 0; i < a.rows(); i++) {
    for (auto j = 0; j < b.size(); ++j) {
      c[i] += a(i, j) * b[j];
    }
  }
  return c;
}

double dotProduct(const RowVector& r, const ColumnVector& c) {
  if (r.size() != c.size()) throw std::logic_error("Row size and column size do not match");
  double res = 0;

  for (auto i = 0; i < c.size(); ++i)
    res += r[i] * c[i];
  
  return res;
}

Matrix outerProduct(const ColumnVector& a, const ColumnVector& b) {
  Matrix res(a.size(), b.size());

  for (auto i = 0; i < a.size(); ++i) 
    for (auto j = 0; j < b.size(); ++j)
      res(i, j) = a[i] * b[j];
  
  return res;
}

Matrix outerProduct(const ColumnVector& a, const RowVector& b) {
  Matrix res(a.size(), b.size());

  for (auto i = 0; i < a.size(); ++i) 
    for (auto j = 0; j < b.size(); ++j)
      res(i, j) = a[i] * b[j];
  
  return res;
}
