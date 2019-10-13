#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <limits>
#include <algorithm>
#include "MatrixColumn.h"
#include "Matrix.h"

MatrixColumn::MatrixColumn(int l) : column(new double [l]), matrix(nullptr), length(l), index(-1), alloc(true) {}

MatrixColumn::MatrixColumn(Matrix& m, int r) : column(new double [m.cols()]), matrix(&m), length(m.cols()), index(r), alloc(true) {
  for (int i = 0; i < length; i++)
    column[i] = matrix->at(i, r);
}

MatrixColumn::MatrixColumn(double * d, int l) : column(d), matrix(nullptr), length(l), index(-1), alloc(false) {}

MatrixColumn::MatrixColumn(const MatrixColumn& m) {
  alloc = true;
  length = m.length;
  matrix = nullptr;
  index = -1;
  column = new double[m.length];
  for (int i = 0; i < length; i++)
    column[i] = m[i];
}

MatrixColumn& MatrixColumn::operator=(const MatrixColumn & m) {
  if (&m == this) return *this;
  if (alloc) delete [] column;

  alloc = true;
  length = m.length;
  index = m.index;
  column = new double[m.length];
  for (int i = 0; i < length; i++)
    column[i] = m[i];
  matrix = m.matrix;
  return *this;
}

MatrixColumn::~MatrixColumn() {
  if (alloc) {
    setMatrix();
    delete [] column;
  }
}

const double& MatrixColumn::operator[](int r) const {
  if(r < 0 || r >= length) throw std::range_error("Subscript out of range");
  return column[r];
}

double& MatrixColumn::operator[](int r) {
  double res = const_cast<double&>(static_cast<const MatrixColumn&>(*this)[r]);
  setMatrix();
  return column[r];
}

void MatrixColumn::setMatrix() {
  if (matrix != nullptr) {
    for (int i = 0; i < length; i++)
      matrix->at(i, index) = column[i];
  }
}

double MatrixColumn::sum() const {
  double sum = 0;
  for (int i = 0; i < length; i++) 
    sum += column[i];
  return sum;
}

double MatrixColumn::absSum() const {
  double sum = 0;
  for (int i = 0; i < length; i++) 
    sum += std::fabs(column[i]);
  return sum;
}

double MatrixColumn::product() const {
  double product = 1;
  for (int i = 0; i < length; i++) 
    product *= column[i];
  return product;
}

MatrixColumn& MatrixColumn::operator+=(double d) {
  for (int i = 0; i < length; i++)
    column[i] += d;
  
  setMatrix();
  return *this;
}

MatrixColumn& MatrixColumn::operator+=(const MatrixColumn& d) {
  if (length != d.length) throw std::logic_error("Lengths do not match");
  
  for (int i = 0; i < length; i++)
    column[i] += d[i];
  
  setMatrix();
  return *this;
}

MatrixColumn operator+(const MatrixColumn& m, double d) {
  MatrixColumn r(m);
  return r += d;
}

MatrixColumn& MatrixColumn::operator-=(double d) {
  for (int i = 0; i < length; i++)
    column[i] -= d;
  
  setMatrix();
  return *this;
}

MatrixColumn& MatrixColumn::operator-=(const MatrixColumn& d) {
  if (length != d.length) throw std::logic_error("Lengths do not match");
  
  for (int i = 0; i < length; i++)
    column[i] -= d[i];
  
  setMatrix();
  return *this;
}

MatrixColumn operator-(const MatrixColumn& m, double d) {
  MatrixColumn r(m);
  return r -= d;
}

MatrixColumn& MatrixColumn::operator*=(double d) {
  for (int i = 0; i < length; i++)
    column[i] *= d;
  
  setMatrix();
  return *this;
}

MatrixColumn operator*(const MatrixColumn& m, double d) {
  MatrixColumn r(m);
  return r *= d;
}

MatrixColumn& MatrixColumn::operator/=(double d) {
  if (!d) throw std::runtime_error("Attempted to divide by zero");
  for (int i = 0; i < length; i++)
    column[i] /= d;
  
  setMatrix();
  return *this;
}

MatrixColumn operator/(const MatrixColumn& m, double d) {
  if (!d) throw std::runtime_error("Attempted to divide by zero");
  MatrixColumn r(m);
  return r /= d;
}

std::ostream& operator<<(std::ostream& os, const MatrixColumn& r) {
  std::cout << "     [";
  for(int i = 0; i < r.getLength(); i++)
    std::cout << std::setw(11) << r[i];
  std::cout << " ]" << std::endl;
  return os;
}
