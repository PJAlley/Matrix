#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <cmath>
#include <limits>
#include <algorithm>
#include "MatrixRow.h"
#include "Matrix.h"

MatrixRow::MatrixRow(int l) : row(new double [l]), length(l), alloc(true) {}

MatrixRow::MatrixRow(Matrix& m, int r) : row(m[r]), length(m.cols()), alloc(false) {}

MatrixRow::MatrixRow(double * d, int l) : row(d), length(l), alloc(false) {}

MatrixRow::MatrixRow(const MatrixRow& m) {
  alloc = true;
  length = m.length;
  row = new double[m.length];
  for (int i = 0; i < length; i++)
    row[i] = m[i];
}

MatrixRow& MatrixRow::operator=(const MatrixRow & m) {
  if (&m == this) return *this;
  if (alloc) delete [] row;

  alloc = true;
  length = m.length;
  row = new double[m.length];
  for (int i = 0; i < length; i++)
    row[i] = m[i];
  return *this;
}

MatrixRow::~MatrixRow() {
  if (alloc) delete [] row;
}

const double& MatrixRow::operator[](int r) const {
  if(r < 0 || r >= length) throw std::range_error("Subscript out of range");
  return row[r];
}

double& MatrixRow::operator[](int r) {
  return const_cast<double&>(static_cast<const MatrixRow&>(*this)[r]);
}

double MatrixRow::sum() const {
  double sum = 0;
  for (int i = 0; i < length; i++) 
    sum += row[i];
  return sum;
}

double MatrixRow::absSum() const {
  double sum = 0;
  for (int i = 0; i < length; i++) 
    sum += std::fabs(row[i]);
  return sum;
}

double MatrixRow::product() const {
  double product = 1;
  for (int i = 0; i < length; i++) 
    product *= row[i];
  return product;
}

MatrixRow& MatrixRow::operator+=(double d) {
  for (int i = 0; i < length; i++)
    row[i] += d;
  
  return *this;
}

MatrixRow& MatrixRow::operator+=(const MatrixRow& d) {
  if (length != d.length) throw std::logic_error("Lengths do not match");
  
  for (int i = 0; i < length; i++)
    row[i] += d[i];
  
  return *this;
}

MatrixRow operator+(const MatrixRow& m, double d) {
  MatrixRow r(m);
  return r += d;
}

MatrixRow& MatrixRow::operator-=(double d) {
  for (int i = 0; i < length; i++)
    row[i] -= d;
  
  return *this;
}

MatrixRow& MatrixRow::operator-=(const MatrixRow& d) {
  if (length != d.length) throw std::logic_error("Lengths do not match");
  
  for (int i = 0; i < length; i++)
    row[i] -= d[i];
  
  return *this;
}

MatrixRow operator-(const MatrixRow& m, double d) {
  MatrixRow r(m);
  return r -= d;
}

MatrixRow& MatrixRow::operator*=(double d) {
  for (int i = 0; i < length; i++)
    row[i] *= d;
  
  return *this;
}

MatrixRow operator*(const MatrixRow& m, double d) {
  MatrixRow r(m);
  return r *= d;
}

MatrixRow& MatrixRow::operator/=(double d) {
  if (!d) throw std::runtime_error("Attempted to divide by zero");
  for (int i = 0; i < length; i++)
    row[i] /= d;
  
  return *this;
}

MatrixRow operator/(const MatrixRow& m, double d) {
  if (!d) throw std::runtime_error("Attempted to divide by zero");
  MatrixRow r(m);
  return r /= d;
}

std::ostream& operator<<(std::ostream& os, const MatrixRow& r) {
  std::cout << "     [";
  for(int i = 0; i < r.getLength(); i++)
    std::cout << std::setw(11) << r[i];
  std::cout << " ]" << std::endl;
  return os;
}
