#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include "Matrix.h"
#include "RowVector.h"

RowVector::RowVector() : length(0) {}

RowVector::RowVector(int r) : row(r), length(r) {}

// Get the r-th row of the matrix and load it in the vector.
// Length is the number of columns in the matrix.
RowVector::RowVector(Matrix& m, int r) : length(m.cols()) {
  if (r >= length) throw std::range_error("Row number out of range.");
  double * m_row = m[r];
  for (int i = 0; i < length; i++)
    row.push_back(m_row[i]);
}

// Input should be a pointer to an array. Length of array MUST be provided.
RowVector::RowVector(double * el, int l) {
  length = l;
  for (int i = 0; i < length; i++)
    row.push_back(el[i]);
}

void RowVector::set_size(int r) {
  if (r <= 0) throw std::range_error("Length cannot be negative.");
  row.resize(r);
  length = r;
}

const double& RowVector::operator[](int r) const {
  if (r < 0 || r >= row.size()) throw std::range_error("Subscript out of range");
  return row[r];
}

double& RowVector::operator[](int r) {
  return const_cast<double&>(static_cast<const RowVector&>(*this)[r]);
}

std::ostream& operator<<(std::ostream& os, const RowVector& r) {
  std::cout << "     [";
  for(const auto i: r.row)
    std::cout << std::setw(11) << i;
  std::cout << " ]" << std::endl;
  return os;
}
