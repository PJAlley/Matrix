#pragma once

#include <vector>
#include <complex>
#include "Matrix.h"

class ColumnVector;

class SquareMatrix : public Matrix {
  public:
    SquareMatrix();
    SquareMatrix(int);
    SquareMatrix(const Matrix&);
    SquareMatrix(Matrix&&);
    SquareMatrix(const SquareMatrix&);
    SquareMatrix(SquareMatrix&&);
    SquareMatrix& operator=(const SquareMatrix&);
    SquareMatrix& operator=(const Matrix&);
    SquareMatrix& operator=(SquareMatrix&&);
    inline const int sides() const { return side; }
    const SquareMatrix& get_inverse();
    virtual ~SquareMatrix();
    double get_determinant();
    double get_trace();
    void balance();
    ColumnVector solve(const ColumnVector&);
    Matrix solve(const Matrix &);
    const std::vector<double> get_characteristic_polynomial();
    void factor();
    Matrix get_eigenvectors();
    void printCoefficients() const;
    void printEigenvalues() const;
  protected:
    int side;
  private:
    double trace;
    SquareMatrix * inverse;
    std::vector<double> coeffs;   // Characteristic Polynomial Coefficients.
    std::vector<std::complex<double>> eigenvalues;
    void char_poly();
};

SquareMatrix operator/(SquareMatrix&, SquareMatrix&);

class IdentityMatrix : public SquareMatrix {
  public:
    IdentityMatrix();
    IdentityMatrix(int);
    virtual ~IdentityMatrix();
};
