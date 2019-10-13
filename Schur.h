#pragma once

#include <complex>

class Schur {
  public:
    Schur(SquareMatrix&);
    ~Schur();
    Schur() = delete;
    Schur(const Schur&) = delete;
    Schur& operator=(const Schur&) = delete;
    void compute();
    void Hessenberg();
    void HQR_Decomposition();
    void printReal() const;
    void printImaginary() const;
    void printEigenvalues() const;
    inline const SquareMatrix * get_matrix() { return matrix; }
    inline SquareMatrix * get_eigenvector() { return evector; }
    std::vector<std::complex<double>> get_eigenvalues() {return eigenvalues; }
    SquareMatrix get_diagonal_matrix();
  private:
    SquareMatrix * matrix;
    SquareMatrix * evector;
    int side;
    std::vector<std::complex<double>> eigenvalues;
};
