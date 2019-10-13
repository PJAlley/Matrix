#pragma once

#include <complex>
#include "SquareMatrix.h"

class Eigen {
  public:
    Eigen(SquareMatrix&);
    ~Eigen();
    Eigen() = delete;
    Eigen(const Eigen&) = delete;
    Eigen& operator=(const Eigen&) = delete;
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
    void transform();
    void balance();
    void back_balance();
    void sort();
    SquareMatrix * matrix;
    SquareMatrix * evector;
    int side;
    std::vector<std::complex<double>> eigenvalues;
    std::vector<double> scale;  // Scaling from balance.
    std::vector<int> perm;      // Permutation from Hessenberg.
};
