#pragma once

#include <vector>

class ColumnVector;

// Single Value Decomposition.
// Decompose matrix M with U*W*Vt.
class SVD {
  public:
    SVD(Matrix &);
    SVD() = delete;
    SVD(const SVD&) = delete;
    SVD& operator=(const SVD&) = delete;
    ~SVD();
    void compute();
    inline const Matrix * get_u() { return u; }
    inline const Matrix * get_v() { return v; }
    inline const Matrix get_vt() { return v->transpose(); }
    const SquareMatrix* get_w();
    int rank();
    int nullity();
    Matrix range();
    Matrix nullspace();
    void printW() const;
    ColumnVector solve(const ColumnVector&);
    Matrix solve(const Matrix &);
  private:
    void decompose();
    void reorder();
    Matrix * u;
    SquareMatrix * v;
    std::vector<double> w;
    int row;
    int col;
};