#pragma once

class ColumnVector;

class LUPMatrix {
  public:
    LUPMatrix(SquareMatrix&);
    ~LUPMatrix();
    LUPMatrix() = delete;
    LUPMatrix(const LUPMatrix&) = delete;
    LUPMatrix& operator=(const LUPMatrix&) = delete;
    void compute();
    double determinant();
    inline const SquareMatrix * get_lower() {return lower;}
    inline const SquareMatrix * get_upper() {return upper;}
    inline const SquareMatrix * get_permutation_matrix() {return p;}
    ColumnVector solve(const ColumnVector&);
    Matrix solve(const Matrix &);
  private:
    int side;
    int sign;
    SquareMatrix * m;
    SquareMatrix * lower;
    SquareMatrix * upper;
    IdentityMatrix * p;
};
