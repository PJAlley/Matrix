#pragma once

class ColumnVector;

class LUMatrix {
  public:
    LUMatrix(SquareMatrix&);
    ~LUMatrix();
    LUMatrix() = delete;
    LUMatrix(const LUMatrix&) = delete;
    LUMatrix& operator=(const LUMatrix&) = delete;
    void compute();
    double determinant();
    inline const SquareMatrix * get_lower() {return lower;}
    inline const SquareMatrix * get_upper() {return upper;}
    ColumnVector solve(const ColumnVector&);
    Matrix solve(const Matrix &);
  private:
    int side;
    SquareMatrix * m;
    SquareMatrix * lower;
    SquareMatrix * upper;
};
