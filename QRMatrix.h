#pragma once

class ColumnVector;

class QRMatrix {
  public:
    QRMatrix(SquareMatrix&);
    ~QRMatrix();
    QRMatrix() = delete;
    QRMatrix(const QRMatrix&) = delete;
    QRMatrix& operator=(const QRMatrix&) = delete;
    void compute();
    void Householder();
    inline const SquareMatrix * get_q() const { return q; }
    inline const SquareMatrix * get_qt() const { return qt; }
    inline const SquareMatrix * get_r() const { return r; }
    ColumnVector solve(const ColumnVector&);
    Matrix solve(const Matrix &);
  private:
    int side;
    SquareMatrix * q;
    SquareMatrix * qt;
    SquareMatrix * r;
};
