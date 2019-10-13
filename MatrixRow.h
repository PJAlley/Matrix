#pragma once

class Matrix;

class MatrixRow {
  public:
    MatrixRow(int);
    MatrixRow(Matrix&, int);
    MatrixRow(double *, int);
    ~MatrixRow();
    MatrixRow() = delete;
    MatrixRow(const MatrixRow&);
    MatrixRow& operator=(const MatrixRow&);
    inline int getLength() const {return length;}
    double& operator[](int);
    const double& operator[](int) const;
    double sum() const;
    double absSum() const;
    double product() const;
    MatrixRow& operator+=(double);            // m[2] += 3
    MatrixRow& operator-=(double);            // m[2] -= 3
    MatrixRow& operator*=(double);            // m[2] *= 3
    MatrixRow& operator/=(double);            // m[2] /= 3
    MatrixRow& operator+=(const MatrixRow&);  // m[2] += m[3]
    MatrixRow& operator-=(const MatrixRow&);  // m[2] -= m[3]
    friend std::ostream& operator<<(std::ostream&, const MatrixRow&);
  private:
    double * row;
    int length;
    bool alloc;
};

MatrixRow operator+(const MatrixRow&, double);  // m[2] = m[3] + 4
MatrixRow operator-(const MatrixRow&, double);  // m[2] = m[3] - 4
MatrixRow operator*(const MatrixRow&, double);  // m[2] = m[3] * 4
MatrixRow operator/(const MatrixRow&, double);  // m[2] = m[3] / 4
