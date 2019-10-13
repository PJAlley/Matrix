#pragma once

class Matrix;

class MatrixColumn {
  public:
    MatrixColumn(int);
    MatrixColumn(Matrix&, int);
    MatrixColumn(double *, int);
    ~MatrixColumn();
    MatrixColumn() = delete;
    MatrixColumn(const MatrixColumn&);
    MatrixColumn& operator=(const MatrixColumn&);
    inline int getLength() const {return length;}
    double& operator[](int);
    const double& operator[](int) const;
    double sum() const;
    double absSum() const;
    double product() const;
    MatrixColumn& operator+=(double);            // m[2] += 3
    MatrixColumn& operator-=(double);            // m[2] -= 3
    MatrixColumn& operator*=(double);            // m[2] *= 3
    MatrixColumn& operator/=(double);            // m[2] /= 3
    MatrixColumn& operator+=(const MatrixColumn&);  // m[2] += m[3]
    MatrixColumn& operator-=(const MatrixColumn&);  // m[2] -= m[3]
    friend std::ostream& operator<<(std::ostream&, const MatrixColumn&);
  private:
    void setMatrix();
    double * column;
    Matrix * matrix;
    int length;
    bool alloc;
    int index;
};

MatrixColumn operator+(const MatrixColumn&, double);  // m[2] = m[3] + 4
MatrixColumn operator-(const MatrixColumn&, double);  // m[2] = m[3] - 4
MatrixColumn operator*(const MatrixColumn&, double);  // m[2] = m[3] * 4
MatrixColumn operator/(const MatrixColumn&, double);  // m[2] = m[3] / 4
