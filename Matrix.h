#pragma once

class SquareMatrix;
class ColumnVector;

class Matrix {
  friend class SquareMatrix;
  public:
    Matrix();
    Matrix(int, int);
    virtual ~Matrix();
    inline int rows() const { return row; }
    inline int cols() const { return col; }
    Matrix(const Matrix&);
    Matrix(Matrix&&);
    virtual Matrix& operator=(const Matrix&);
    virtual Matrix& operator=(Matrix&&);
    Matrix transpose();
    Matrix& operator+=(const Matrix&);
    Matrix& operator+=(double);
    Matrix& operator-=(const Matrix&);
    Matrix& operator-=(double);
    Matrix& operator*=(double);
    Matrix& operator/=(double);
    Matrix& swap(Matrix&);
    void swapRows(int, int);
    void swapColumns(int, int);
    double * operator[](int);
    const double * operator[](int) const ;
    double& operator()(int, int);
    const double& operator()(int, int) const;
    double& at(int, int);
    const double& at(int, int) const;
    void reducedRowEchelonForm();
    ColumnVector gauss(ColumnVector&);
    friend std::ostream& operator<<(std::ostream&, const Matrix&);
    friend std::istream& operator>>(std::istream&, Matrix&);
  protected:
    int row;
    int col;
    double ** matrix;
    void clear();
    void check(int, int) const;
    void check(int) const;
};

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator+(const Matrix&, double);
Matrix operator+(double, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, double);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Matrix operator*(Matrix&, Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
double dotProduct(const Matrix&, const Matrix&);
Matrix operator/(const Matrix&, double);

void swap_elements(double&, double&);
