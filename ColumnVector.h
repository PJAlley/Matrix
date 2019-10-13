#pragma once

#include <vector>

class RowVector;

class ColumnVector {
  public:
    ColumnVector();
    ColumnVector(int);
    ColumnVector(const ColumnVector&);
    ColumnVector(const Matrix&, int);
    void set_size(int);
    inline int size() const { return length; }
    double& operator[](int);
    const double& operator[](int) const ;
    friend std::ostream& operator<<(std::ostream&, const ColumnVector&);
    friend std::istream& operator>>(std::istream&, ColumnVector&);
  private:
    std::vector<double> col;
    int length;
};

ColumnVector operator*(const Matrix&, const ColumnVector&);
double dotProduct(const RowVector&, const ColumnVector&);
Matrix outerProduct(const ColumnVector&, const ColumnVector&);
Matrix outerProduct(const ColumnVector&, const RowVector&);

