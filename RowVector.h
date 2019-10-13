#pragma once

class RowVector {
  public:
    RowVector();
    RowVector(int);
    RowVector(Matrix&, int);
    RowVector(double *, int);
    void set_size(int);
    inline int size() const { return length; }
    double& operator[](int);
    const double& operator[](int) const ;
    friend std::ostream& operator<<(std::ostream&, const RowVector&);
  private:
    std::vector<double> row;
    int length;
};
