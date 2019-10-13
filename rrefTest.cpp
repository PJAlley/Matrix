#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include "Matrix.h"
#include "SquareMatrix.h"
//#include "MatrixRow.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 21 - 10;
    }
  }
}

int main() {
  srand(time(NULL));
  //int size = 4;
  //SquareMatrix sm(size);
  //load_random_values(sm);

  int r;
  int c;
  std::cout << "Enter number of rows: ";
  std::cin >> r;
  std::cout << "Enter number of columns: ";
  std::cin >> c;
  Matrix rc(r, c);
  std::cin >> rc;

  std::cout << "Matrix:" << std::endl << rc;

  rc.reducedRowEchelonForm();

  std::cout << "Reduced Row Echelon Form:" << std::endl << rc;
  return 0;
}

//1 2 -1 -4 2 3 -1 -11 -2 0 -3 22

//3 -3 3 3 -3 3 6 -6 6
