#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
//#include <omp.h>
#include <chrono>
#include "Matrix.h"

void load_random_values(Matrix& matrix) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      matrix(i, j) = rand() % 21 - 10;
    }
  }
}

int main() {
  srand(time(NULL));
  std::cout << "Testing Matrix..." << std::endl;
  
  SquareMatrix sm(10);
  load_random_values(sm);
  std::cout << "sm:" << std::endl << sm;
  
  //double start = omp_get_wtime();
  auto start = std::chrono::high_resolution_clock::now();
  double det = sm.get_determinant();
  auto end = std::chrono::high_resolution_clock::now();
  //double end = omp_get_wtime();
  
  std::cout << "Determinant: " << det << std::endl;
  //double diff = end - start;
  auto elapsed = end - start;
  //std::cout << "Took " << diff << std::endl;
  std::cout << "Elapsed time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";

  // LU Decomposition
  //start = omp_get_wtime();
  //start = std::chrono::high_resolution_clock::now();
  //det = sm.get_det();
  //end = std::chrono::high_resolution_clock::now();
  //end = omp_get_wtime();
  //std::cout << "Determinant: " << det << std::endl;
  //diff = end - start;
  //auto elapsed2 = end - start;
  //std::cout << "Took " << diff << std::endl;
  //std::cout << "Elapsed time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
  std::cout << "Upper: " << std::endl << sm.LUupper();
  std::cout << "Lower: " << std::endl << sm.LUlower();

  return 0;
}