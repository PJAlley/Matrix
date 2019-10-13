#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
#include <limits>
#include "SquareMatrix.h"
//#include "Matrix.h"
#include "Schur.h"


Schur::Schur(SquareMatrix& s) : matrix(new SquareMatrix(s)), evector(new IdentityMatrix(s.sides())), 
  side(s.sides()), eigenvalues(s.sides()) { compute(); }

Schur::~Schur() {
  delete matrix;
  delete evector;
}

void Schur::printReal() const {
  std::cout << "Real:"<< std::endl;
  std::cout << "     [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << real(i);
  std::cout << " ]" << std::endl;
}

void Schur::printImaginary() const {
  std::cout << "Imaginary:"<< std::endl;
  std::cout << "     [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << imag(i);
  std::cout << " ]" << std::endl;
}

void Schur::printEigenvalues() const {
  std::cout << "Eigenvalues:"<< std::endl;
  std::cout << std::setw(6) << "R [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << real(i);
  std::cout << " ]" << std::endl;

  std::cout << std::setw(6) << "I [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << imag(i);
  std::cout << " ]" << std::endl;
}

SquareMatrix Schur::get_diagonal_matrix () {
  SquareMatrix diag(side);

  for (int i = 0; i < side; i++) {
    diag(i, i) = real(eigenvalues[i]);
    if (imag(eigenvalues[i]) > 0) {
      diag(i, i + 1) = imag(eigenvalues[i]);
    }
    else if (imag(eigenvalues[i]) < 0) {
      diag(i, i - 1) = imag(eigenvalues[i]);
    }
  }
  return diag;
}

void Schur::compute() {
  Hessenberg();
  HQR_Decomposition();
}

void Schur::Hessenberg() {
  using std::fabs;

  int low = 0;
  int high = side - 1;

  std::vector<double> ort(side);

  for (int m = low + 1; m <= high - 1; m++) {
    // Scale column.
    double scale = 0;
    for (int i = m; i <= high; i++) {
      scale += fabs(matrix->at(i, m - 1));
    }

    if (scale) {
      // ompute Householder transformation.
      double h = 0;
      for (int i = high; i >= m; i--) {
        ort[i] = matrix->at(i, m - 1) / scale;
        h += ort[i] * ort[i];
      }
      double g = std::sqrt(h);
      if (ort[m] > 0) g = -g;
      h -= ort[m] * g;
      ort[m] -= g;

      // Apply Householder similarity transformation

      for (int j = m; j < side; j++) {
        double f = 0;
        for (int i = high; i >= m; i--) {
          f += ort[i] * matrix->at(i, j);
        }
        f /= h;
        for (int i = m; i <= high; i++) {
          matrix->at(i, j) -= f * ort[i];
        }
      }

      for (int i = 0; i <= high; i++) {
        double f = 0;
        for (int j = high; j >= m; j--) {
          f += ort[j] * matrix->at(i, j);
        }
        f /= h;
        for (int j = m; j <= high; j++) {
          matrix->at(i, j) -= f * ort[j];
        }
      }
      ort[m] = scale * ort[m];
      matrix->at(m, m - 1) = scale * g;
    }
  }

  // Accumulate transformations (Algol's ortran).
  for (int m = high - 1; m >= low + 1; m--) {
    if (matrix->at(m, m-1) != 0) {
      for (int i = m + 1; i <= high; i++) {
        ort[i] = matrix->at(i, m - 1);
      }
      for (int j = m; j <= high; j++) {
        double g = 0;
        for (int i = m; i <= high; i++) {
          g += ort[i] * evector->at(i, j);
        }
        // Double division avoids possible underflow
        g = (g / ort[m]) / matrix->at(m, m - 1);
        for (int i = m; i <= high; i++) {
          evector->at(i, j) += g * ort[i];
        }
      }
    }
  }
  for (int i = 0; i < side; i++) {
    for (int j = 0; j < i; ++j) {
      if (j < i - 1) matrix->at(i, j) = 0;
    }
  }
}

void Schur::HQR_Decomposition() {
  using std::pow;
  using std::min;
  using std::fmax;
  using std::fabs;
  using std::sqrt;
  using std::complex;

  //  This is derived from the Algol procedure hqr2,
  //  by Martin and Wilkinson, Handbook for Auto. Comp.,
  //  Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.

  // Initialize

  int n = side - 1;
  int low = 0;
  int high = side - 1;
  double eps = pow(2.0, -52.0);
  double shift = 0.0;
  double p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;

  // Store roots isolated by balanc and compute matrix norm

  double norm = 0.0;
  for (int i = 0; i < side; i++) {
    if ((i < low) || (i > high)) {
      eigenvalues[i] = matrix->at(i, i);
    }
    for (int j = std::max(i - 1, 0); j < side; j++) {
      norm = norm + fabs(matrix->at(i, j));
    }
  }

  // Outer loop over eigenvalue index

  int iter = 0;
  while (n >= low) {
    // Look for single small sub-diagonal element

    int l = n;
    while (l > low) {
      s = fabs(matrix->at(l - 1, l - 1)) + fabs(matrix->at(l, l));
      if (s == 0.0) {
        s = norm;
      }
      if (fabs(matrix->at(l, l - 1)) < eps * s) {
        break;
      }
      l--;
    }
   
    // Check for convergence
    // One root found

    if (l == n) {
      matrix->at(n, n) = matrix->at(n, n) + shift;
      eigenvalues[n] = matrix->at(n, n);
      n--;
      iter = 0;
    } 
    else if (l == n - 1) {
      // Two roots found
      w = matrix->at(n, n - 1) * matrix->at(n - 1, n);
      p = (matrix->at(n - 1, n - 1) - matrix->at(n, n)) / 2.0;
      q = p * p + w;
      z = sqrt(fabs(q));
      matrix->at(n, n) = matrix->at(n, n) + shift;
      matrix->at(n - 1, n - 1) = matrix->at(n - 1, n - 1) + shift;
      x = matrix->at(n, n);

      // Real pair
      if (q >= 0) {
        if (p >= 0) {
          z = p + z;
        } 
        else {
          z = p - z;
        }
        eigenvalues[n - 1] = x + z;
        eigenvalues[n] = eigenvalues[n - 1];
        if (z != 0.0) {
          eigenvalues[n] = x - w / z;
        }
        x = matrix->at(n, n - 1);
        s = fabs(x) + fabs(z);
        p = x / s;
        q = z / s;
        r = sqrt(p * p + q * q);
        p = p / r;
        q = q / r;

        // Row modification
        for (int j = n - 1; j < side; j++) {
          z = matrix->at(n - 1, j);
          matrix->at(n - 1, j) = q * z + p * matrix->at(n, j);
          matrix->at(n, j) = q * matrix->at(n, j) - p * z;
        }

        // Column modification
        for (int i = 0; i <= n; i++) {
          z = matrix->at(i, n - 1);
          matrix->at(i, n - 1) = q * z + p * matrix->at(i, n);
          matrix->at(i, n) = q * matrix->at(i, n) - p * z;
        }

        // Accumulate transformations
        for (int i = low; i <= high; i++) {
          z = evector->at(i, n - 1);
          evector->at(i, n - 1) = q * z + p * evector->at(i, n);
          evector->at(i, n) = q * evector->at(i, n) - p * z;
        }

      } 
      else {
        // Complex pair
        eigenvalues[n] = complex(x + p, -z);
        eigenvalues[n - 1] = conj(eigenvalues[n]);
      }
      n = n - 2;
      iter = 0;

     // No convergence yet

    } 
    else {
      // No roots yet. Form shift.
      x = matrix->at(n, n);
      y = 0.0;
      w = 0.0;
      if (l < n) {
        y = matrix->at(n - 1, n - 1);
        w = matrix->at(n, n - 1) * matrix->at(n - 1, n);
      }

      // Wilkinson's original ad hoc shift
      if (iter == 10) {
        shift += x;
        for (int i = low; i <= n; i++) {
          matrix->at(i, i) -= x;
        }
        s = fabs(matrix->at(n, n - 1)) + fabs(matrix->at(n - 1, n - 2));
        x = y = 0.75 * s;
        w = -0.4375 * s * s;
      }

      // MATLAB's new ad hoc shift
      if (iter == 30) {
        s = (y - x) / 2.0;
        s = s * s + w;
        if (s > 0) {
          s = sqrt(s);
          if (y < x) {
            s = -s;
          }
          s = x - w / ((y - x) / 2.0 + s);
          for (int i = low; i <= n; i++) {
            matrix->at(i, i) -= s;
          }
          shift += s;
          x = y = w = 0.964;
        }
      }

      iter = iter + 1;   // (Could check iteration count here.)

      // Look for two consecutive small sub-diagonal elements
      int m = n - 2;
      while (m >= l) {
        z = matrix->at(m, m);
        r = x - z;
        s = y - z;
        p = (r * s - w) / matrix->at(m + 1, m) + matrix->at(m, m + 1);
        q = matrix->at(m + 1, m + 1) - z - r - s;
        r = matrix->at(m + 2, m + 1);
        s = fabs(p) + fabs(q) + fabs(r);
        p = p / s;
        q = q / s;
        r = r / s;
        if (m == l) {
          break;
        }
        if (fabs(matrix->at(m, m - 1)) * (fabs(q) + fabs(r)) <
          eps * (fabs(p) * (fabs(matrix->at(m - 1, m - 1)) + fabs(z) + fabs(matrix->at(m + 1, m + 1))))) {
            break;
          }
        m--;
      }

      for (int i = m + 2; i <= n; i++) {
        matrix->at(i, i - 2) = 0.0;
        if (i > m + 2) {
          matrix->at(i, i - 3) = 0.0;
        }
      }

      // Double QR step involving rows l:n and columns m:n
      for (int k = m; k <= n - 1; k++) {
        int notlast = (k != n - 1);
        if (k != m) {
          p = matrix->at(k, k - 1);
          q = matrix->at(k + 1, k - 1);
          r = (notlast ? matrix->at(k + 2, k - 1) : 0.0);
          x = fabs(p) + fabs(q) + fabs(r);
          if (x != 0.0) {
            p = p / x;
            q = q / x;
            r = r / x;
          }
        }
        if (x == 0.0) {
          break;
        }
        s = sqrt(p * p + q * q + r * r);
        if (p < 0) {
          s = -s;
        }
        if (s != 0) {
          if (k != m) {
            matrix->at(k, k - 1) = -s * x;
          } 
          else if (l != m) {
            matrix->at(k, k - 1) = -matrix->at(k, k - 1);
          }
          p = p + s;
          x = p / s;
          y = q / s;
          z = r / s;
          q = q / p;
          r = r / p;

          // Row modification
          for (int j = k; j < side; j++) {
            p = matrix->at(k, j) + q * matrix->at(k + 1, j);
            if (notlast) {
              p = p + r * matrix->at(k + 2, j);
              matrix->at(k + 2, j) = matrix->at(k + 2, j) - p * z;
            }
            matrix->at(k, j) = matrix->at(k, j) - p * x;
            matrix->at(k + 1, j) = matrix->at(k + 1, j) - p * y;
          }

          // Column modification
          for (int i = 0; i <= min(n, k + 3); i++) {
            p = x * matrix->at(i, k) + y * matrix->at(i, k + 1);
            if (notlast) {
              p = p + z * matrix->at(i, k + 2);
              matrix->at(i, k + 2) = matrix->at(i, k + 2) - p * r;
            }
            matrix->at(i, k) = matrix->at(i, k) - p;
            matrix->at(i, k + 1) = matrix->at(i, k + 1) - p * q;
          }

          // Accumulate transformations
          for (int i = low; i <= high; i++) {
            p = x * evector->at(i, k) + y * evector->at(i, k + 1);
            if (notlast) {
              p = p + z * evector->at(i, k + 2);
              evector->at(i, k + 2) = evector->at(i, k + 2) - p * r;
            }
            evector->at(i, k) = evector->at(i, k) - p;
            evector->at(i, k + 1) = evector->at(i, k + 1) - p * q;
          }
        }  // (s != 0)
      }  // k loop
    }  // check convergence
  }  // while (n >= low)
  
  // All roots found. Backsubstitute to find vectors of upper triangular form.
  if (norm == 0.0) {
    return;
  }

  for (n = side - 1; n >= 0; n--) {
    p = real(eigenvalues[n]);
    q = imag(eigenvalues[n]);
    
    if (q == 0) {
      // Real vector
      int l = n;
      matrix->at(n, n) = 1.0;
      for (int i = n - 1; i >= 0; i--) {
        w = matrix->at(i, i) - p;
        r = 0.0;
        for (int j = l; j <= n; j++) {
          r = r + matrix->at(i, j) * matrix->at(j, n);
        }
        if (imag(eigenvalues[i]) < 0.0) {
          z = w;
          s = r;
        } 
        else {
          l = i;
          if (imag(eigenvalues[i]) == 0.0) {
            if (w != 0.0) {
              matrix->at(i, n) = -r / w;
            } 
            else {
              matrix->at(i, n) = -r / (eps * norm);
            }
          } 
          else {
            // Solve real equations
            x = matrix->at(i, i + 1);
            y = matrix->at(i + 1, i);
            q = (real(eigenvalues[i]) - p) * (real(eigenvalues[i]) - p) + imag(eigenvalues[i]) * imag(eigenvalues[i]);
            t = (x * s - z * r) / q;
            matrix->at(i, n) = t;
            if (fabs(x) > fabs(z)) {
              matrix->at(i + 1, n) = (-r - w * t) / x;
            } 
            else {
              matrix->at(i + 1, n) = (-s - y * t) / z;
            }
          }

          // Overflow control
          t = fabs(matrix->at(i, n));
          if ((eps * t) * t > 1) {
            for (int j = i; j <= n; j++) {
              matrix->at(j, n) = matrix->at(j, n) / t;
            }
          }
        }
      }
    } 
    else if (q < 0) {
      // Complex vector
      int l = n - 1;

      // Last vector component imaginary so matrix is triangular
      if (fabs(matrix->at(n, n - 1)) > fabs(matrix->at(n - 1, n))) {
        matrix->at(n - 1, n - 1) = q / matrix->at(n, n - 1);
        matrix->at(n - 1, n) = -(matrix->at(n, n) - p) / matrix->at(n, n - 1);
      } 
      else {
        complex<double> temp = complex<double>(0, -matrix->at(n - 1, n)) / complex<double>(matrix->at(n - 1, n - 1) - p, q);
        matrix->at(n - 1, n - 1) = real(temp);
        matrix->at(n - 1, n) = imag(temp);
      }
      matrix->at(n, n - 1) = 0.0;
      matrix->at(n, n) = 1.0;
      for (int i = n - 2; i >= 0; i--) {
        double ra, sa, vr, vi;
        ra = 0.0;
        sa = 0.0;
        for (int j = l; j <= n; j++) {
          ra = ra + matrix->at(i, j) * matrix->at(j, n - 1);
          sa = sa + matrix->at(i, j) * matrix->at(j, n);
        }
        w = matrix->at(i, i) - p;

        if (imag(eigenvalues[i]) < 0.0) {
          z = w;
          r = ra;
          s = sa;
        } 
        else {
          l = i;
          if (imag(eigenvalues[i]) == 0) {
            complex<double> temp = complex(-ra, -sa) / complex(w, q);
            matrix->at(i, n - 1) = real(temp);
            matrix->at(i, n) = imag(temp);
          } 
          else {
            // Solve complex equations
            x = matrix->at(i, i + 1);
            y = matrix->at(i + 1, i);
            vr = (real(eigenvalues[i]) - p) * (real(eigenvalues[i]) - p) + imag(eigenvalues[i]) * imag(eigenvalues[i]) - q * q;
            vi = (real(eigenvalues[i]) - p) * 2.0 * q;
            if ((vr == 0.0) && (vi == 0.0)) {
              vr = eps * norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z));
            }
            complex<double> temp = complex(x * r - z * ra + q * sa, x * s - z * sa - q * ra) / complex(vr, vi);
            matrix->at(i, n - 1) = real(temp);
            matrix->at(i, n) = imag(temp);
            if (fabs(x) > (fabs(z) + fabs(q))) {
              matrix->at(i + 1, n - 1) = (-ra - w * matrix->at(i, n - 1) + q * matrix->at(i, n)) / x;
              matrix->at(i + 1, n) = (-sa - w * matrix->at(i, n) - q * matrix->at(i, n - 1)) / x;
            } 
            else {
              complex<double> temp = complex(-r - y * matrix->at(i, n - 1), -s - y * matrix->at(i, n)) / complex(z, q);
              matrix->at(i + 1, n - 1) = real(temp);
              matrix->at(i + 1, n) = imag(temp);
            }
          }

          // Overflow control
          t = fmax(fabs(matrix->at(i, n - 1)), fabs(matrix->at(i, n)));
          if ((eps * t) * t > 1) {
            for (int j = i; j <= n; j++) {
              matrix->at(j, n - 1) = matrix->at(j, n - 1) / t;
              matrix->at(j, n) = matrix->at(j, n) / t;
            }
          }
        }
      }
    }
  }

  // Vectors of isolated roots
  for (int i = 0; i < side; i++) {
    if (i < low || i > high) {
      for (int j = i; j < side; j++) {
        evector->at(i, j) = matrix->at(i, j);
      }
    }
  }

  // Back transformation to get eigenvectors of original matrix
  for (int j = side - 1; j >= low; j--) {
    for (int i = low; i <= high; i++) {
      z = 0.0;
      for (int k = low; k <= min(j, high); k++) {
        z = z + evector->at(i, k) * matrix->at(k, j);
      }
      evector->at(i, j) = z;
    }
  }
}
