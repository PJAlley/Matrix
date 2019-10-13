#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
#include <limits>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "Eigen.h"

// Returns a with the sign of b.
static inline double sign(double a, double b) {
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

Eigen::Eigen(SquareMatrix& s) : matrix(new SquareMatrix(s)), evector(new IdentityMatrix(s.sides())),
  side(s.sides()), eigenvalues(s.sides()), scale(s.sides(), 1), perm(s.sides()) { compute(); }

Eigen::~Eigen() {
  delete matrix;
  delete evector;
}

void Eigen::printReal() const {
  std::cout << "Real:"<< std::endl;
  std::cout << "     [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << real(i);
  std::cout << " ]" << std::endl;
}

void Eigen::printImaginary() const {
  std::cout << "Imaginary:"<< std::endl;
  std::cout << "     [";
  for(const auto i: eigenvalues)
    std::cout << std::setw(11) << imag(i);
  std::cout << " ]" << std::endl;
}

void Eigen::printEigenvalues() const {
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

SquareMatrix Eigen::get_diagonal_matrix () {
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

void Eigen::compute() {
  balance();
  Hessenberg();
  transform();
  HQR_Decomposition();
  back_balance();
  sort();
}

void Eigen::balance() {
  const double radix = std::numeric_limits<double>::radix;
  double sqrRadix = radix * radix;
  bool done = false;

  while (!done) {
    done = true;
    for (auto i = 0; i < side; ++i) {
      double r = 0;
      double c = 0;
      for (auto j = 0; j < side; ++j) {
        if (j != i) {
          r += std::abs(matrix->at(i, j));
          c += std::abs(matrix->at(j, i));
        }
      }
      if (c != 0 && r != 0) {
        double g = r / radix;
        double f = 1;
        double s = c + r;
        while (c < g) {
          f *= radix;
          c *= sqrRadix;
        }
        g = r * radix;
        while (c > g) {
          f /= radix;
          c /= sqrRadix;
        }

        if ((c + r) / f < 0.95 * s) {
          done = false;
          g = 1 / f;
          scale[i] *= f;
          for (auto j = 0; j < side; ++j) {
            matrix->at(i, j) *= g;
            matrix->at(j, i) *= f;
          }
        }
      }
    }
  }
}

void Eigen::back_balance() {
  for (int i = 0; i < side; i++)
    for (int j = 0; j < side; j++)
      evector->at(i, j) *= scale[i];
}

void Eigen::Hessenberg() {
  for (auto r = 1; r < side - 1; ++r) {
    double pivot = 0;
    int i = r;
    for (auto j = r; j < side; ++j) {
      double x = matrix->at(j, r - 1);
      if (std::abs(x) > std::abs(pivot)) {
        pivot = x;
        i = j;
      }
    }
    perm[r] = i;
    if (i != r) {
      for (auto j = r - 1; j < side; ++j)
        swap_elements(matrix->at(i, j), matrix->at(r, j));
      matrix->swapColumns(i, r);
    }

    if (pivot) {
      for (auto i = r + 1; i < side; ++i) {
        double y = matrix->at(i, r - 1);
        if (y) {
          y /= pivot;
          matrix->at(i, r - 1) = y;
          for (auto j = r; j < side; ++j) matrix->at(i, j) -= y * matrix->at(r, j);
          for (auto j = 0; j < side; ++j) matrix->at(j, r) += y * matrix->at(j, i);
        }
      }
    }
  }
  
  for (auto i = 0; i < side; i++) {
    for (auto j = 0; j < i; ++j) {
      if (j < i - 1) matrix->at(i, j) = 0;
    }
  }
  
}

void Eigen::transform() {
  for (auto p = side - 2; p > 0; --p) {
    for (auto k = p + 1; k < side; ++k)
      evector->at(k, p) = matrix->at(k, p - 1);
    auto i = perm[p];
    if (i != p) {
      for (auto j = p; j < side; ++j) {
        evector->at(p, j) = evector->at(i, j);
        evector->at(i, j) = 0;
      }
      evector->at(i, p) = 1;
    }
  }
}

void Eigen::sort() {
  using std::complex;
  std::vector<double> temp(side);

  int i;
  for (int j = 1; j < side; j++) {
    complex<double> x = eigenvalues[j];
    for (int k = 0; k < side; k++) {
      temp[k] = evector->at(k, j);
    }
    for (i = j - 1; i >= 0; i--) {
      if (real(eigenvalues[i]) >= real(x)) break;
      eigenvalues[i + 1] = eigenvalues[i];
      for(int k = 0; k < side; k++) {
        evector->at(k, i + 1) = evector->at(k, i);
      }
    }
    eigenvalues[i + 1] = x;
    for (int k = 0; k < side; k++)
      evector->at(k, i + 1) = temp[k];
  }
}

void Eigen::HQR_Decomposition() {
  using std::pow;
  using std::fmin;
  using std::fmax;
  using std::abs;
  using std::sqrt;
  using std::complex;

  const double epsilon = std::numeric_limits<double>::epsilon();

  // Compute matrix norm. Used in locating a small diagonal element.
  double norm = 0;
  for (auto i = 0; i < side; ++i) {
    for (auto j = fmax(i - 1, 0); j < side; ++j)
      norm += abs(matrix->at(i, j));
  }

  double p = 0.0, q = 0.0, r = 0.0, s = 0.0;
  double w = 0.0, x = 0.0, y = 0.0, z = 0.0;

  int nn = side - 1;
  int l;

  //  Gets changed only by an exceptional shift.
  double shift = 0;

  // Outer loop over eigenvalue index.
  while (nn >= 0) {     // Begin search for next eigenvalue.
    // Begin iteration: Look for single small sub-diagonal element
    int iter = 0;
    do {
      for (l = nn; l > 0; l--) {
        s = abs(matrix->at(l - 1, l - 1)) + abs(matrix->at(l, l));
        if (!s) s = norm;
        if (abs(matrix->at(l, l - 1)) < epsilon * s) {
          matrix->at(l, l - 1) = 0;
          break;
        }
      }

      x = matrix->at(nn, nn);

      if (l == nn) {
        // One root found.
        eigenvalues[nn] = matrix->at(nn, nn) = x + shift;
        nn--;
        iter = 0;
      }
      else {
        y = matrix->at(nn - 1, nn - 1);
        w = matrix->at(nn, nn - 1) * matrix->at(nn - 1, nn);
        if (l == nn - 1) {
          // Two roots found.
          p = 0.5 * (y - x);
          q = p * p + w;
          z = sqrt(abs(q));
          x += shift;
          matrix->at(nn, nn) = x;
          matrix->at(nn - 1, nn - 1) = y + shift;
          if (q >= 0) {
            // A real pair
            z = p + sign(z, p);
            eigenvalues[nn - 1] = eigenvalues[nn] = x + z;
            if (z) eigenvalues[nn] = x - w / z;
            x = matrix->at(nn, nn - 1);
            s = abs(x) + abs(z);
            p = x / s;
						q = z / s;
            r = sqrt(p * p + q * q);
            p /= r;
            q /= r;

            // Row modification.
            for (auto rr = nn - 1; rr < side; ++rr) {
              z = matrix->at(nn - 1, rr);
              matrix->at(nn - 1, rr) = q * z + p * matrix->at(nn, rr);
              matrix->at(nn, rr) = q * matrix->at(nn, rr) - p * z;
            }
            // Column modification.
            for (auto cc = 0; cc <= nn; ++cc) {
              z = matrix->at(cc, nn - 1);
              matrix->at(cc, nn - 1) = q * z + p * matrix->at(cc, nn);
              matrix->at(cc, nn) = q * matrix->at(cc, nn) - p * z;
            }
            // Accumulate transformations
            for (auto ii = 0; ii < side; ++ii) {
              z = evector->at(ii, nn - 1);
              evector->at(ii, nn - 1) = q * z + p * evector->at(ii, nn);
              evector->at(ii, nn) = q * evector->at(ii, nn) - p * z;
            }
          }
          else {
            // Complex pair.
            eigenvalues[nn] = complex(x + p, -z);
            eigenvalues[nn - 1] = conj(eigenvalues[nn]);
          }
          nn -= 2;
          iter = 0;
        }
        else {
          // No roots yet. Shift.
          // Exceptional shift at 10x. (Wilkinson's.)
          if (iter) {
            if (iter == 10) {
              shift += x;
              for (auto ii = 0; ii < nn; ++ii) matrix->at(ii, ii) -= x;
              s = abs(matrix->at(nn, nn - 1)) + abs(matrix->at(nn - 1, nn - 2));
              x = y = 0.75 * s;
              w = -0.4375 * s * s;
            }

            // Another exceptional shift at 30x. (MATLAB's ad hoc shift.)
            if (iter == 30) {
              s = (y - x) / 2;
              s = s * s + w;
              if (s > 0) {
                s = sqrt(s);
                if (y < x) s = -s;
                s = x - w / ((y - x) / 2.0 + s);
                for (auto ii = 0; ii <= nn; ++ii) matrix->at(ii, ii) = -s;
                shift += s;
                x = y = w = 0.964;
              }
            }
          }

          // Just stop.
          if (iter == 200 * side) 
            throw std::runtime_error("Too many iterations.");

          iter++;
          int m;
          for (m = nn - 2; m >= l; m--) {
            // Look for two consecutive small sub-diagonal elements
            z = matrix->at(m, m);
            r = x - z;
						s = y - z;
            p = (r * s - w) / matrix->at(m + 1, m) + matrix->at(m, m + 1);
            q = matrix->at(m + 1, m + 1) - z - r - s;
            r = matrix->at(m + 2, m + 1);
            // Scale to prevent overﬂow or underflow.
            s = abs(p) + abs(q) + abs(r);
            p /= s;
            q /= s;
            r /= s;

            if (m == l) break;
            double u = abs(matrix->at(m, m - 1)) * (abs(q) + abs(r));
            double v = abs(p) * (abs(matrix->at(m - 1, m - 1)) + abs(z) + abs(matrix->at(m + 1, m + 1)));
            if (u <= epsilon * v) break;
          }
          for (int i = m; i < nn - 1; i++) {
            matrix->at(i + 2, i) = 0;
            if (i != m) matrix->at(i + 2, i - 1) = 0;
          }
          for (int k = m; k < nn; k++) {
            // Double QR step on rows l to nn and columns m to nn.
            if (k != m) {
              p = matrix->at(k, k - 1); // Begin setup of Householder vector
              q = matrix->at(k + 1, k - 1);
              r = 0;
              if (k + 1 != nn) r = matrix->at(k + 2, k - 1);
              if ((x = abs(p) + abs(q) + abs(r)) != 0) {
								p /= x;
								q /= x;
								r /= x;
							}
            }
            if ((s = sign(sqrt(p * p + q * q + r * r), p)) != 0.0) {
              if (k == m) {
                if (l != m) matrix->at(k, k - 1) = -matrix->at(k, k - 1);
              }
              else matrix->at(k, k - 1) = -s * x;
              p += s;
							x = p / s;
							y = q / s;
							z = r / s;
							q /= p;
							r /= p;
              // Row modification.
              for (auto rr = k; rr < side; ++rr) {
                p = matrix->at(k, rr) + q * matrix->at(k + 1, rr);
                if (k + 1 != nn) {
                  p += r * matrix->at(k + 2, rr);
                  matrix->at(k + 2, rr) -= p * z;
                }
                matrix->at(k + 1, rr) -= p * y;
                matrix->at(k, rr) -= p * x;
              }
              // Column modification.
              int minn = fmin(nn, k + 3);
              for (auto cc = 0; cc <= minn; ++cc) {
                p = x * matrix->at(cc, k) + y * matrix->at(cc, k + 1);
                if (k + 1 != nn) {
                  p += z * matrix->at(cc, k + 2);
                  matrix->at(cc, k + 2) -= p * r;
                }
                matrix->at(cc, k + 1) -= p * q;
                matrix->at(cc, k) -= p;
              }
              // Accumulate transformations
              for (auto ii = 0; ii < side; ++ii) {
                p = x * evector->at(ii, k) + y * evector->at(ii, k + 1);
                if (k + 1 != nn) {
                  p += z * evector->at(ii, k + 2);
                  evector->at(ii, k + 2) -= p * r;
                }
                evector->at(ii, k + 1) -= p * q;
                evector->at(ii, k) -= p;
              }
            }
          }
        }
      }
    } while (l + 1 < nn);
  }

  // All roots found. Backsubstitute to find vectors of upper triangular form.
  if (norm) {
    for (nn = side - 1; nn >= 0; nn--) {
      p = real(eigenvalues[nn]);
      q = imag(eigenvalues[nn]);
      int na = nn - 1;
      if (!q) {
        // Real vector.
        int m = nn;
        matrix->at(nn, nn) = 1;
        for (int i = nn - 1; i >= 0; --i) {
          w = matrix->at(i, i) - p;
          r = 0;
          for (int j = m; j <= nn; j++) 
            r += matrix->at(i, j) * matrix->at(j, nn);
          if (imag(eigenvalues[i]) < 0) {
            z = w;
            s = r;
          }
          else {
            m = i;
            if (imag(eigenvalues[i]) == 0) {
              if (!w)
                matrix->at(i, nn) = -r / (norm * epsilon);
              else 
                matrix->at(i, nn) = -r / w;
            }
            else {
              // Solve real equations.
              x = matrix->at(i, i + 1);
              y = matrix->at(i + 1, i);
              q = pow(real(eigenvalues[i]) - p, 2) + pow(imag(eigenvalues[i]), 2);
              double t = (x * s - z * r) / q;
              matrix->at(i, nn) = t;
              matrix->at(i + 1, nn) = (abs(x) > abs(z)) ? (-r - w * t) / x : (-s - y * t) / z;
            }
            // Overflow control.
            double tt = abs(matrix->at(i, nn));
            if (tt * tt * epsilon > 1) {
              for (int j = i; j <= nn; ++j)
                matrix->at(j, nn) /= tt;
            }
          }
        }
      }
      else if (q < 0) {
        int m = na;
        if (abs(matrix->at(nn, na)) > abs(matrix->at(na, nn))) {
          matrix->at(na, na) = q / matrix->at(nn, na);
          matrix->at(na, nn) = -(matrix->at(nn, nn) - p) / matrix->at(nn, na);
        }
        else {
          //complex<double> a(0, -matrix->at(na, nn));
          //complex<double> b(matrix->at(na, na) - p, q);
          //complex<double> temp = a / b;
          complex<double> temp = complex<double>(0, -matrix->at(na, nn)) / complex<double>(matrix->at(na, na) - p, q);
          matrix->at(na, na) = real(temp);
          matrix->at(na, nn) = imag(temp);
        }
        matrix->at(nn, na) = 0;
        matrix->at(nn, nn) = 1;
        for (int i = nn - 2; i >= 0; i--) {
          double ra = 0;
          double sa = 0;
          w = matrix->at(i, i) - p;
          for (int j = m; j < nn; j++) {
            ra += matrix->at(i, j) * matrix->at(j, na);
            sa += matrix->at(i, j) * matrix->at(j, nn);
          }
          if (imag(eigenvalues[i]) < 0) {
            z = w;
						r = ra;
						s = sa;
          }
          else {
            m = i;
            if (imag(eigenvalues[i]) == 0) {
              complex<double> temp = complex(-ra, -sa) / complex(w, q);
              matrix->at(i, na) = real(temp);
              matrix->at(i, nn) = imag(temp);
            }
            else {
              // Solve complex equations.
              x = matrix->at(i, i + 1);
	            y = matrix->at(i + 1, i);
              double vr = pow(real(eigenvalues[i]) - p, 2) + pow(imag(eigenvalues[i]), 2) - pow(q, 2);
              double vi = (real(eigenvalues[i]) - p) * 2 * q;
              if (!vr && !vi)
                vr = epsilon * norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z));

              complex<double> temp = complex(x * r - z * ra + q * sa, x * s - z * sa - q * ra) / complex(vr, vi);
              matrix->at(i, na) = real(temp);
              matrix->at(i, nn) = imag(temp);
              if (abs(x) > (abs(z) + abs(q))) {
                matrix->at(i + 1, na) = (-ra - w * matrix->at(i, na) + q * matrix->at(i, nn)) / x;
                matrix->at(i + 1, nn) = (-sa - w * matrix->at(i, nn) + q * matrix->at(i, na)) / x;
              }
              else {
                complex<double> temp = complex(-r - y * matrix->at(i, na), -s - y * matrix->at(i, nn)) / complex(z, q);
                matrix->at(i + 1, na) = real(temp);
                matrix->at(i + 1, nn) = imag(temp);
              }
            }
          }
          // Overflow control.
          double tt = fmax(abs(matrix->at(i, na)), abs(matrix->at(i, nn)));
          if (tt * tt * epsilon > 1) {
            for (auto j = i; j <= nn; ++j) {
              matrix->at(j, na) /= tt;
              matrix->at(j, nn) /= tt;
            }
          }
        }
      }
    }
    // Back transformation to get eigenvectors of original matrix
    // Multiply by transformation matrix to give vectors of original full matrix.
    for (int j = side - 1; j >= 0; j--) {
      for (int i = 0; i < side; ++i) {
        z = 0;
        for (int k = 0; k <= j; ++k)
          z += evector->at(i, k) * matrix->at(k, j);
        evector->at(i, j) = z;
      }
    }
  }
}
