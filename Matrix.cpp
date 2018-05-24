#include "Matrix.h"
#include <cmath>
#include <iomanip>

#define EPS 1e-10

/* PUBLIC MEMBER FUNCTIONS
 ********************************/

Matrix::Matrix(int cols, int rows, int value) : rows_(rows), cols_(cols) {
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = value;
        }
    }
}

Matrix::Matrix(int cols, int rows) : Matrix(cols, rows, 0) {}

Matrix::Matrix() : rows_(1), cols_(1) {
    allocSpace();
    p[0][0] = 0;
}

Matrix::~Matrix() {
    for (int i = 0; i < rows_; ++i) {
        delete[] p[i];
    }
    delete[] p;
}

Matrix::Matrix(const Matrix& m) : rows_(m.rows_), cols_(m.cols_) {
    allocSpace();
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
}

Matrix& Matrix::operator=(const Matrix& m) {
    if (this == &m) {
        return *this;
    }

    if (rows_ != m.rows_ || cols_ != m.cols_) {
        for (int i = 0; i < rows_; ++i) {
            delete[] p[i];
        }
        delete[] p;

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocSpace();
    }

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& m) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] += m.p[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& m) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] -= m.p[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(const Matrix& m) {
    Matrix temp(rows_, m.cols_);
    for (int i = 0; i < temp.rows_; ++i) {
        for (int j = 0; j < temp.cols_; ++j) {
            for (int k = 0; k < cols_; ++k) {
                temp.p[i][j] += (p[i][k] * m.p[k][j]);
            }
        }
    }
    return (*this = temp);
}

Matrix& Matrix::operator*=(double num) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] *= num;
        }
    }
    return *this;
}

Matrix& Matrix::operator/=(double num) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            p[i][j] /= num;
        }
    }
    return *this;
}

/* HELPER FUNCTIONS */
void Matrix::allocSpace() {
    p = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
        p[i] = new double[cols_];
    }
}

void Matrix::initRibbon(double a1, double a2, double a3) {
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            switch (i - j) {
                case -2:
                case 2:
                    p[i][j] = a3;
                    break;
                case -1:
                case 1:
                    p[i][j] = a2;
                    break;
                case 0:
                    p[i][j] = a1;
                    break;
                default:
                    break;
            }
        }
    }
}

void Matrix::initSin(int f) {
    for (int i = 0; i < rows_; ++i) {
        p[i][0] = sin((i + 1) * (f + 1));
    }
}
/* NON-MEMBER FUNCTIONS */

Matrix operator+(const Matrix& m1, const Matrix& m2) {
    Matrix temp(m1);
    return (temp += m2);
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {
    Matrix temp(m1);
    return (temp -= m2);
}

Matrix operator*(const Matrix& m1, const Matrix& m2) {
    Matrix temp(m1);
    return (temp *= m2);
}

Matrix operator*(const Matrix& m, double num) {
    Matrix temp(m);
    return (temp *= num);
}

Matrix operator*(double num, const Matrix& m) { return (m * num); }

Matrix operator/(const Matrix& m, double num) {
    Matrix temp(m);
    return (temp /= num);
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    for (int i = 0; i < m.rows_; ++i) {
        os << std::setw(3) << m.p[i][0];
        for (int j = 1; j < m.cols_; ++j) {
            os << " " << std::setw(3) << m.p[i][j];
        }
        os << std::endl;
    }
    return os;
}

double calculate_residuum(Matrix& A, Matrix& x, Matrix& b) {
    Matrix res = A * x - b;
    double norm = 0;
    for (int i = 0; i < res.y(); i++) {
        norm += res(0, i) * res(0, i);
    }
    return sqrt(norm);
}

int solveJacobi(Matrix& A, Matrix& x, Matrix& b) {
    int iter = 0;
    double res = 1;
    Matrix x1(1, x.y());
    while (res > EPS) {
        for (int i = 0; i < x.y(); i++) {
            double sum = 0;
            for (int j = 0; j < x.y(); j++) {
                if (i == j) continue;
                sum += A(j, i) * x(0, j);
            }
            x1(0, i) = (b(0, i) - sum) / A(i, i);
        }
        x = x1;
        res = calculate_residuum(A, x, b);
        iter++;
    }
    return iter;
}

int solveGaussSeidel(Matrix& A, Matrix& x, Matrix& b) {
    int iter = 0;
    double res = 1;
    Matrix x1(1, x.y());
    while (res > EPS) {
        for (int i = 0; i < x.y(); i++) {
            double sum = 0;
            for (int j = 0; j < x.y(); j++) {
                if (i == j) continue;
                sum += A(j, i) * x(0, j);
            }
            x(0, i) = (b(0, i) - sum) / A(i, i);
        }
        res = calculate_residuum(A, x, b);
        iter++;
    }
    return iter;
}

int solveLU(Matrix& A, Matrix& x, Matrix& b) {
    int iter = 0;
    Matrix U = A;
    Matrix L = A;
    L.initRibbon(1, 0, 0);
    for (int k = 0; k < A.y() - 1; ++k) {
        for (int j = k + 1; j < A.y(); ++j) {
            L(j, k) = U(j, k) / U(k, k);
            for (int i = k; i < A.y(); ++i) {
                U(j, i) -= L(j, k) * U(k, i);
            }
        }
    }
    Matrix y(1, x.y());
    // L*y=b
    y(0, 0) = b(0, 0);
    for (int i = 1; i < y.y(); ++i) {
        double sum = 0;
        for (int j = 0; j < i - 1; ++j) {
            sum += L(i, j) * y(0, j);
        }
        y(0, i) = b(0, i) - sum;
    }

    // Ux=y;
    x(0, x.y() - 1) = y(0, y.y() - 1) / U(U.x() - 1, U.y() - 1);

    for (int i = x.y() - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < x.y(); ++j) {
            sum += U(i, j) * x(0, j);
        }
        x(0, i) = (y(0, i) - sum) / U(i, i);
    }
    return iter;  // 11626
}
