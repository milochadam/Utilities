#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>

class Matrix {
   public:
    Matrix(int, int, int);
    Matrix(int, int);
    Matrix();
    ~Matrix();
    Matrix(const Matrix&);
    Matrix& operator=(const Matrix&);

    double& operator()(int x, int y) { return p[y][x]; }

    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator*=(const Matrix&);
    Matrix& operator*=(double);
    Matrix& operator/=(double);

    friend std::ostream& operator<<(std::ostream&, const Matrix&);

   private:
    int rows_, cols_;
    double** p;

    void allocSpace();

   public:
    int x() const { return cols_; }
    int y() const { return rows_; }
    void initRibbon(double, double, double);
    void initSin(int);
};

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Matrix operator/(const Matrix&, double);

int solveJacobi(Matrix& A, Matrix& x, Matrix& b);
int solveGaussSeidel(Matrix& A, Matrix& x, Matrix& b);
int solveLU(Matrix& A, Matrix& x, Matrix& b);

double calculate_residuum(Matrix& A, Matrix& x, Matrix& b);
#endif
