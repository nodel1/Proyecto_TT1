#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
    double **data;

    // Parameterized constructor
    Matrix(const int n_row, const int n_column);
    Matrix(const int n);
    
    // Member operators
	double& operator () (const int row, const int column);
    double& operator () (const int n);
    Matrix& operator + (Matrix &m);
    Matrix& operator - (Matrix &m);
    Matrix& operator * (Matrix &m);
    Matrix& operator / (Matrix &m);
    Matrix& operator = ( Matrix &m);
    Matrix& operator * (double scalar);
    Matrix& operator / (double scalar);
    Matrix& operator + (double scalar);
    Matrix& operator - (double scalar);
    
    // Non-member operators
    friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Matrix operations
Matrix& inv(Matrix &m);
Matrix& transpose(Matrix &m);
double norm(Matrix &m);
double dot(Matrix &a, Matrix &b);
Matrix& cross(Matrix &a, Matrix &b);
Matrix& extract_vector(Matrix &m);
Matrix& union_vector(Matrix &v, int rows, int cols);
Matrix& extract_row(Matrix &m, int row);
Matrix& extract_column(Matrix &m, int col);
void assign_row(Matrix &m, int row, Matrix &v);
void assign_column(Matrix &m, int col, Matrix &v);

// Special matrices
Matrix& zeros(const int n_row, const int n_column);
Matrix& zeros(const int n);
Matrix& eye(const int n);

#endif