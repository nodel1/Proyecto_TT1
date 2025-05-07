// $Header$
//--------------------------------------------------
// Matrix
//--------------------------------------------------
// Proyecto_TT1: Matrix Operations header
// 
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/11
//
/** @file matrix.hpp
 *  @brief Cabecera principal para operaciones de matrix
 *  
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

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
	
	//constructor vacio
	Matrix();

    // Parameterized constructor
    Matrix(const int n_row, const int n_column);
    Matrix(const int n);
    
    // Member operators
	double& operator () (const int row, const int column);
	double operator()(const int row, const int column) const;
    double& operator () (const int n);
    Matrix& operator + (Matrix &m);
    Matrix& operator - (Matrix &m);
    Matrix& operator * (Matrix &m);
    Matrix& operator / (Matrix &m);
    Matrix& operator = (const Matrix &m);
    Matrix& operator * (double scalar);
    Matrix& operator / (double scalar);
    Matrix& operator + (double scalar);
    Matrix& operator - (double scalar);
    
    // Non-member operators
    friend ostream& operator << (ostream &o, Matrix &m);
};

ostream& operator << (ostream &o, Matrix &m);

//--------------------------------------------------
// inv (Matrix &m)
//--------------------------------------------------
/**
 * @brief Calculates the inverse of a matrix.
 *
 * @param [in,out] m The input matrix to be inverted.
 * @return Reference to the inverted matrix.
 */
Matrix& inv(Matrix &m);

//--------------------------------------------------
// transpose (Matrix &m)
//--------------------------------------------------
/**
 * @brief Calculates the transpose of a matrix.
 *
 * @param [in,out] m The input matrix to be transposed.
 * @return Reference to the transposed matrix.
 */
Matrix& transpose(Matrix &m);

//--------------------------------------------------
// norm (Matrix &m)
//--------------------------------------------------
/**
 * @brief Calculates the norm of a matrix.
 *
 * @param [in] m The input matrix.
 * @return The norm of the matrix as a double value.
 */
double norm(Matrix &m);

//--------------------------------------------------
// dot (Matrix &a, Matrix &b)
//--------------------------------------------------
/**
 * @brief Calculates the dot product of two matrices.
 *
 * @param [in] a The first input matrix.
 * @param [in] b The second input matrix.
 * @return The dot product as a double value.
 */
double dot(Matrix &a, Matrix &b);

//--------------------------------------------------
// cross (Matrix &a, Matrix &b)
//--------------------------------------------------
/**
 * @brief Calculates the cross product of two matrices.
 *
 * @param [in] a The first input matrix.
 * @param [in] b The second input matrix.
 * @return Reference to the resulting matrix from the cross product.
 */
Matrix& cross(Matrix &a, Matrix &b);

//--------------------------------------------------
// extract_vector (Matrix &m)
//--------------------------------------------------
/**
 * @brief Extracts a vector from a matrix.
 *
 * @param [in] m The input matrix.
 * @return Reference to the extracted vector as a matrix.
 */
Matrix& extract_vector(Matrix &m);

//--------------------------------------------------
// union_vector (Matrix &v, int rows, int cols)
//--------------------------------------------------
/**
 * @brief Constructs a matrix by reshaping a vector.
 *
 * @param [in] v The input vector.
 * @param [in] rows The number of rows for the resulting matrix.
 * @param [in] cols The number of columns for the resulting matrix.
 * @return Reference to the resulting matrix.
 */
Matrix& union_vector(Matrix &v, int rows, int cols);

//--------------------------------------------------
// extract_row (Matrix &m, int row)
//--------------------------------------------------
/**
 * @brief Extracts a specific row from a matrix.
 *
 * @param [in] m The input matrix.
 * @param [in] row The index of the row to extract.
 * @return Reference to the extracted row as a matrix.
 */
Matrix& extract_row(Matrix &m, int row);

//--------------------------------------------------
// extract_column (Matrix &m, int col)
//--------------------------------------------------
/**
 * @brief Extracts a specific column from a matrix.
 *
 * @param [in] m The input matrix.
 * @param [in] col The index of the column to extract.
 * @return Reference to the extracted column as a matrix.
 */
Matrix& extract_column(Matrix &m, int col);

//--------------------------------------------------
// assign_row (Matrix &m, int row, Matrix &v)
//--------------------------------------------------
/**
 * @brief Assigns a vector to a specific row in a matrix.
 *
 * @param [in,out] m The input matrix to modify.
 * @param [in] row The index of the row to assign.
 * @param [in] v The vector to assign to the row.
 */
void assign_row(Matrix &m, int row, Matrix &v);

//--------------------------------------------------
// assign_column (Matrix &m, int col, Matrix &v)
//--------------------------------------------------
/**
 * @brief Assigns a vector to a specific column in a matrix.
 *
 * @param [in,out] m The input matrix to modify.
 * @param [in] col The index of the column to assign.
 * @param [in] v The vector to assign to the column.
 */
void assign_column(Matrix &m, int col, Matrix &v);



// Special matrices
Matrix& zeros(const int n_row, const int n_column);
Matrix& zeros(const int n);
Matrix& eye(const int n);

#endif