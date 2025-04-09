#include "..\include\matrix.h"
#include <iostream>
#include <cmath>

using namespace std;

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
        cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }
    
    this->n_row = n_row;
    this->n_column = n_column;
    this->data = (double **) malloc(n_row * sizeof(double *));
    
    if (this->data == NULL) {
        cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < n_row; i++) {
        this->data[i] = (double *) malloc(n_column * sizeof(double));
    }
}


double& Matrix::operator () (const int row, const int column) {
    if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
        cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
    }
    
    return this->data[row - 1][column - 1];
}

Matrix& Matrix::operator + (Matrix &m) {
    if (this->n_row != m.n_row || this->n_column != m.n_column) {
        cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) + m(i, j);
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
    if (this->n_row != m.n_row || this->n_column != m.n_column) {
        cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
    }
    
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);
    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) - m(i, j);
        }
    }
    
    return *m_aux;
}

Matrix& Matrix::operator * (Matrix &m) {
    if (this->n_column != m.n_row) {
        cout << "Matrix multiplication: dimension mismatch\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(this->n_row, m.n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            double sum = 0.0;
            for (int k = 1; k <= this->n_column; k++) {
                sum += (*this)(i, k) * m(k, j);
            }
            (*m_aux)(i, j) = sum;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator / (Matrix &m) {
    if (m.n_row != m.n_column || this->n_column != m.n_row) {
        cout << "Matrix division: dimensions must be square and compatible\n";
        exit(EXIT_FAILURE);
    }

    Matrix inv_m = inv(m);
    Matrix *m_aux = new Matrix(this->n_row, inv_m.n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= inv_m.n_column; j++) {
            double sum = 0.0;
            for (int k = 1; k <= this->n_column; k++) {
                sum += (*this)(i, k) * inv_m(k, j);
            }
            (*m_aux)(i, j) = sum;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator = (Matrix &m) {
    if (this == &m)
        return *this;

    if (this->n_row != m.n_row || this->n_column != m.n_column) {
        this->n_row = m.n_row;
        this->n_column = m.n_column;

        this->data = new double*[n_row];
        for (int i = 0; i < n_row; ++i)
            this->data[i] = new double[n_column];
    }

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*this)(i, j) = m(i, j);
        }
    }

    return *this;
}

ostream& operator << (ostream &o, Matrix &m) {
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
            printf("%5.20lf ", m(i,j));
        o << "\n";
    }
    
    return o;
}

Matrix& zeros(const int n_row, const int n_column) {
    Matrix *m_aux = new Matrix(n_row, n_column);
    
    for (int i = 1; i <= n_row; i++) {
        for (int j = 1; j <= n_column; j++) {
            (*m_aux)(i, j) = 0;
        }
    }
    
    return *m_aux;
}

Matrix& eye(const int n) {
    Matrix *result = new Matrix(n, n);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*result)(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    return *result;
}

Matrix& transpose(Matrix &m) {
    Matrix *result = new Matrix(m.n_column, m.n_row);

    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*result)(j, i) = m(i, j);
        }
    }

    return *result;
}

Matrix& inv(Matrix &m) {
    if (m.n_row != m.n_column) {
        cout << "Matrix inverse: matrix must be square\n";
        exit(EXIT_FAILURE);
    }

    int n = m.n_row;
    Matrix *result = new Matrix(n, n);
    Matrix temp(n, n);
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            temp(i, j) = m(i, j);
            (*result)(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int k = 1; k <= n; k++) {
        int max_row = k;
        for (int i = k + 1; i <= n; i++) {
            if (fabs(temp(i, k)) > fabs(temp(max_row, k))) {
                max_row = i;
            }
        }

        if (max_row != k) {
            for (int j = 1; j <= n; j++) {
                double temp_val = temp(k, j);
                temp(k, j) = temp(max_row, j);
                temp(max_row, j) = temp_val;

                temp_val = (*result)(k, j);
                (*result)(k, j) = (*result)(max_row, j);
                (*result)(max_row, j) = temp_val;
            }
        }

        double pivot = temp(k, k);
        if (fabs(pivot) < 1e-10) {
            cout << "Matrix inverse: singular matrix (no inverse exists)\n";
            exit(EXIT_FAILURE);
        }

        for (int j = 1; j <= n; j++) {
            temp(k, j) /= pivot;
            (*result)(k, j) /= pivot;
        }

        for (int i = 1; i <= n; i++) {
            if (i != k) {
                double factor = temp(i, k);
                for (int j = 1; j <= n; j++) {
                    temp(i, j) -= factor * temp(k, j);
                    (*result)(i, j) -= factor * (*result)(k, j);
                }
            }
        }
    }

    return *result;
}

Matrix& Matrix::operator + (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) + scalar;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator - (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) - scalar;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator * (double scalar) {
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) * scalar;
        }
    }

    return *m_aux;
}

Matrix& Matrix::operator / (double scalar) {
    if (scalar == 0.0) {
        cout << "Matrix division: division by zero\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(i, j) = (*this)(i, j) / scalar;
        }
    }

    return *m_aux;
}

Matrix::Matrix(const int n) {
    if (n <= 0) {
        cout << "Matrix create: error in n\n";
        exit(EXIT_FAILURE);
    }
    
    this->n_row = n;
    this->n_column = n;
    this->data = (double **) malloc(n * sizeof(double *));
    
    if (this->data == NULL) {
        cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < n; i++) {
        this->data[i] = (double *) malloc(n * sizeof(double));
    }
}

double& Matrix::operator () (const int n) {
    if (this->n_row != 1 && this->n_column != 1) {
        cout << "Matrix vector access: matrix must be a row or column vector\n";
        exit(EXIT_FAILURE);
    }
    if (n <= 0 || n > (this->n_row * this->n_column)) {
        cout << "Matrix vector access: index out of bounds\n";
        exit(EXIT_FAILURE);
    }

    if (this->n_row == 1) {
        return this->data[0][n - 1];
    } else {
        return this->data[n - 1][0];
    }
}

Matrix& zeros(const int n) {
    Matrix *result = new Matrix(n, n);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*result)(i, j) = 0.0;
        }
    }

    return *result;
}

double norm(Matrix &m) {
    double sum = 0.0;
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            sum += m(i, j) * m(i, j);
        }
    }
    return sqrt(sum);
}

double dot(Matrix &a, Matrix &b) {
    if ((a.n_row != 1 && a.n_column != 1) || (b.n_row != 1 && b.n_column != 1) ||
        (a.n_row * a.n_column != b.n_row * b.n_column)) {
        cout << "Dot product: matrices must be vectors of same length\n";
        exit(EXIT_FAILURE);
    }

    double result = 0.0;
    int length = a.n_row * a.n_column;
    for (int i = 1; i <= length; i++) {
        result += a(i) * b(i);
    }
    return result;
}

Matrix& cross(Matrix &a, Matrix &b) {
    if ((a.n_row != 3 || a.n_column != 1) && (a.n_row != 1 || a.n_column != 3) ||
        (b.n_row != 3 || b.n_column != 1) && (b.n_row != 1 || b.n_column != 3)) {
        cout << "Cross product: matrices must be 3x1 or 1x3 vectors\n";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(3, 1);
    double a1 = a(1), a2 = a(2), a3 = a(3);
    double b1 = b(1), b2 = b(2), b3 = b(3);

    (*result)(1, 1) = a2 * b3 - a3 * b2;
    (*result)(2, 1) = a3 * b1 - a1 * b3;
    (*result)(3, 1) = a1 * b2 - a2 * b1;

    return *result;
}

Matrix& extract_vector(Matrix &m) {
    Matrix *result = new Matrix(1, m.n_row * m.n_column);
    int k = 1;
    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*result)(1, k++) = m(i, j);
        }
    }
    return *result;
}

Matrix& union_vector(Matrix &v, int rows, int cols) {
    if (v.n_row * v.n_column != rows * cols) {
        cout << "Union vector: vector length must match matrix size\n";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(rows, cols);
    int k = 1;
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            (*result)(i, j) = v(k++);
        }
    }
    return *result;
}

Matrix& extract_row(Matrix &m, int row) {
    if (row <= 0 || row > m.n_row) {
        cout << "Extract row: invalid row index\n";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(1, m.n_column);
    for (int j = 1; j <= m.n_column; j++) {
        (*result)(1, j) = m(row, j);
    }
    return *result;
}

Matrix& extract_column(Matrix &m, int col) {
    if (col <= 0 || col > m.n_column) {
        cout << "Extract column: invalid column index\n";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(m.n_row, 1);
    for (int i = 1; i <= m.n_row; i++) {
        (*result)(i, 1) = m(i, col);
    }
    return *result;
}

void assign_row(Matrix &m, int row, Matrix &v) {
    if (row <= 0 || row > m.n_row || v.n_row != 1 || v.n_column != m.n_column) {
        cout << "Assign row: invalid row index or vector size\n";
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j <= m.n_column; j++) {
        m(row, j) = v(1, j);
    }
}

void assign_column(Matrix &m, int col, Matrix &v) {
    if (col <= 0 || col > m.n_column || v.n_column != 1 || v.n_row != m.n_row) {
        cout << "Assign column: invalid column index or vector size\n";
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= m.n_row; i++) {
        m(i, col) = v(i, 1);
    }
}