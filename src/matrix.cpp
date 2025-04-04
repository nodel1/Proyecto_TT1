#include "..\include\matrix.h"

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
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
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
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
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
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
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& Matrix::operator * (Matrix &m) {  // Multiplicación de matrices
    if (this->n_column != m.n_row) {
        cout << "Matrix multiplication: dimension mismatch\n";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(this->n_row, m.n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            double sum = 0.0;
            for (int k = 1; k <= this->n_column; k++) {
                sum += (*this)(i, k) * m(k, j);
            }
            (*result)(i, j) = sum;
        }
    }

    return *result;
}


Matrix& inv(Matrix &m) {                    //DE INTERNET
    if (m.n_row != m.n_column) {
        cout << "Matrix inverse: matrix must be square\n";
        exit(EXIT_FAILURE);
    }

    int n = m.n_row;
    Matrix *result = new Matrix(n, n); // Matriz inversa a devolver
    Matrix temp = m;                  // Copia de la matriz original

    // Inicializar 'result' como matriz identidad
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*result)(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    // Eliminación Gauss-Jordan
    for (int k = 1; k <= n; k++) {
        // Pivoteo parcial (para evitar división por cero)
        int max_row = k;
        for (int i = k + 1; i <= n; i++) {
            if (fabs(temp(i, k)) > fabs(temp(max_row, k))) {
                max_row = i;
            }
        }

        // Intercambiar filas si es necesario
        if (max_row != k) {
            for (int j = 1; j <= n; j++) {
                std::swap(temp(k, j), temp(max_row, j));
                std::swap((*result)(k, j), (*result)(max_row, j));
            }
        }

        double pivot = temp(k, k);
        if (fabs(pivot) < 1e-10) { // Verificar si la matriz es singular
            cout << "Matrix inverse: singular matrix (no inverse exists)\n";
            exit(EXIT_FAILURE);
        }

        // Normalizar la fila del pivote
        for (int j = 1; j <= n; j++) {
            temp(k, j) /= pivot;
            (*result)(k, j) /= pivot;
        }

        // Eliminación hacia adelante y atrás
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

Matrix& Matrix::operator / (Matrix &m) {     //División de matrices (A/B = A * inv(B))
    if (m.n_row != m.n_column || this->n_column != m.n_row) {
        cout << "Matrix division: dimensions must be square and compatible\n";
        exit(EXIT_FAILURE);
    }

    Matrix inv_m = inv(m); 
    return (*this) * inv_m; 
}

Matrix& Matrix::operator * (double scalar) {    // Multiplicación por escalar
    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*result)(i, j) = (*this)(i, j) * scalar;
        }
    }

    return *result;
}

Matrix& Matrix::operator / (double scalar) {     //División por escalar
    if (scalar == 0.0) {
        cout << "Matrix division: division by zero\n";
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*result)(i, j) = (*this)(i, j) / scalar;
        }
    }

    return *result;
}

Matrix& Matrix::operator + (double scalar) { //Suma de escalar a matriz
    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*result)(i, j) = (*this)(i, j) + scalar;
        }
    }

    return *result;
}

Matrix& Matrix::operator - (double scalar) {  //Resta de escalar a matriz
    Matrix *result = new Matrix(this->n_row, this->n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*result)(i, j) = (*this)(i, j) - scalar;
        }
    }

    return *result;
}

Matrix& eye(const int n) {  //Matriz identidad de tamaño n
    Matrix *result = new Matrix(n, n);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*result)(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    return *result;
}


Matrix& zeros(const int n) {  //Matriz cuadrada de ceros de tamaño n
    Matrix *result = new Matrix(n, n);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*result)(i, j) = 0.0;
        }
    }

    return *result;
}

Matrix& transpose(Matrix &m) {   //Matriz transpuesta
    Matrix *result = new Matrix(m.n_column, m.n_row);

    for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*result)(j, i) = m(i, j);
        }
    }

    return *result;
}







