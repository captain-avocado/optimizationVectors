//
// Created by Артем Калоев on 05/12/2018.
//

#include "Matrix.h"

using std::cout;
using std::cin;
using std::endl;

void Matrix::createMatrix() {
    values.resize(rows);
    for (int i = 0; i < rows; i++) {
        values[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            values[i][j] = 0;
        }
    }
}

Matrix::Matrix() {
    rows = cols = 5;
    createMatrix();
}

Matrix::Matrix(int newRows, int newCols) {
    rows = newRows;
    cols = newCols;
    createMatrix();
}

Matrix::Matrix(int newRows) {
    rows = newRows;
    cols = 1;
    createMatrix();
}

Matrix::Matrix(Matrix &matrixToCopy) {
    rows = matrixToCopy.rows;
    cols = matrixToCopy.cols;
    values.resize(rows);
    for (int i = 0; i < rows; i++) {
        values[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            values[i][j] = matrixToCopy.values[i][j];
        }
    }
}

int Matrix::getRows() const {   return rows;    }
int Matrix::getCols() const {   return cols;    }

void Matrix::printMatrix() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << values[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Matrix::setMatrix() {
    cout << "Cols = " << cols << ", rows = " << rows << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << "[" << i << "][" << j << "] = ";
            cin >> values[i][j];
        }
    }
    printMatrix();
}

double Matrix::getByIndex(int i, int j) const {
    return values[i][j];
}

double Matrix::setByIndex(int i, int j, double val) {
    values[i][j] = val;
}

Matrix Matrix::operator+ (const Matrix &right) {
    if (rows == right.rows && cols == right.cols) {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.values[i][j] = values[i][j] + right.values[i][j];
            }
        }
        return result;
    } else {
        cout << "Матрицы разного размера. Сложение невозможно" << endl;
    }
}

Matrix Matrix::operator- (const Matrix &right) {
    if (rows == right.rows && cols == right.cols) {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.values[i][j] = values[i][j] - right.values[i][j];
            }
        }
        return result;
    } else {
        cout << "Матрицы разного размера. Вычитание невозможно" << endl;
    }
}

const Matrix& Matrix::operator=(const Matrix &right) {
    if (&right != this) {
        rows = right.rows;
        cols = right.cols;
        values.resize(rows);
        for (int i = 0; i < rows; i++) {
            values[i].resize(cols);
            for (int j = 0; j < cols; j++) {
                values[i][j] = right.values[i][j];
            }
        }
    }
    return *this;
}

bool Matrix::operator== (const Matrix &right) const {
    if (rows != right.rows || cols != right.cols) {
        return false;
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (values[i][j] != right.values[i][j]) {
                return false;
            }
        }
    }

    return true;
}



