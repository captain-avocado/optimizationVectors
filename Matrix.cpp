//
// Created by Артем Калоев on 05/12/2018.
//

#include "Matrix.h"

using std::cout;
using std::cin;
using std::endl;

void Matrix::createMatrix() {
    values.resize(cols);
    for (int i = 0; i < cols; i++) {
        values[i].resize(rows);
        for (int j = 0; j < rows; j++) {
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

Matrix::Matrix(vector<double> x) {
    rows = x.size();
    cols = 1;
    values.resize(cols);
    for (int i = 0; i < cols; i++) {
        values[i].resize(rows);
        for (int j = 0; j < rows; j++) {
            values[i][j] = x[i];
        }
    }
}

Matrix::Matrix(Matrix &matrixToCopy) {
    rows = matrixToCopy.rows;
    cols = matrixToCopy.cols;
    values.resize(cols);
    for (int i = 0; i < cols; i++) {
        values[i].resize(rows);
        for (int j = 0; j < rows; j++) {
            values[i][j] = matrixToCopy.values[i][j];
        }
    }
}

int Matrix::getRows() const {   return rows;    }
int Matrix::getCols() const {   return cols;    }

void Matrix::printMatrix() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << values[j][i] << " ";
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
            cin >> values[j][i];
        }
    }
    printMatrix();
}

double Matrix::getByIndex(int i, int j) const {
    if (i < rows && j < cols) {
        return values[j][i];
    } else {
        cout << "Ошибка: индексы превышают границы матрицы" << endl;
        exit(1);
    }
}

double Matrix::setByIndex(int i, int j, double val) {
    if (i < rows && j < cols) {
        values[j][i] = val;
    } else {
        cout << "Ошибка: индексы превышают границы матрицы" << endl;
        exit(1);
    }
}

Matrix Matrix::transpose() {
//    int tmpCols = cols;
//    cols = roes
}

vector<double> Matrix::getVector(int i) {
    return values[i];
}

Matrix Matrix::operator+ (const Matrix &right) {
    if (rows == right.rows && cols == right.cols) {
        Matrix result(rows, cols);
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
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
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                result.values[i][j] = values[i][j] - right.values[i][j];
            }
        }
        return result;
    } else {
        cout << "Матрицы разного размера. Вычитание невозможно" << endl;
    }
}

Matrix Matrix::operator*(const double &alpha) {
    Matrix result(rows, cols);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            result.values[i][j] = values[i][j] * alpha;
        }
    }
    return result;
}

const Matrix& Matrix::operator=(const Matrix &right) {
    if (&right != this) {
        rows = right.rows;
        cols = right.cols;
        values.resize(cols);
        for (int i = 0; i < cols; i++) {
            values[i].resize(rows);
            for (int j = 0; j < rows; j++) {
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
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            if (values[i][j] != right.values[i][j]) {
                return false;
            }
        }
    }

    return true;
}



