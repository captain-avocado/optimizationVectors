//
// Created by Артем Калоев on 05/12/2018.
//

#ifndef MO_MATRIX_H
#define MO_MATRIX_H


#include <iostream>
#include <vector>
#include <math.h>

using std::cout;
using std::cin;
using std::vector;

class Matrix {

private:
    bool isVector;

//    vector< vector<double> > values;
    void createMatrix();
public:
    int rows, cols;
    
    vector< vector<double> > values;

    Matrix();
    Matrix(int);
    Matrix(vector<double>);
    Matrix(int, int);
    Matrix(Matrix &matrixToCopy);
//    ~Matrix(); //не нужен

    int getRows() const;
    int getCols() const;
    void setMatrix();
    void printMatrix();
    double getByIndex(int, int) const;
    double setByIndex(int, int, double);
    double norma();

    Matrix transpose();
    vector<double> getVector(int);

    Matrix operator+ (const Matrix &right);
    Matrix operator- (const Matrix &right);
    Matrix divide(double divider);

    //оператор умножения на скаляр –– скаляр всегда справа
    Matrix operator* (const double &alpha);
    Matrix operator* (const Matrix &right);
    Matrix operator/ (const double &right);
    const Matrix &operator= (const Matrix &right);
    bool operator== (const Matrix &right) const;

};


#endif //MO_MATRIX_H
