//
// Created by Артем Калоев on 19/12/2018.
//

#include "Vector.h"

using std::endl;
using std::cin;
using std::cout;

Vector::Vector() {
    n = 5;
    values.resize(5);
    for (int i = 0; i < n; i++) {
        values[i] = 0;
    }
}

Vector::Vector(int size) {
    n = size;
    values.resize(size);
    for (int i = 0; i < n; i++) {
        values[i] = 0;
    }
}

Vector::Vector(vector<double> x) {
    n = x.size();
    values.resize(n);
    for (int i = 0; i < n; i++) {
        values[i] = x[i];
    }
}

Vector::Vector(Vector &x) {
    n = x.n;
    values.resize(n);
    for (int i = 0; i < n; i++) {
        values[i] = x.values[i];
    }
}

void Vector::printVector() {
    for (int i = 0; i < n; i++) {
        cout << values[i] << " ";
    }
    cout << endl;
}

void Vector::setVector() {
    cout << "Size = " << n << endl;
    for (int i = 0; i < n; i++) {
        cin >> values[i];
    }
    printVector();
}

double Vector::norma() {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += pow(values[i], 2);
    }
    return sqrt(sum);
}

Vector Vector::operator+ (const Vector &right) {
    if (n == right.n) {
        Vector result(n);
        for (int i = 0; i < n; i++) {
            result.values[i] = values[i] + right.values[i];
        }
        return result;
    } else {
        cout << "Векторы разного размера. Сложение невозможно" << endl;
    }
}

Vector Vector::operator- (const Vector &right) {
    if (n == right.n) {
        Vector result(n);
        for (int i = 0; i < n; i++) {
            result.values[i] = values[i] - right.values[i];
        }
        return result;
    } else {
        cout << "Векторы разного размера. Вычитание невозможно" << endl;
    }
}

Vector Vector::operator*(const double &alpha) {
    Vector result(n);
    for (int i = 0; i < n; i++) {
        result.values[i] = values[i] * alpha;
    }
    return result;
}

double Vector::operator*(const Vector &right) {
    //транспонированное умножение
    double res = 0;
    if (n == right.n) {
        for (int i = 0; i < n; i++) {
            res += values[i] * right.values[i];
        }
    } else {
        cout << "Перемножение векторов невозможно. Векторы разных размеров" << endl;
    }
    return res;
}


const Vector& Vector::operator=(const Vector &right) {
    if (&right != this) {
        n = right.n;
        values.resize(n);
        for (int i = 0; i < n; i++) {
            values[i] = right.values[i];
        }
    }
    return *this;
}

bool Vector::operator== (const Vector &right) const {
    if (n != right.n) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        if (values[i] != right.values[i]) {
            return false;
        }
    }
    return true;
}