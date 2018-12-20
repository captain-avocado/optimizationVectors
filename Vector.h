//
// Created by Артем Калоев on 19/12/2018.
//

#ifndef MO_VECTOR_H
#define MO_VECTOR_H

#include <iostream>
#include <vector>
#include <math.h>

using std::cout;
using std::cin;
using std::vector;

class Vector {
public:
    int n;
    vector<double> values;

    Vector();
    Vector(int);
    Vector(Vector &vectorToCopy);
    Vector(vector<double>);

    void setVector();
    void printVector();

    double norma();

    Vector operator+ (const Vector &right);
    Vector operator- (const Vector &right);

    //оператор умножения на скаляр –– скаляр всегда справа
    Vector operator* (const double &alpha);
    double operator* (const Vector &right);
    const Vector &operator= (const Vector &right);
    bool operator== (const Vector &right) const;
};


#endif //MO_VECTOR_H
