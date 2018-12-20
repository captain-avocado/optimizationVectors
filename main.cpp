#include <iostream>
#include "Matrix.h"
#include "Optimize.h"
#include "Vector.h"

using std::cout;
using std::cin;
using std::endl;

int main() {
    std::string expression = "(x[0] - 1)^2 + (x[1] - 3)^2 + 4*(x[2] + 5)^2";
    Vector x0(3);
    x0.values[0] = 4;
    x0.values[1] = -1;
    x0.values[2] = 2;

//    x0.setVector();
    Optimize op(expression, x0);

    op.Koshi();


//    Matrix x;
//    x = p.transpose();
//    x.printMatrix();
//    x = x0.transpose() * p;
//    x.printMatrix();
//    op.CGM();
//    op.MPK();
//    for (int i = 0; i < 7; i++) {
//        for (int j = 0; j < 5; j++) {
//            cout << "Связка: " << op.getOneDimMethodByIndex(i) << " + " << op.getInterpMethodByIndex(j) << endl;
//            op.setIndexQueue(i, j);
//            op.CGM();
//            x = op.newPoint(op.alphaSearch());
//            cout << "Результат: " <<   x.values[0][0] << " " << x.values[0][1]   << endl;
//        }
//
//    }

//    Vector h(2);
//    h.setVector();
//    Vector m(2);
//    m.setVector();
//    cout << m * h << endl;


    return 0;
}