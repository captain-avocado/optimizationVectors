#include <iostream>
#include "Matrix.h"
#include "Optimize.h"

using std::cout;
using std::cin;
using std::endl;

int main() {
    std::string expression;
//    std::cout << "Введите строку:" << std::endl;
//    cin >> expression;

    expression = "100*(x[1] - x[0]^2)^2 + (1 - x[0])^2";
    Optimize op(expression);

    Matrix p(2), x0(2);
    x0.values[0][0] = 10; x0.values[0][1] = 10;
    p.values[0][0] = 2; //p.values[1][0] = 1;
    p.values[0][1] = 3; //p.values[1][1] = 4;
    op.setX0(x0);
    op.setP(p);


    Matrix x;
//    x = p.transpose();
//    x.printMatrix();
//    x = x0.transpose() * p;
//    x.printMatrix();
//    op.CGM();
    op.MPK();
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


    return 0;
}
