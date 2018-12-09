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

    expression = "x[0]^2+3*x[1]^2+2*x[0]*x[1]";
    Optimize op(expression);

    vector<double> p(2), x0(2);
    x0[0] = 1; x0[1] = 1;
    p[0] = 2; p[1] = 3;
    op.setX0(x0);
    op.setP(p);

    vector<double> x;

    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 5; j++) {
            cout << "Связка: " << op.getOneDimMethodByIndex(i) << " + " << op.getInterpMethodByIndex(j) << endl;
            op.setIndexQueue(i, j);
            x = op.newPoint(op.alphaSearch());
            cout << "Результат: " <<   x[0] << " " << x[1]   << endl;
        }

    }


    return 0;
}
