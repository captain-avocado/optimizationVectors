#include <iostream>
#include <sstream>
#include <cmath>
#include "Matrix.h"
#include "Parser.h"

using std::cout;
using std::cin;
using std::endl;

std::string expression;

double eval(const Expression& e) {
    switch (e.args.size()) {
        case 2: {
            auto a = eval(e.args[0]);
            auto b = eval(e.args[1]);
            if (e.token == "+") return a + b;
            if (e.token == "-") return a - b;
            if (e.token == "*") return a * b;
            if (e.token == "/") return a / b;
            if (e.token == "^") return pow(a, b);
            if (e.token == "mod") return (int)a % (int)b;
            throw std::runtime_error("Unknown binary operator");
        }

        case 1: {
            auto a = eval(e.args[0]);
            if (e.token == "+") return +a;
            if (e.token == "-") return -a;
            if (e.token == "abs") return abs(a);
            if (e.token == "sin") return sin(a);
            if (e.token == "cos") return cos(a);
            throw std::runtime_error("Unknown unary operator");
        }

        case 0:
            return strtod(e.token.c_str(), nullptr);
    }

    throw std::runtime_error("Unknown expression type");
}

// Вычислить значение функции
double y(std::vector<double> x) {
    std::string tExpr = expression; // Временное выражение для обработки парсером
    bool flag = true;
    while (flag) {
        std::ostringstream stream;
        for (unsigned int j = 0; j < tExpr.length(); j++) {
            if (tExpr[j] == 'x') {
                std::string substr = "";
                int shift;
                for (int k = j + 2; tExpr[k] != ']'; k++) {
                    substr.push_back(tExpr[k]);
                    shift = k - j;
                }
                int number = stoi(substr);
                if (number >= x.size()) {
                    std::cout << "Invalid expression, no x[" << number << "]" << std::endl;
                    exit(1);
                }
                stream << x[number];
                tExpr.erase(j, 2 + shift);
                tExpr.insert(j, stream.str());
                break;
            } else if (j == tExpr.length() - 1) {
                flag = false;
            }
        }
    }
    char *c = new char[tExpr.size() + 1];
    std::copy(tExpr.begin(), tExpr.end(), c);
    c[tExpr.size()] = '\0';
    Parser p(c);
    return eval(p.parse());
}

int main() {

    std::cout << "Введите строку:" << std::endl;
    cin >> expression;

    vector<double> x(2);
    x[0] = 1;
    x[1] = 2;

    cout << "Результат: " << y(x) << endl;

//    Matrix matr;
//    matr.printMatrix();

    return 0;
}
