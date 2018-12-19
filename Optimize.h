//
// Created by Артем Калоев on 06/12/2018.
//

#ifndef MO_OPTIMIZE_H
#define MO_OPTIMIZE_H

#include <string>
#include <vector>
#include <iostream>
#include "Parser.h"
#include "Matrix.h"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

class Optimize {
private:
    string expression;
    Matrix x0;
    Matrix p;
    Matrix g;
    double alpha;

    //по умолчанию: ЗС-1 + Пауэлл
    int oneDimIndex = 6;
    int interpIndex = 3;
    string oneDimMethods[7] = {
            "ЗС-1",
            "ЗС-2",
            "Фибоначчи-1",
            "Фибоначчи-2",
            "Дихотомия",
            "Трехточеченый поиск",
            "Больцано"
    };
    string interpMethods[5] = {
            "Пауэлл",
            "Давидон",
            "ДСК",
            "Кубическая интерполяция",
            "Квадратичная интерполяция"
    };

    double e = 0.001;
    int kMaxSwann = 40;
    int kMaxOneDim = 5;
    int alphaCounter = 0;

    int dfOption = 1;
    double df(double);
    void df();

    void swann(double x, double &a, double &b, int &k);

    double ZS1(double &a, double &b, int &k);
    double ZS2(double &a, double &b, int &k);
    double Fib(int n);
    double fibonacci1(double &a, double &b, int &k);
    double fibonacci2(double &a, double &b, int &k);
    double dichtomy(double &a, double &b, int &k);
    double TPS(double &a, double &b, int &k);
    double bolcano(double &a, double &b, int &k);


    double Powell(double a, double b, int &k);
    double Davidon(double a, double b, int &k);
    void swannDSK(double x, double &a, double &b, double &c);
    double DSK(double x, double a, double b, int &k);
    double cubicInterpolation(double a, double b, int &k);
    double quadraticInterpolation(double a, double b, int &k);



    double eval(const Expression&);
    double f(double);

public:
    Optimize(const string& expr)    {   expression = expr;   }
    double y(Matrix x);
    void setOneDimIndex(int n)          {   oneDimIndex = n;    }
    void setInterpIndex(int n)          {   interpIndex = n;    }
    void setIndexQueue(int n, int m)    {   oneDimIndex = n; interpIndex = m;   }
    void getOneDimMethodsList();
    void getInterpMethodsList();
    string getOneDimMethodByIndex(int i)    {   return oneDimMethods[i];    };
    string getInterpMethodByIndex(int i)    {   return interpMethods[i];    };

    Matrix newPoint(double);
    double alphaSearch();

    void MPK();
    void CGM();
    void setMaxSwann(int n)         {   kMaxSwann = n;  }
    void setMaxOneDim(int n)        {   kMaxOneDim = n; }
    void setDfOption(int n)         {   dfOption = n;   }
    void setP(Matrix);
    void setX0(Matrix);
};


#endif //MO_OPTIMIZE_H
