//
// Created by Артем Калоев on 06/12/2018.
//

#include "Optimize.h"
#include <cmath>
#include <sstream>

Matrix Optimize::newPoint(double h) {
    return x0 + (p * h);
}

double Optimize::eval(const Expression& e) {
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

double Optimize:: y(Matrix x) {

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
                if (number >= x.getRows()) {
                    std::cout << "Invalid expression, no x[" << number << "]" << std::endl;
                    exit(1);
                }
                stream << x.values[0][number];
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

double Optimize::f(double h) {
    Matrix x = newPoint(h);
    return y(x);
}

double Optimize::df(double alpha){
    Matrix x = newPoint(alpha);
    g.values[0].resize(x.getRows());
    switch(dfOption){
        case 0: {
            for (int i = 0; i < x.getRows(); i++) {
                x.values[0][i] += e;
                g.values[0][i] = y(x);

                x.values[0][i] -= e;
                g.values[0][i] -= y(x);

                g.values[0][i] /= e;
            }
            break;
        }
        case 1: {
            for (int i = 0; i < x.getRows(); i++) {
                x.values[0][i] += e;
                g.values[0][i] = y(x);

                x.values[0][i] -= 2 * e;
                g.values[0][i] -= y(x);

                x.values[0][i] += e;
                g.values[0][i] /= 2 * e;
            }
            break;
        }
        case 2: {
            for (int i = 0; i < x.getRows(); i++) {
                x.values[0][i] -= e;
                g.values[0][i] = y(x);
                x.values[0][i] += e;
                g.values[0][i] -= 4 * y(x);

                x.values[0][i] += e;
                g.values[0][i] += 3 * y(x);

                x.values[0][i] -= e;
                g.values[0][i] /= 2 * e;
            }
            break;
        }
        case 3: {
            for (int i = 0; i < x.getRows(); i++) {
                x.values[0][i] += 2 * e;
                g.values[0][i] = -y(x);

                x.values[0][i] -= e;
                g.values[0][i] += 8 * y(x);

                x.values[0][i] -= 2 * e;
                g.values[0][i] -= 8 * y(x);

                x.values[0][i] -= e;
                g.values[0][i] += y(x);

                x.values[0][i] += 2 * e;
                g.values[0][i] /= 12 * e;
            }
            break;
        }
    }
    double result = 0;
    for (int i = 0; i < g.getRows(); i++) {
        result += g.values[0][i] * p.values[0][i];
    }
    return result;
}

void Optimize::swann(double x, double &a, double &b, int &k) {
    double x1 = x, x2, h = e;

    //изменить направление шага в случае возрастания функции
    if (f(x - h) < f(x)) { h = -h; }
    x2 = x + h;

    //пока функция убывает, перейти дальше, удвоив размер шага
    for (a = b = 0, k = 0; f(x1) > f(x2) && k < kMaxSwann; k++) {
        k++;
        h *= 2;
        x1 = x2;
        x2 += h;
    }

    //обозначить границы промежутка
    double x0 = x1 - h/2;
    if (x0 < x2) { a = x0; b = x2; }
    else { a = x2; b = x0; }
}

double Optimize::ZS1(double &a, double &b, int &k) {
    double phi = sqrt(5)/2 - 0.5,
            x1 = a + (1 - phi) * (b - a),
            x2 = a + phi * (b - a);
    for (k = 0; (b - a) > e && k < kMaxOneDim; k++) {
        if (f(x1) > f(x2)) {
            a = x1; x1= x2; x2 = a + phi * (b - a);
        } else {
            b = x2; x2 = x1; x1 = a + (1 - phi) * (b - a);
        }
    }
    return 0.5 * (a + b);
}

double Optimize::ZS2(double &a, double &b, int &k) {
    //найти первую точку деления промежутка
    double phi = sqrt(5)/2 - 0.5,
            x1 = a + phi * (b - a),
            x2;

    //пока промежуток недостаточно маленький
    for (k = 0; (b - a) > e && k < kMaxOneDim; k++) {

        //взять симметричную точку деления промежутка
        x2 = a + b - x1;

        //рассмотреть 4 ситуации, переместить одну из границ промежутка и пересчитать первую точку деления промежутка
        if (x2 <= x1) {
            if (f(x2) <= f(x1)) { b = x1; x1 = x2; }
            else { a = x2;  x1 = a + phi * (b - a); }
        } else {
            if (f(x2) <= f(x1)) { b = x2; x1 = a + phi * (b - a); }
            else { a = x1;  x1 = x2; }
        }
    }

    //вернуть минимум
    return 0.5 * (b + a);
}

double Optimize::Fib(int n) {
    if (n < 2) {
        if (n == 0) return 1;
        else return n;
    }
    return Fib(n - 1) + Fib(n - 2);
}

double Optimize::fibonacci1(double &a, double &b, int &k) {
    double x1, x2, Ln = 0.1 * e;

    int n = 0;
    while (Fib(n) < (b - a) / Ln) // вычисление числа n
    { n++; }


    x1 = a + (Fib(n - 2)*(b - a)) / Fib(n);
    x2 = a + (Fib(n - 1)*(b - a)) / Fib(n);

    for (k = 1; k < n - 1 && k < kMaxOneDim; k++) {
        if (f(x1) < f(x2)) {
            b = x2;
            x2 = x1;
            x1 = a + (Fib(n - k - 2)*(b - a)) / Fib(n - k);
        } else {
            a = x1;
            x1 = x2;
            x2 = a + (Fib(n - k - 1)*(b - a)) / Fib(n - k);
        }
    }

    if (f(x1) > f(x2)) {
        return 0.5 * (x1 + b);
    } else {
        return 0.5 * (a + x2);
    }
}

double Optimize::fibonacci2(double &a, double &b, int &k) {
    double x1, x2, Ln = e;

    double n = 0;
    while (Fib(n) < (b - a) / Ln) // вычисление числа n
    { n++; }

    x1 = a + (Fib(n - 1) / Fib(n)) * (b - a) + (pow((-1), n) / Fib(n)) * e;

    for (k = 1; k != n && k < kMaxOneDim; k++) {
        x2 = a + b - x1;
        if (f(x1) < f(x2)) {
            if (x1 < x2) {
                b = x2;
            }
            else {
                a = x2;
            }
        } else {
            if (x1 < x2) {
                a = x1;
                x1 = x2;
            }
            else {
                b = x1;
                x1 = x2;
            }
        }
    }
    return x2;
}

double Optimize::dichtomy(double &a, double &b, int &k) {
    double x1, x2, x3, eps = e / 2;

    for (k = 0; (b - a) > e && k < kMaxOneDim; k++) {
        x2 = (a + b) / 2;
        x1 = x2 - eps;
        x3 = x2 + eps;

        if (f(x1) < f(x3)) {
            b = x3;
        } else if (f(x1) > f(x3)) {
            a = x1;
        } else {
            return 0.5 * (a + b);
        }
    }

    return 0.5 * (a + b);
}

double Optimize::TPS(double &a, double &b, int &k) {
    double x2 = 0.5 * (a + b),
            x1 = 0.5 * (a + x2),
            x3 = 0.5 * (b + x2);

    for (k = 0; (b - a) > e && k < kMaxOneDim; k++) {
        if (f(x2) > f(x3)) {
            a = x2; x2 = x3;
        } else {
            if (f(x1) < f(x2)) {
                b = x2; x2 = x1;
            } else {
                a = x1; b = x3;
            }
        }

        x1 = 0.5 * (a + x2);
        x3 = 0.5 * (b + x2);
    }

    return 0.5 * (a + b);
}

double Optimize::bolcano(double &a, double &b, int &k) {
    double x1 = (a + b) / 2;
    k = 1;
    while (!((df(x1) <= e) && (abs(b - a) <= e))) {
        if (df(x1) > 0)
            b = x1;
        else a = x1;
        k++;
        if (k >= kMaxOneDim) {
            break;
        }
        x1 = (a + b) / 2;
    }
    return (a + b) / 2;
}

double Optimize::Powell(double a, double c, int &k) {
    double e1, e2;
    e1 = e2 = e;
    //найти точку b, как середину промежутка, и d по формуле
    double b = 0.5 * (a + c), d;
    d = 0.5 * (f(a) * (b * b - c * c) + f(b) * (c * c - a * a) + f(c) * (a * a - b * b))
        /(f(a) * (b - c) + f(b) * (c - a) + f(c) * (a - b));
    k = 1;
    //пока КОП не удовлетворен
    while(!(abs((b - d) / b) <= e1 &&
            abs((f(b) - f(d)) / f(b)) <= e2))
    {
        k++;

        //рассмотреть 4 ситуации, взять новую границу промежутка и одну из точек деления промежутка
        if (f(b) < f(d))
        {
            if (d < b) { a = d; }
            else { c = d; }
        } else {
            if (b < d) { a = b;
                b = d; }
            else { c = b; b = d; }
        }

        //найти d по формуле
        d = 0.5 * (a + b) + 0.5 *
                            ((f(a) - f(b)) * (b - c) * (c - a))
                            /(f(a) * (b - c) + f(b) * (c - a) + f(c) * (a - b));
    }

    //вернуть минимум
    return 0.5 * (b + d);
}

double Optimize::Davidon(double a, double b, int &k) {
    double gam, w, z, d = 0;
    k = 0;
    do {
        z = df(a) + df(b) +(3 * (f(a) - f(b)))/(b - a);
        w = sqrt(z * z - df(a) * df(b));
        gam = (w - df(a) + z) / (2 * w - df(a) + df(b));
        d = a + gam * (b - a);

        if (df(d) > 0) {
            b = d;
        } else {
            a = d;
        }

        k++;
    } while (df(d) > e);
    return d;
}

void Optimize::swannDSK(double x, double &a, double &b, double &c) {
    double x1 = x, x2, h = e;

    //изменить направление шага в случае возрастания функции
    if (f(x - h) < f(x)) { h = -h; }
    x2 = x + h;

    //пока функция убывает, перейти дальше, удвоив размер шага
    while (f(x1) > f(x2)) {
        h *= 2;
        x1 = x2;
        x2 += h;
    }

    //обозначить границы промежутка
    double x0 = x1 - h/2;
    if (x0 < x2) { a = x0; b = x2; }
    else { a = x2; b = x0; }

    c = b;
    b = (a + b) / 2;
    if (f(a) < f(b)) {
        b = a;
        a -= h;
    }
}

double Optimize::DSK(double x, double a, double b, int &k) {
    double c, d, h = e;
    swannDSK(x, a, b, c);
    d = b + 0.5*(pow(b - a, 2.0)*(f(b) - f(c)) - pow(b - c, 2.0)*(f(b) - f(a))) / ((b - a)*(f(b) - f(c)) - (b - c)*(f(b) - f(a)));

    for (k = 1; (abs(b - d) / b > e) || (abs(f(b) - f(d)) / f(b) > e); k++) {
        if (f(b) < f(d)) { x = b; }
        else { x = d; }
        h = h / 2;
        swannDSK(x, a, b, c);
        d = b + 0.5*(pow(b - a, 2.0)*(f(b) - f(c)) - pow(b - c, 2.0)*(f(b) - f(a))) / ((b - a)*(f(b) - f(c)) - (b - c)*(f(b) - f(a)));
    }
    return (b + d) / 2;
}

double Optimize::cubicInterpolation(double a, double b, int &k) {
    double gam, w, z, x;
    k = 0;
    do {
        z = df(a) + df(b) + (3 * (f(a) - f(b))) / (b - a);
        w = sqrt(z * z - df(a) * df(b));
        gam = (w - df(a) + z) / (2 * w - df(a) + df(b));

        if (gam > 1) { x = b; }
        else if (gam < 0) { x = a; }
        else x = a + gam * (b - a);

        if (df(x) > 0) { b = x; }
        else { a = x; }

        k++;

    } while ((b != x) && (a != x) && (df(x) > e));
    return x;
}

double Optimize::quadraticInterpolation(double a, double b, int &k) {
    double c = b, h, d;
    b = 0.5 * (a + c);
    if (c != 0) { h = e * fabs(c); }
    else { h = e; }
    d = 0.5 * (f(a)*(b*b - c*c) + f(b)*(c*c - a*a) + f(c)*(a*a - b*b)) / (f(a)*(b - c) + f(b)*(c - a) + f(c)*(a-b));
    for (k = 1; ((fabs((b - d) / b) > e) || (fabs((f(b) - f(d)) / f(b)) > e)); k++) {
        if (f(d) < f(b)) { b = d; }
        a = b - h/2;
        c = b + h/2;
        d = 0.5 * (f(a)*(b*b - c*c) + f(b)*(c*c - a*a) + f(c)*(a*a - b*b)) / (f(a)*(b - c) + f(b)*(c - a) + f(c)*(a-b));
    }

    return 0.5 * (b + d);
}

double Optimize::alphaSearch() {
    double a, b;
    int k = 0;
    alphaCounter = 0;

    swann(0, a, b, k);
    cout << "k1 = " << k << endl;
    alphaCounter += k;

    switch (oneDimIndex) {
        case 0: {
            alpha = ZS1(a, b, k);
            break;
        }
        case 1: {
            alpha = ZS2(a, b, k);
            break;
        }
        case 2: {
            alpha = fibonacci1(a, b, k);
            break;
        }
        case 3: {
            alpha = fibonacci2(a, b, k);
            break;
        }
        case 4: {
            alpha = dichtomy(a, b, k);
            break;
        }
        case 5: {
            alpha = TPS(a, b, k);
            break;
        }
        case 6: {
            alpha = bolcano(a, b, k);
            break;
        }
    }
    cout << "k2 = " << k << endl;
    alphaCounter += k;

    switch (interpIndex) {
        case 0: {
            alpha = Powell(a, b, k);
            break;
        }
        case 1: {
            alpha = Davidon(a, b, k);
            break;
        }
        case 2: {
            alpha = DSK(alpha, a, b, k);
            break;
        }
        case 3: {
            alpha = cubicInterpolation(a, b, k);
            break;
        }
        case 4: {
            alpha = quadraticInterpolation(a, b, k);
            break;
        }
    }
    cout << "k3 = " << k << endl;
    alphaCounter += k;

    cout << "Alpha = " << alpha << endl;
    cout << "AlphaCounter = " << alphaCounter << endl;
    return  alpha;
}

void Optimize::df() {
    for (int i = 0; i < x0.getRows(); i++) {
        x0.values[0][i] += e;
        g.values[0][i] = y(x0);

        x0.values[0][i] -= 2 * e;
        g.values[0][i] -= y(x0);

        x0.values[0][i] += e;
        g.values[0][i] /= 2 * e;
    }
}

void Optimize::MPK() {
    x0.values[0][0] = -2; x0.values[0][1] = -2;
    x0.values[0].resize(2);
    x0.values.resize(1);
    x0.cols = 1;
    x0.rows = 2;
    g.cols = 1;
    g.rows = 2;
    g.values[0].resize(2);
    g.values.resize(1);
    int j = 0, k;
    double alpha = 0, a, b;
    vector<double> vp = p.values[0], vg = g.values[0], vx0 = x0.values[0], vx = x0.values[0], prevx = x0.values[0], d = x0.values[0];
    double norma;
    do {
        for (int i = 0; i < vx0.size(); i++) {
            vx0[i] += e;
            vg[i] = y(Matrix(vx0));

            vx0[i] -= 2 * e;
            vg[i] -= y(Matrix(vx0));

            vx0[i] += e;
            vg[i] /= 2 * e;
        }

        x0 = Matrix(vx0);
        g = Matrix(vg);

        // Нахожу p, как антиградиент g
        for (int i = 0; i < vg.size(); i++) {
            vp[i] = -vg[i];
        }

        // Записываю p, чтобы другие методы работали нормально
        p = Matrix(vp);

        swann(0, a, b, k);
        alpha = ZS1(a, b, k);
        alpha = Powell(a, b, k);
//        alpha = alphaSearch();

        // При первой итерации нахожу вторую точку через вектор p и альфв
        if (j == 0) {
            for (int i = 0; i < vx0.size(); i++) {
                vx[i] = vx0[i] + vp[i] * alpha;
            }
        }

        for (int i = 0; i < vx0.size(); i++) {
            vx0[i] += e;
            vg[i] = y(Matrix(vx0));

            vx0[i] -= 2 * e;
            vg[i] -= y(Matrix(vx0));

            vx0[i] += e;
            vg[i] /= 2 * e;
        }

        x0 = Matrix(vx0);
        g = Matrix(vg);

        for (int i = 0; i < vg.size(); i++) {
            vp[i] = -vg[i];
        }

        p = Matrix(vp);

        swann(0, a, b, k);
        alpha = ZS1(a, b, k);
        alpha = Powell(a, b, k);

        // Запоминаю текущую точку, чтобы она стала пред-пред последней, чтобы потом захватить ее в x0
        prevx = vx;

        // Нахожу следующую (третью) точку для нахождения d
        for (int i = 0; i < vx0.size(); i++) {
            vx[i] = prevx[i] + vp[i] * alpha;
        }

        // Нахожу d, как разность первой и третьей точки
        for (int i = 0; i < vx.size(); i++) {
            d[i] = vx[i] - vx0[i];
        }

        swann(0, a, b, k);
        alpha = ZS1(a, b, k);
        alpha = Powell(a, b, k);

        // Нахожу четвертую точку, используя d, как нормирующее значение.
        for (int i = 0; i < d.size(); i++) {
            vx[i] = vx[i] + d[i] * alpha;
        }

        j++;

        norma = 0;
        for (int i = 0; i < d.size(); i++) {
            norma += pow(d[i], 2);
        }
        norma = sqrt(norma);

        vx0 = prevx;
    } while (norma > 0.001);

    x0.printMatrix();
    cout << "i = " << j << endl;
}

void Optimize::CGM() {

    //численное дифференцирование для n-переменных
    //выбирались альфа-методы
    //изменялись начальные точки в зависимости от функции
    //разные бетта-функции
    //можно сделать дифференцирование 2 и 4 способом

    x0.values[0][0] = -2; x0.values[0][1] = -2;
    x0.values[0].resize(2);
    x0.values.resize(1);
    x0.cols = 1;
    x0.rows = 2;
    g.cols = 1;
    g.rows = 2;
    g.values[0].resize(2);
    g.values.resize(1);
    int j = 0;
    double alpha = 0, betta1, betta2, betta3, betta4, a, b;
    int k = 0;
    Matrix prevG(x0.getRows()), gamma(x0.getRows());
    vector<double> vp = p.values[0], vg = g.values[0], vx0 = x0.values[0];
    double norma;
    do {
//        df();
        for (int i = 0; i < vx0.size(); i++) {
            vx0[i] += e;
            vg[i] = y(Matrix(vx0));

            vx0[i] -= 2 * e;
            vg[i] -= y(Matrix(vx0));

            vx0[i] += e;
            vg[i] /= 2 * e;
        }

        x0 = Matrix(vx0);
        g = Matrix(vg);

//         p = Matrix(g.getRows()) - g;

        for (int i = 0; i < vg.size(); i++) {
            vp[i] = -vg[i];
        }
        p = Matrix(vp);
//        if (j == 0) {
//            for (int i = 0; i < vg.size(); i++) {
//                vp[i] = -vg[i];
//            }
//            p = Matrix(vp);
////            p = Matrix(g.getRows()) - g;
//        } else {
//            betta1 = pow(g.norma(), 2) / pow(prevG.norma(), 2);
//            cout << "betta1 = " << betta1 << endl;
//            gamma = g - prevG;
//            betta2 = (g.transpose() * gamma).values[0][0] / (p.transpose() * gamma).values[0][0];
//            cout << "betta2 = " << betta2 << endl;
//            betta3 = (g.transpose() * gamma).values[0][0] / (prevG.transpose() * prevG).values[0][0];
//            cout << "betta3 = " << betta3 << endl;
//            betta4 = (g.transpose() * g).values[0][0] / (prevG.transpose() * prevG).values[0][0];
//            cout << "betta4 = " << betta4 << endl;
//            for (int i = 0; i < vg.size(); i++) {
//                vp[i] = -vg[i];
//            }
//            p = Matrix(vp);
////            p = Matrix(g.getRows()) - g + p * betta3;
//
//        }

//        prevG = g;

//        swann(0, a, b, k);
//        alpha = ZS1(a, b, k);
//        alpha = Powell(a, b, k);
        alpha = alphaSearch();

        for (int i = 0; i < vx0.size(); i++) {
            vx0[i] = vx0[i] + vp[i] * alpha;
        }
//        x0 = x0 + p * alpha;
        j++;

        norma = 0;
        for (int i = 0; i < vg.size(); i++) {
            norma += pow(vg[i], 2);
        }
        norma = sqrt(norma);
    } while (norma > 0.00001);
    x0.printMatrix();
    cout << "i = " << j << endl;

}

void Optimize::setX0(Matrix x) {
    x0.values[0].resize(x.getRows());
    for (int i = 0; i < x.getRows(); i++) {
        x0.values[0][i] = x.values[0][i];
    }
}

void Optimize::setP(Matrix x) {
    p.values[0].resize(x.getRows());
    for (int i = 0; i < x.getRows(); i++) {
        p.values[0][i] = x.values[0][i];
    }
}

void Optimize::getOneDimMethodsList() {
    cout << "Одномерные методы: " << endl;
    for (int i = 0; i < oneDimMethods->size(); i++) {
        cout << oneDimMethods[i] << endl;
    }
}

void Optimize::getInterpMethodsList() {
    cout << "Интерполяционные методы: " << endl;
    for (int i = 0; i < interpMethods->size(); i++) {
        cout << interpMethods[i] << endl;
    }
}

