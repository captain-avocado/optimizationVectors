//
// Created by Артем Калоев on 05/12/2018.
//

#ifndef MO_PARSER_H
#define MO_PARSER_H

#include <vector>
#include <iostream>
#include <iomanip>
using std::setw;

// Структура выражения
struct Expression {
    Expression(std::string token) : token(token) {} // Знаки
    Expression(std::string token, Expression a) : token(token), args{a} {} // Выражения с одним значением, например, -5
    Expression(std::string token, Expression a, Expression b) : token(token), args{a, b} {} // Обычные выражения

    std::string token;
    std::vector<Expression> args;
};

class Parser {
public:
    explicit Parser(const char* input) : input(input) {} // Конструктор, принимает строку с выражением
    Expression parse(); // Основная функция парсинга
private:
    int get_priority(const std::string& token);
    std::string parse_token(); // Парсит один токен
    Expression parse_simple_expression(); // Парсит простое выражение
    Expression parse_binary_expression(int min_priority); // Парсит бинарное выражение

    const char* input; // Кусок строки, который еще не распарсили
};


#endif //MO_PARSER_H
