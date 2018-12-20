cmake_minimum_required(VERSION 3.8.2)
project(MO)

set(CMAKE_CXX_STANDARD 17)

set(Sources matrix.cpp parser.cpp optimize.cpp Vector.h)
set(Headers matrix.h parser.h optimize.h Vector.cpp)
add_executable(MO main.cpp ${Sources} ${Headers})