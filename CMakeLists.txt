cmake_minimum_required(VERSION 3.9)
project(MO)

set(CMAKE_CXX_STANDARD 17)

set(Sources matrix.cpp parser.cpp)
set(Headers matrix.h parser.h)
add_executable(MO main.cpp ${Sources} ${Headers})