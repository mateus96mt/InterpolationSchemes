cmake_minimum_required(VERSION 3.17)
project(c_plu_plus_solver)

set(CMAKE_CXX_STANDARD 14)

add_executable(c_plu_plus_solver src/main.cpp include/UpwindSchemes.h include/Tools.h src/Tools.cpp include/Solver.h src/Solver.cpp include/Aux.h)