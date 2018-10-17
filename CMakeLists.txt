cmake_minimum_required(VERSION 3.12)
project(ISM_Lab5)

set(CMAKE_CXX_STANDARD 17)

add_executable(ISM_Lab5 main.cpp
                        Cloneable.h
                        PRNG.h
                        MultiplicativePRNG.h
                        MultiplicativePRNG.cpp)