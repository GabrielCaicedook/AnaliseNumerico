cmake_minimum_required(VERSION 3.14)

# Especifica el estándar de C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


file(GLOB common_headers *.h json.hpp)
file(GLOB common_sources *.cpp)

add_executable( AnaliseNumerico main.cpp ${common_sources} ${common_headers})
