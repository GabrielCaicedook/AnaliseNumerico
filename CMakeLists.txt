cmake_minimum_required(VERSION 3.14)

# Especifica el estándar de C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

  project(GabrielAnalisis)
  
file(GLOB common_headers *.h json.hpp)
file(GLOB common_sources *.cpp)
#enabling Boost library
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)

# Especifica la versión mínima requerida de Boost
find_package(Boost 1.65.0 REQUIRED)

# Incluye los directorios de encabezado de Boost
include_directories(${Boost_INCLUDE_DIRS})



add_executable( AnaliseNumerico main.cpp ${common_sources} ${common_headers} )



#add_custom_target( GlobalFiles SOURCES ${PZ_SRC})



target_link_libraries(AnaliseNumerico ${Boost_LIBRARIES})



