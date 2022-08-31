cmake_minimum_required(VERSION 3.23)
project(bivariate_interpolation)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(boost)
include_directories(boost/math)
include_directories(boost/math/interpolators)
include_directories(boost/math/interpolators/detail)

add_executable(bivariate_interpolation
        boost/math/interpolators/detail/triangulation.hpp)
