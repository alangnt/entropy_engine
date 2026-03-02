cmake_minimum_required(VERSION 3.16)        # minimum cmake version required
project(entropy_engine LANGUAGES CXX)       # project name + skip C compiler check

set(CMAKE_CXX_STANDARD 17)                  # C++ version
set(CMAKE_CXX_STANDARD_REQUIRED ON)         # requires specific version

# help CMake find Homebrew's libomp on macOS
if(APPLE)
    execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
endif()

add_executable(universe main.cpp)           # builds an exe named "universe" from main.cpp
find_package(OpenMP REQUIRED)               # requires OpenMP installed
target_link_libraries(universe PRIVATE OpenMP::OpenMP_CXX)  # link OpenMP target
target_compile_options(universe PRIVATE -Wall -Wextra -Wpedantic -O3) # multithreading