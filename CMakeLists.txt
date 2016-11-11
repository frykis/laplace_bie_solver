cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(BIE_laplace)

enable_language(C)
enable_language(Fortran)


set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}${GCC_COVERAGE_LINK_FLAGS}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(GSL_DIR ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/pdc/vol/gsl/1.16/GNU/4.9/bin")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/api)

find_package(vecLib REQUIRED)
#find_package(GSL REQUIRED)

include_directories(vecLib_INCLUDE_DIR)
include_directories(${PROJECT_SOURCE_DIR}/api)

add_executable(
main_exe
src/main.c
src/init.c
src/initFunction.c
src/solveDensity_lapack.c
src/computeSolution.c
src/specialquadlapl.c
src/computeError.c
)


target_link_libraries(main_exe ${vecLib_LINKER_LIBS})
target_link_libraries(main_exe gsl ${GSL_LIBRARIES})
