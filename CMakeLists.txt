cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(BIE_laplace)

enable_language(C)
enable_language(Fortran)


set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}${GCC_COVERAGE_LINK_FLAGS}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/api)

find_package(vecLib REQUIRED)

include_directories(vecLib_INCLUDE_DIR)
include_directories(${PROJECT_SOURCE_DIR}/api)

add_executable(
main_exe
src/init.c
src/main.c
src/solveDensity.c
src/computeSolution.c
)


target_link_libraries(main_exe ${vecLib_LINKER_LIBS})
target_link_libraries(main_exe gsl ${GSL_LIBRARIES})
