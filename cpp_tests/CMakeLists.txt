cmake_minimum_required (VERSION 2.6)
project (SolverCuda)
enable_language(Fortran)

add_library(read_matrix read_matrix.f90)
add_executable(solver_cuda solver_cuda.cpp read_matrix.f90) 
#target_link_libraries(solver_cuda read_matrix)

