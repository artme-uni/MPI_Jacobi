cmake_minimum_required(VERSION 3.16)
project(MPI_Jacobi)

set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

include_directories(/usr/local/opt/mpich/include)
link_directories(/usr/local/opt/mpich/lib)

add_executable(MPI_Jacobi main.cpp Jacobi.cpp Jacobi.h baseline_data.cpp baseline_data.h)
target_link_libraries(MPI_Jacobi mpi)




