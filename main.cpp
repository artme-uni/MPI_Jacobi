#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "Jacobi.h"

int main(int argc, char **argv)
{
    int pr_rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &pr_rank);

    if (NX % comm_size != 0)
    {
        std::cerr << ("Grid size should be multiple of number of processes\n");
        MPI_Finalize();
        return 0;
    }

    int pr_cells_count = (NX / comm_size);
    int pr_cells_shift = pr_rank * pr_cells_count - 1;

    double *grid = (double *) malloc(sizeof(double) * (pr_cells_count + 2) * NY * NZ);

    double t_start = MPI_Wtime();
    int iteration_count = jacobi_method(grid, pr_cells_count, pr_cells_shift);
    double t_end = MPI_Wtime();

    if (pr_rank == 0)
    {
        printf("Iteration count: %d\n", iteration_count);
        printf("Time: %f\n", t_end - t_start);
        grid_print(1, grid);
    }

    free(grid);

    MPI_Finalize();
    return 0;
}


