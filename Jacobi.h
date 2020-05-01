#ifndef MPI_JACOBI_JACOBI_H
#define MPI_JACOBI_JACOBI_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "baseline_data.h"

int jacobi_method(double *grid, int pr_cells_count, int pr_cells_shift);

void update_grid_cell(int index, double *grid, double *previous_grid, int pr_cells_shift, double *pr_diff);

void update_grid_center(double *grid, double *previous_grid, int pr_cells_count, int pr_cells_shift, double *pr_diff);

void update_grid_bound(double *grid, double *previous_grid, int pr_cells_count, int pr_cells_shift, double *pr_diff);

void grid_init(double *grid, double *previous_grid, int pr_cells_count, int pr_cells_shift);

void send_grid_bound(double *grid, int pr_rank, int comm_size, int pr_cells_count, MPI_Request *request_prev, MPI_Request *request_next);

void recieve_grid_bound(int pr_rank, int comm_size, MPI_Request *request_prev, MPI_Request *request_next);

double iteration_func(int i, int j, int k, double* grid, int pr_cells_shift);

void grid_print(int layer_index, double *grid);

#endif //MPI_JACOBI_JACOBI_H
