#ifndef MPI_JACOBI_BASELINE_DATA_H
#define MPI_JACOBI_BASELINE_DATA_H

//Область моделирования
#define X0 -1
#define X1 1
#define Y0 -1
#define Y1 1
#define Z0 -1
#define Z1 1

//Параметр уравнения
#define A (1e5)

//Порог сходимости
#define E (1e-8)

//Начальное приближение
#define PHI0 0

//Размеры сетки
#define NX 100
#define NY 100
#define NZ 100

//Шаги сетки
#define hx ((double)(X1-X0)/(NX-1))
#define hy ((double)(Y1-Y0)/(NY-1))
#define hz ((double)(Z1-Z0)/(NZ-1))

//Искомая функция (зависимость от i,j,k)
double phi(double i, double j, double k);

//Правая часть уравнения (зависимость от i,j,k)
double ro(double i, double j, double k);

#define I(i, j, k) ((i)*NY*NZ+(j)*NZ+(k))

#endif //MPI_JACOBI_BASELINE_DATA_H
