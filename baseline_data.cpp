#include "baseline_data.h"

double phi(double i, double j, double k)
{
    double x = X0 + (i) * hx;
    double y = X0 + (j) * hy;
    double z = X0 + (k) * hz;
    return x * x + y * y + z * z;
}

double ro(double i, double j, double k)
{
    return 77 - A * phi(i, j, k);
}
