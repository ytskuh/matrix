#include <string.h>
#include <stdio.h>
#include "matrix.h"

int main() {
    mat_t a, b, c;
    mat_inits(4, 4, a, b, c, NULL);
    double f[4][4] = {
        {2, -1, 1, 3},
        {4, -1, 0, 5},
        {-4, 3, -3, -5},
        {2, 2, -3, 2}
    };
    mat_set_arr(a, (double*)f);
    __print_mat("%6.3f ", a);
    printf("\n");
    mat_inv_gauss(b, a);
    __print_mat("%6.3f ", b);
    printf("\n");
    mat_mul(c, a, b);
    __print_mat("%6.3f ", c);
    printf("\n");
    mat_lu(b, c, a);
    __print_mat("%6.3f ", b);
    printf("\n");
    __print_mat("%6.3f ", c);
    printf("\n");
    mat_mul(a, b, c);
    __print_mat("%6.3f ", a);
    printf("\n");

    double br[4] = {13, 17, -25, 0};
    double br2[4];
    mat_solve(br2, a, br);

    mat_clears(a, b, c, NULL);
}
