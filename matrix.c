#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include "matrix.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

void mat_init(mat_t rop, m_uint rows, m_uint cols) {
    rop->cols = cols;
    rop->rows = rows;
    rop->data = (double**)malloc(rows * sizeof(double*));
    for (m_uint i = 0; i < rows; i++) {
        rop->data[i] = (double*)malloc(cols * sizeof(double));
        memset(rop->data[i], 0, cols * sizeof(double));
    }
}

void mat_inits(m_uint rows, m_uint cols, mat_ptr rop, ...) {
    va_list ap;
    va_start(ap, rop);
    do {
        mat_init(rop, rows, cols);
        rop = va_arg(ap, mat_ptr);
    } while(rop != NULL);
    va_end(ap);
}

void mat_set_zero(mat_t rop) {
    for(m_uint i = 0; i < rop->rows; i++)
        memset(rop->data[i], 0, rop->cols * sizeof(double));
}

void mat_set_eye(mat_t rop) {
    mat_set_zero(rop);
    for(m_uint i = 0; i < min(rop->cols, rop->rows); i++)
        rop->data[i][i] = 1.0;
}

void mat_set_arr(mat_t rop, const double *op) {
    for (m_uint i = 0; i < rop->rows; i++) {
        memcpy(rop->data[i], op +  i*rop->cols, rop->cols * sizeof(double));
    }
}

// Function to destroy a rop and free its memory
void mat_clear(mat_t rop) {
    for (m_uint i = 0; i < rop->rows; i++)
        free(rop->data[i]);
    free(rop->data);
}

void mat_clears(mat_ptr rop, ...) {
    va_list ap;
    va_start(ap, rop);
    do {
        mat_clear(rop);
        rop = va_arg(ap, mat_ptr);
    } while(rop != NULL);
    va_end(ap);
}

void __print_mat(const char *fmt, const mat_t op) {
    for (m_uint i = 0; i < op->rows; i++) {
        for (m_uint j = 0; j < op->cols; j++) {
            printf(fmt, op->data[i][j]);
        }
        printf("\n");
    }
}

void mat_copy_into(mat_t rop, const mat_t op) {
    for(m_uint i = 0; i < op->rows; i++)
        memcpy(rop->data[i], op->data[i], op->cols * sizeof(double));
}

void mat_copy_from_place(mat_t rop, const mat_t op, m_uint x, m_uint y) {
    unsigned int width = min(rop->cols, op->cols - y);
    unsigned int height = min(rop->rows, op->rows - x);
    for (m_uint i = 0; i < height; i++)
        memcpy(rop->data[i], op->data[i+x] + y, width * sizeof(double));
}

// Swap two rows of a matrix
void mat_swap_rows(mat_t rop, m_uint row1, m_uint row2) {
    double* temp = rop->data[row1];
    rop->data[row1] = rop->data[row2];
    rop->data[row2] = temp;
}

int mat_transpose(mat_t rop, const mat_t op) {
    if(rop->rows != op-> cols || rop->cols != op->rows) return 1;
    for(m_uint i = 0; i < op->rows; i++)
        for(m_uint j = 0; j < op->cols; j++)
            rop->data[j][i] = op->data[i][j];
    return 0;
}

int mat_transpose_inplace(mat_t rop) {
    if(rop->rows != rop-> cols) return 1;
    for (m_uint i = 0; i < rop->rows; i++) {
        for (m_uint j = i + 1; j < rop->cols; j++) {
            double temp = rop->data[i][j];
            rop->data[i][j] = rop->data[j][i];
            rop->data[j][i] = temp;
        }
    }
    return 0;
}

// Inverse of a matrix using Gauss-Jordan elimination
int mat_inv_gauss(mat_t rop, const mat_t op) {
    // Check the size
    if (op->rows != op->cols)
        return 1;
    if (op->rows != rop->rows || op->cols != rop->cols)
        return 2;
    
    m_uint n = op->rows;
    mat_t temp;
    mat_init(temp, n, 2 * n);
    
    // Initialize the augmented matrix with the original matrix and identity matrix
    mat_copy_into(temp, op);
    for(m_uint i=0; i<n; i++)
        temp->data[i][i + n] = 1.0;
    
    // Apply Gauss-Jordan elimination
    for (m_uint i = 0; i < n; i++) {
        // Find the pivot element
        m_uint pivot_row = i;
        double pivot = temp->data[i][i];
        for (m_uint j = i + 1; j < n; j++) {
            if (fabs(temp->data[j][i]) > fabs(pivot)) {
                pivot_row = j;
                pivot = temp->data[j][i];
            }
        }
        
        // Swap rows if necessary
        if (pivot_row != i)
            mat_swap_rows(temp, i, pivot_row);
        
        // Divide the pivot row by the pivot element
        for (m_uint j = 0; j < 2 * n; j++) {
            temp->data[i][j] /= pivot;
        }
        
        // Eliminate other rows
        for (m_uint j = 0; j < n; j++) {
            if (j != i) {
                double factor = temp->data[j][i];
                for (m_uint k = 0; k < 2 * n; k++) {
                    temp->data[j][k] -= factor * temp->data[i][k];
                }
            }
        }
    }
    // Extract the inverse matrix from the augmented matrix
    mat_copy_from_place(rop, temp, 0, n);
    
    mat_clear(temp);
    return 0;
}

int mat_lu(mat_t rop1, mat_t rop2, const mat_t op) {
    // Check if the input matrix is square
    if (op->rows != op->cols)
        return 1;  // LU decomposition is only applicable to square matrices
    if (rop1->rows != rop1->cols || rop2->rows != rop2->cols)
        return 2;  // rop1 and rop2 must be square matrices
    if (rop1->rows != op->rows || rop2->rows != op->rows)
        return 3;  // rop1 and rop2 must have the same size as op

    m_uint n = op->rows;
    mat_set_zero(rop1);
    mat_set_zero(rop2);

    // Perform LU decomposition
    for (m_uint k = 0; k < n; k++) {
        // Compute the elements of the lower triangular matrix
        for (m_uint i = k; i < n; i++) {
            double sum = 0.0;
            for (m_uint p = 0; p < k; p++)
                sum += rop1->data[i][p] * rop2->data[k][p];
            rop1->data[i][k] = op->data[i][k] - sum;
        }

        // Compute the elements of the upper triangular matrix
        for (m_uint j = k + 1; j < n; j++) {
            double sum = 0.0;
            for (m_uint p = 0; p < k; p++)
                sum += rop1->data[k][p] * rop2->data[j][p];
            rop2->data[j][k] = (op->data[k][j] - sum) / rop1->data[k][k];
        }
        rop2->data[k][k] = 1.0;
    }
    mat_transpose_inplace(rop2);
    return 0;  // LU decomposition successful
}

double __norm(const double* vec, const m_uint n) {
    double sum = 0.0;
    for (m_uint i = 0; i < n; i++)
        sum += vec[i] * vec[i];
    return sqrt(sum);
}

// QR decomposition method
int mat_qr(mat_t rop1, mat_t rop2, const mat_t op) {
    m_uint m = op->rows;
    m_uint n = op->cols;
    if (rop1->rows != m || rop1->cols != m)
        return 1;
    if (rop2->rows != m || rop2->cols != n)
        return 2;
    
    mat_t temp1, temp2, temp3;
    mat_init(temp1, n, n);
    mat_init(temp2, n, m);
    mat_init(temp3, n, m);
    mat_set_eye(rop1);
    mat_transpose(temp3, op);
    
    double *u = (double*)malloc(n * sizeof(double));
    double s;

    for (m_uint k=0; k<n; k++) {
        memset(u, 0, n * sizeof(double));
        mat_set_eye(temp1);
        u[k] = temp3->data[k][k] + __norm(temp3->data[k] + k, m - k);
        memcpy(u+k+1, temp3->data[k]+k+1, (m-k-1) * sizeof(double));
        s = __norm(u+k, m-k);
        for (m_uint i=k; i<m; i++)
            u[i] /= s;
        for (m_uint i=k; i<n; i++) {
            temp1->data[i][i] -= 2*u[i]*u[i];
            for (m_uint j=i+1; j<n; j++) {
                temp1->data[i][j] -= 2*u[i]*u[j];
                temp1->data[j][i] -= 2*u[i]*u[j];
            }
        }
        mat_mul(temp2, temp3, temp1);
        mat_copy_into(temp3, temp2);
        mat_mul(temp2, temp1, rop1);
        mat_copy_into(rop1, temp2);
    }
    mat_transpose(rop2, temp3);
    mat_transpose_inplace(rop1);
    mat_clears(temp1, temp2, temp3, NULL);
    return 0;
}

int __mat_solve_assume_l(double *rop, const mat_t op1, const double *op2) {
    m_uint n = op1->rows;

    // Solve the linear system assuming op1 is a lower triangular matrix
    for (m_uint i = 0; i < n; i++) {
        double sum = 0.0;
        for (m_uint j = 0; j < i; j++)
            sum += op1->data[i][j] * rop[j];
        rop[i] = (op2[i] - sum) / op1->data[i][i];
    }

    return 0;  // Linear system solved successfully
}

int __mat_solve_assume_u(double *rop, const mat_t op1, const double *op2) {
    m_uint n = op1->rows;

    // Solve the linear system assuming op1 is an upper triangular matrix
    for (m_uint i = n - 1; i < n; i--) {
        double sum = 0.0;
        for (m_uint j = i + 1; j < n; j++)
            sum += op1->data[i][j] * rop[j];
        rop[i] = (op2[i] - sum) / op1->data[i][i];
    }

    return 0;  // Linear system solved successfully
}

int mat_solve(double *rop, const mat_t op1, const double *op2) {
    m_uint n = op1->rows;

    mat_t l, u;
    mat_init(l, n, n);
    mat_init(u, n, n);

    // Perform LU decomposition
    mat_lu(l, u, op1);

    // Solve the equation
    __mat_solve_assume_l(rop, l, op2);
    __mat_solve_assume_u(rop, u, rop);

    mat_clear(l);
    mat_clear(u);
    return 0;
}

// Function to perform matrix multiplication
int mat_mul(mat_t rop, const mat_t op1, const mat_t op2) {
    // Check if the matrices can be multiplied
    if (op1->cols != op2->rows) return 1;
    if (rop->rows != op1->rows || rop->cols != op2->cols) return 2;

    mat_t op2t;
    mat_init(op2t, op2->cols, op2->rows);
    mat_transpose(op2t, op2);

    for (m_uint i = 0; i < op1->rows; i++) {
        for (m_uint j = 0; j < op2->cols; j++) {
            double sum = 0.0;
            for (m_uint k = 0; k < op1->cols; k++)
                sum += op1->data[i][k] * op2t->data[j][k];
            rop->data[i][j] = sum;
        }
    }
    mat_clear(op2t);
    return 0;
}

void mat_mul_vec(double *rop, const mat_t op1, const double *op2) {
    m_uint rows1 = op1->rows;
    m_uint cols1 = op1->cols;

    for (m_uint i = 0; i < rows1; i++) {
        double sum = 0.0;
        for (m_uint j = 0; j < cols1; j++)
            sum += op1->data[i][j] * op2[j];
        rop[i] = sum;
    }
}

#undef min