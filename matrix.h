/*
 * mat.h
 * 
 * This header file contains functions for matrix operations, including initialization, destruction, inversion, and multiplication.
 * 
 */

typedef unsigned int m_uint;

typedef struct {
    double **data;
    m_uint rows;
    m_uint cols;
} __mat_struct;

typedef __mat_struct mat_t[1];
typedef __mat_struct *mat_ptr;

// Function to initialize a matrix
void mat_init(mat_t rop, m_uint rows, m_uint cols);
void mat_inits(m_uint rows, m_uint cols, mat_ptr rop, ...);

void mat_set_zero(mat_t rop);
void mat_set_eye(mat_t rop);

void mat_set_arr(mat_t rop, const double *op);

// Function to destroy a matrix and free its memory
void mat_clear(mat_t rop);
void mat_clears(mat_ptr rop, ...);

void __print_mat(const char *fmt, const mat_t op);
int mat_printf(const char *fmt, ...);

void mat_copy_into(mat_t rop, const mat_t op);
void mat_copy_from_place(mat_t rop, const mat_t op, m_uint x, m_uint y);
void mat_swap_rows(mat_t rop, m_uint row1, m_uint row2);

int mat_transpose(mat_t rop, const mat_t op);
int mat_transpose_inplace(mat_t rop);

// Function to calculate the inverse of a matrix
int mat_inv_gauss(mat_t rop, const mat_t op);

int mat_lu(mat_t rop1, mat_t rop2, const mat_t op);
int mat_qr(mat_t rop1, mat_t rop2, const mat_t op);

int __mat_solve_assume_l(double *rop, const mat_t op1, const double *op2);
int __mat_solve_assume_u(double *rop, const mat_t op1, const double *op2);

int mat_solve(double *rop, const mat_t op1, const double *op2);

// Function to perform matrix multiplication
int mat_mul(mat_t rop, const mat_t op1, const mat_t op2);
void mat_mul_vec(double *rop, const mat_t op1, const double *op2);
