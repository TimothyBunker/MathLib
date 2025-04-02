#ifndef MATHLIB_H
#define MATHLIB_H

#include <stddef.h>

// Vector Struct
typedef struct {
	size_t size;
	double *data;
} Vector;

typedef struct {
	double x, y, z, w;
} Vector4;

typedef struct {
	double x, y, z;
} Vector3;

typedef struct {
	double x, y;
} Vector2;

// Matrix Struct

typedef struct {
	size_t rows;
	size_t cols;
	double *data;
} Matrix;

// General Math

double sqrt_approx(double n);
double mysin(double x, int rad);
double mycos(double x, int rad);
double mytan(double x, int rad);
double myasin(double x);
double myacos(double x);
double myatan(double x);
double expo(double base, int n);
double fact(int n);
double rad_to_theta(double r);
double theta_to_rad(double t);
double reduce_angle(double r);
int my_abs(int value);

// Vector Operations
Vector create_vector(size_t size);
void free_vector(Vector *v);
void print_vector(const Vector *v);
Vector vector_add(const Vector *v1, const Vector *v2);
Vector vector_sub(const Vector *v1, const Vector *v2);
double dot(const Vector *v1, const Vector *v2);
Vector zero_vector(const Vector *v);
Vector proj_vector(const Vector *v1, const Vector *v2);
Vector mult_vector(double n, const Vector *v);
Vector vector_from_coords(double x, double y, double z);
double magnitude(const Vector *v);
Vector normalize(const Vector *v);
Vector cross(const Vector *v1, const Vector *v2);
Vector row_vector(const Matrix *m, int row);

// Matrix Operations
Matrix create_matrix(size_t rows, size_t cols);
void free_matrix(Matrix *m);
void print_matrix(const Matrix *m);
Matrix matrix_add(const Matrix *m1, const Matrix *m2);
Matrix matrix_sub(const Matrix *m1, const Matrix *m2);
Matrix matrix_mult(const Matrix *m1, const Matrix *m2);
Matrix zero_matrix(const Matrix *m);
Matrix transpose(const Matrix *m);
Matrix get_sub_matrix(const Matrix *m, size_t row_ex, size_t col_ex);
double get_cofactor(const Matrix *m, size_t row, size_t col);
double get_sub_matrix_determinant(const Matrix *m);
double get_determinant(const Matrix *m);
Matrix invert(const Matrix *m);
Matrix get_cofactor_matrix(const Matrix *m);
Matrix scalar_mult(const Matrix *m, double n);
Matrix adjugate(const Matrix *m);
Matrix ortho_proj_matrix(const Matrix *a);
Matrix onto_proj(const Matrix *a, const Matrix *b);
Matrix create_identity_matrix(size_t size);
Matrix create_translation_matrix(double tx, double ty, double tz);
Matrix create_rotation_matrix(double ux, double uy, double uz, double theta);
Matrix create_scaling_matrix(double sx, double sy, double sz);
Matrix create_perspective_matrix(double fov, double aspect, double near, double far);
Matrix create_orthographic_matrix(double left, double right, double bottom, double top, double near, double far);
Matrix transformation(const Matrix *t, const Matrix *s, const Matrix *r);
Matrix matrix_from_vector(const Vector *v);                                  // There functions are essentially all vectors just stored in a Matrix struct lol
Matrix create_matrix_vector(double x, double y, double z, double w);
Matrix create_view_matrix(const Matrix *r, const Matrix *u, const Matrix *f, const Matrix *c);
Matrix forward_matrix(const Matrix *t, const Matrix *c);
Matrix right_matrix(const Matrix *f, const Matrix *u);
Matrix up_matrix(const Matrix *r, const Matrix *f);
Matrix matrix_cross(const Matrix *m1, const Matrix *m2);
Matrix matrix_from_coords(double x, double y, double z);
double matrix_magnitude(const Matrix *m);
double matrix_dot(const Matrix *m1, const Matrix *m2);
Matrix create_pndc(const Matrix *pclip);
Matrix create_pscreen(const Matrix *pndc, double width, double height);
void set_matrix_value(const Matrix *m, size_t row, size_t col, double value);
double get_matrix_value(const Matrix *m, size_t row, size_t col);
Matrix sub_matrix_homogenous(const Matrix *m1, const Matrix *m2);
double homo_magnitude(const Matrix *m);
double homo_matrix_dot(const Matrix *m1, const Matrix *m2);
Matrix matrix_normalize(const Matrix *m);
void init_matrix(Matrix *m, size_t rows, size_t cols);

#endif // MATHLIB_H
