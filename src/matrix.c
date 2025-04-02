#include <stdio.h>
#include <stdlib.h>
#include <mathlib.h>
#include <constants.h>

Matrix
create_matrix(size_t rows, size_t cols) {
	Matrix m;
	m.rows = rows;
	m.cols = cols;
	m.data = (double *)malloc(rows * cols * sizeof(double));
	return m;
}

void
set_matrix_value(const Matrix *m, size_t row, size_t col, double value) {
	m->data[row * m->cols + col] = value;
}

double
get_matrix_value(const Matrix *m, size_t row, size_t col) {
	return m->data[row * m->cols + col];
}

void 
init_matrix(Matrix *m, size_t rows, size_t cols) {
	size_t i;
	m->rows = rows;
	m->cols = cols;
	m->data = (double *)malloc(rows * cols * sizeof(double));

	if (m->data == NULL) {
		fprintf(stderr, "Error: unable to allocate memory for matrix\n");
		exit(1);
	}

	for (i = 0; i < rows * cols; ++i) {
		m->data[i] = 0.0;
	}
}

void 
free_matrix(Matrix *m) {
	if (m->data != NULL) {
		free(m->data);
		m->data = NULL;
	}
	m->rows = 0;
	m->cols = 0;
}

Matrix
zero_matrix(const Matrix *m) {
	size_t i, j;
	Matrix result;

	result = create_matrix(m->rows, m->cols);
	for (i = 0; i < m->rows; ++i) {
		for (j = 0; j < m->cols; ++j) {
			result.data[i * m->cols + j] = 0.0; 
		}
	}
	return result;
}

double
matrix_magnitude(const Matrix *m){
	double result;
	size_t i;

	for (i = 0; i < m->cols; ++i) {
		result += m->data[i] * m->data[i];
	}

	return sqrt_approx(result); 
}

double
homo_magnitude(const Matrix *m){
	double result;
	size_t i;

	for (i = 0; i < m->cols - 1; ++i) {
		result += m->data[i] * m->data[i];
	}

	return sqrt_approx(result); 
}

Matrix
scalar_mult(const Matrix *m, double n){
	size_t i, j;
	Matrix result;

	result = create_matrix(m->rows, m->cols);

	for (i = 0; i < m->rows; ++i) {
		for (j = 0; j < m->cols; ++j) {
			result.data[i * m->cols + j] = n * m->data[i * m->cols + j];
			if (result.data[i * m->cols + j] == -0.0f)
				result.data[i * m->cols + j] = 0.0f;
		}
	}
	return result;	
}

void
print_matrix(const Matrix *m) {
	size_t i, j;
	for (i = 0; i < m->rows; ++i) {
		for(j = 0; j <m->cols; j++) {
			printf("| %f |", m->data[i * m->cols + j]);
		}
		printf("\n");
	}
}

Matrix
matrix_add(const Matrix *m1, const Matrix *m2) {
	size_t i, j;
	Matrix result;
	if (m1->rows != m2->rows || m1->cols != m2->cols) {
		fprintf(stderr, "Error: Matrix dimensions do not match for addition.\n");
		exit(EXIT_FAILURE);
	}

	result = create_matrix(m1->rows, m1->cols);
	for (i = 0; i < m1->rows; ++i) {
		for (j = 0; j < m1->cols; ++j) {
			result.data[i * m1->cols + j] = m1->data[i * m1->cols + j] + m2->data[i* m1->cols + j];
		}
	}
	return result;
}


Matrix
matrix_sub(const Matrix *m1, const Matrix *m2) {
	size_t i, j;
	Matrix result;
	if (m1->rows != m2->rows || m1->cols != m2->cols) {
		fprintf(stderr, "Error: Matrix dimensions do not match for subtraction.\n");
		exit(EXIT_FAILURE);
	}

	result = create_matrix(m1->rows, m1->cols);
	for (i = 0; i < m1->rows; ++i) {
		for (j = 0; j < m1->cols; ++j) {
			result.data[i * m1->cols + j] = m1->data[i * m1->cols + j] - m2->data[i* m1->cols + j];
		}
	}
	return result;
}

Matrix
matrix_mult(const Matrix *m1, const Matrix *m2) {
	size_t i, j, k;
	Matrix result;

	if (m1->cols != m2->rows) {
		fprintf(stderr, "Error: The number of columns in the first matrix must equal the number of rows in the second. matrix_mult\n");
		exit(EXIT_FAILURE);
	}

	result = create_matrix(m1->rows, m2->cols);
	result = zero_matrix(&result); 
	for (i = 0; i < result.rows; ++i) {
		for (j = 0; j < result.cols; ++j) {
			for (k = 0; k < m1->cols; ++k) {
				result.data[i * result.cols + j] += m1->data[i * m1->cols + k] * m2->data[k * m2->cols + j];
			}
		}
	}						
	return result;
}


double
matrix_dot(const Matrix *m1, const Matrix *m2) {
	size_t i;
	double result;

	
	if(m1->cols != m2->cols) {
		fprintf(stderr, "Error: Matrix do not match size for dot product.\n");
		exit(EXIT_FAILURE);
	}	

	result = 0.0;	
	for (i = 0; i < m1->cols; ++i) {
		result = result + m1->data[i] * m2->data[i];
	}
	return result;
}

double
homo_matrix_dot(const Matrix *m1, const Matrix *m2) {
	size_t i;
	double result;

	
	if(m1->cols != m2->cols) {
		fprintf(stderr, "Error: Matrix do not match size for dot product.\n");
		exit(EXIT_FAILURE);
	}	

	result = 0.0;	
	for (i = 0; i < m1->cols - 1; ++i) {
		result = result + m1->data[i] * m2->data[i];
	}
	return result;
}

Matrix
transpose(const Matrix *m) {
	size_t i, j;
	Matrix result;

 	result = create_matrix(m->cols, m->rows);
	for (i = 0; i < m->rows; ++i) {
		for (j = 0; j < m->cols; ++j) {
			result.data[j * result.cols + i] = m->data[i * m->cols + j];
		}
	}
	return result;
}

Matrix
get_sub_matrix(const Matrix *m, size_t row_ex, size_t col_ex) {
	size_t i, j;
	int sub_i, sub_j;
	Matrix result;
	
	sub_i = 0;
	result = create_matrix(m->rows - 1, m->cols - 1);

	for (i = 0; i < m->rows; ++i) {
		if (i == row_ex) continue;

		sub_j = 0;
		for (j = 0; j < m->cols; ++j) {
			if (j == col_ex) continue;
			
			result.data[sub_i * result.cols + sub_j] = m->data[i * m->cols + j];
			sub_j++;
		} 
		sub_i++;
	}

	return result;
}

double 
get_sub_matrix_determinant(const Matrix *m) {
	double det;
	size_t j;

	if(m->rows != m->cols) {
		fprintf(stderr, "Matrix must be square to compute determinant.\n");
		exit(EXIT_FAILURE);
	}
	if (m->rows == 1) {
		return m->data[0];
	}
	if (m->rows == 2) {
		return m->data[0] * m->data[3] - m->data[1] * m->data[2];
	}

	det = 0.0;
	for (j = 0; j < m->cols; ++j) {
		Matrix cofactor_matrix = get_sub_matrix(m, 0, j);
		double cofactor = get_sub_matrix_determinant(&cofactor_matrix);
		free_matrix(&cofactor_matrix);
		det += (j % 2 == 0 ? 1 : -1) * m->data[j] * cofactor;
	}	
	return det;
}

double 
get_cofactor(const Matrix *m, size_t row, size_t col) {
	Matrix cofactor_matrix;
	double cofactor_value;
	
	cofactor_matrix = get_sub_matrix(m, row, col);
	cofactor_value = get_sub_matrix_determinant(&cofactor_matrix);
	free_matrix(&cofactor_matrix);
	return (row + col) % 2 == 0 ? cofactor_value : -cofactor_value;
}

double
get_determinant(const Matrix *m) {
	double det;
	size_t j, i;

	if(m->rows != m->cols) {
		fprintf(stderr, "Matrix must be square to compute determinant.\n");
		exit(EXIT_FAILURE);
	}
	if (m->rows == 1) {
		return m->data[0];
	}
	if (m->rows == 2) {
		return m->data[0] * m->data[3] - m->data[1] * m->data[2];
	}
	
	det = 0.0;
	for (j = 0; j < m->rows; ++j) {
		for (i = 0; i < m->cols; ++i) { 
			double cofactor_value = get_cofactor(m, j, i);
			det += m->data[j * m->cols + i] * cofactor_value;
		}
	}

	return det;
}

Matrix
invert(const Matrix *m) {
	double det;
	Matrix adj, result;

	det = get_determinant(m);
	if (det == 0) {
		fprintf(stderr, "Matrix is singular and cannot be inverted.\n");
		exit(EXIT_FAILURE);
	}

	adj = adjugate(m);
	result = scalar_mult(&adj, 1.0 / det);
	free_matrix(&adj);

	return result;
}

Matrix
get_cofactor_matrix(const Matrix *m){
	size_t i, j;
	double minor;
	Matrix result, sub;
	
	result = create_matrix(m->rows, m->cols);
	for (i = 0; i < result.rows; ++i) {
		for (j = 0; j < result.cols; ++j) {
			sub = get_sub_matrix(m, i, j);
			minor = get_determinant(&sub);
			result.data[i * m->cols + j] = ((i+j) % 2 == 0 ? 1 : -1) * minor;
		}
	}
	free_matrix(&sub);

	return result;
}

Matrix
adjugate(const Matrix *m) {
	Matrix result, cofac;
	cofac = get_cofactor_matrix(m);
	result = transpose(&cofac);
	free_matrix(&cofac);
	return result;
}

Matrix
ortho_proj_matrix(const Matrix *a) {
	Matrix gram, at, igram, result;
	
	at = transpose(a);
	gram = matrix_mult(&at, a);
	igram = invert(&gram);
	result = matrix_mult(a, &igram);
	result = matrix_mult(&result, &at);

	free_matrix(&at);
	free_matrix(&gram);
	free_matrix(&igram);
	return result;
}

Matrix
onto_proj(const Matrix *a, const Matrix *b) {
	Matrix result, ortho;

	ortho = ortho_proj_matrix(a);
	result = matrix_mult(&ortho, b);
	free_matrix(&ortho);

	return result;
}

Matrix
create_identity_matrix(size_t size) {
	Matrix result;
	size_t i, j;

	result = create_matrix(size, size);
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			if (i == j) {
				result.data[i * result.cols + j] = 1;
			} else {
				result.data[i * result.cols + j] = 0;
			}
		}
	}

	return result;
}

Matrix
create_translation_matrix(double tx, double ty, double tz) {
	Matrix result;
	result = create_identity_matrix(4);
	
	result.data[3] = tx;
	result.data[7] = ty;
	result.data[11] = tz;

	return result;

}

Matrix
create_scaling_matrix(double sx, double sy, double sz) {
	Matrix result;
	result = create_identity_matrix(4);

	result.data[0] = sx;
	result.data[5] = sy;
	result.data[10] = sz;

	return result;

}


Matrix 
matrix_normalize(const Matrix *m) {
	Matrix result;
	size_t i;
	double mag;

	result = create_matrix(1, 3);
	mag = matrix_magnitude(m); 
	
	for (i = 0; i < m->cols; ++i) {
		result.data[i] = m->data[i]/mag;
	}

	return result;
}

Matrix
create_rotation_matrix(double ux, double uy, double uz, double theta) {
	Matrix result;
       	Matrix norm_temp, norm;
	result = create_identity_matrix(4);
	double nux, nuy, nuz, uxs, uys, uzs, cos_theta, sin_theta, one_minus_cos; // squared values and theta values
	
	norm_temp = matrix_from_coords(ux, uy, uz);
	norm = matrix_normalize(&norm_temp);

	nux = norm.data[0];
	nuy = norm.data[1];
	nuz = norm.data[2];

	uxs = nux * nux;
	uys = nuy * nuy;
	uzs = nuz * nuz;

	cos_theta = mycos(theta, 0);
	sin_theta = mysin(theta, 0);
	one_minus_cos = 1 - cos_theta;
	
	free_matrix(&norm_temp);
	free_matrix(&norm);

	result.data[0] = cos_theta + uxs * one_minus_cos;
	result.data[1] = nux * nuy * one_minus_cos - nuz * sin_theta;
	result.data[2] = nux * nuz * one_minus_cos + nuy * sin_theta;
	result.data[4] = nuy * nux * one_minus_cos + nuz * sin_theta;
	result.data[5] = cos_theta + uys * one_minus_cos;
	result.data[6] = nuy * nuz * one_minus_cos - nux * sin_theta;
	result.data[8] = nuz * nux * one_minus_cos - nuy * sin_theta;
	result.data[9] = nuz * nuy * one_minus_cos + nux * sin_theta;
	result.data[10] = cos_theta + uzs * one_minus_cos;

	return result;

}

Matrix
transformation(const Matrix *t, const Matrix *s, const Matrix *r) {
	Matrix trs, result;

	trs = matrix_mult(r, s);
	result = matrix_mult(&trs, t);
	free_matrix(&trs);

	return result;
}

Matrix
matrix_from_vector(const Vector *v) {
	Matrix result;
	size_t i;
	result = create_matrix(1, v->size);
	
	for (i = 0; i < v->size; ++i) {
		result.data[i] = v->data[i];
	}

	return result;
}


Matrix
matrix_cross(const Matrix *m1, const Matrix *m2) {
	Matrix result;
	result = create_matrix(1, 3);

	result.data[0] = m1->data[1] * m2->data[2] - m1->data[2] * m2->data[1];
	result.data[1] = m1->data[2] * m2->data[0] - m1->data[0] * m2->data[2];
       	result.data[2] = m1->data[0] * m2->data[1] - m1->data[1] * m2->data[0];
	
	if (result.data[0] == -0.0f) 
		result.data[0] = 0.0f;
	if (result.data[1] == -0.0f) 
		result.data[1] = 0.0f;
	if (result.data[2] == -0.0f)
		result.data[2] = 0.0f;

	return result; 
}


Matrix
matrix_from_coords(double x, double y, double z){
	Matrix result;
	result = create_matrix(1, 3);

	result.data[0] = x;
	result.data[1] = y;
	result.data[2] = z;

	return result;
}


Matrix
create_matrix_vector(double x, double y, double z, double w){
	Matrix vec;
	vec = create_matrix(1, 4);

	vec.data[0] = x;
	vec.data[1] = y;
	vec.data[2] = z;
	vec.data[3] = w;

	return vec;
}

Matrix 
forward_matrix(const Matrix *t, const Matrix *c) {
        Matrix c_coords, t_coords, forward, temp, final; 
	double mag;
	
	c_coords = matrix_from_coords(c->data[0], c->data[1], c->data[2]);
	t_coords = matrix_from_coords(t->data[0], t->data[1], t->data[2]);
	
	forward = matrix_sub(&t_coords, &c_coords);
	
	mag = matrix_magnitude(&forward);
	temp = scalar_mult(&forward, 1 / mag);
	
	final = create_matrix_vector(temp.data[0], temp.data[1], temp.data[2], 0.0);

	free_matrix(&forward);
	free_matrix(&temp);
	free_matrix(&c_coords);
	free_matrix(&t_coords);

	return final;
}

Matrix
right_matrix(const Matrix *f, const Matrix *u) {
	double mag;
	Matrix right, temp, final, f_coords, u_coords;

	u_coords = matrix_from_coords(u->data[0], u->data[1], u->data[2]);
	f_coords = matrix_from_coords(f->data[0], f->data[1], f->data[2]);

	right = matrix_cross(&f_coords, &u_coords);
	mag = homo_magnitude(&right);
	temp = scalar_mult(&right, 1 / mag);

	final = create_matrix_vector(temp.data[0], temp.data[1], temp.data[2], 0.0);
	free_matrix(&right);
	free_matrix(&temp);
	free_matrix(&u_coords);
	free_matrix(&f_coords);

	return final; 
}

Matrix
up_matrix(const Matrix *r, const Matrix *f) {
	// updates up vector
	Matrix up, final, r_coords, f_coords;

	r_coords = matrix_from_coords(r->data[0], r->data[1], r->data[2]);
	f_coords = matrix_from_coords(f->data[0], f->data[1], f->data[2]);
	
	up = matrix_cross(&r_coords, &f_coords);
	final = create_matrix_vector(up.data[0], up.data[1], up.data[2], 0.0);
	
	free_matrix(&r_coords);
	free_matrix(&f_coords);
	free_matrix(&up);
	return final;
}

Matrix
sub_matrix_homogenous(const Matrix *m1, const Matrix *m2) {
	size_t i;
	Matrix result;
	
	result = create_matrix(1, 4);
	for (i = 0; i < m1->cols - 1; ++i) {
		result.data[i] = m1->data[i] - m2->data[i];
		if (result.data[i] == -0.0f) 
			result.data[i] = 0.0f;
	}
	result.data[3] = m1->data[3]; 
	return result;
}

Matrix
create_view_matrix(const Matrix *r, const Matrix *u, const Matrix *f, const Matrix *c) {
	Matrix view, r_mat, u_mat, f_mat; 
	double rc, uc, fc;
	
	r_mat = scalar_mult(r, -1);
	rc = homo_matrix_dot(&r_mat, c);
	
	u_mat = scalar_mult(u, -1);
	uc = homo_matrix_dot(&u_mat, c);
	
	f_mat = scalar_mult(f, -1);
	fc = homo_matrix_dot(f, c);

	view = create_identity_matrix(4);
	view.data[0] = r->data[0];
	view.data[1] = r->data[1];
	view.data[2] = r->data[2];
	view.data[3] = rc;
	view.data[4] = u->data[0];
	view.data[5] = u->data[1];
	view.data[6] = u->data[2];
	view.data[7] = uc;
	view.data[8] = f_mat.data[0];
	view.data[9] = f_mat.data[1];
	view.data[10] = f_mat.data[2];
	view.data[11] = fc;

	return view;
}

Matrix
create_perspective_matrix(double fov, double aspect, double near, double far) {
	Matrix perspective;
	double tan_fov_div_2;

	tan_fov_div_2 = mytan(fov / 2, 1);
	perspective = create_identity_matrix(4);

	perspective.data[0] = 1 / (aspect * tan_fov_div_2);
	perspective.data[5] = 1 / tan_fov_div_2;
	perspective.data[10] = (far + near) / (near - far);
	perspective.data[11] = (2 * far * near) / (near - far);
	perspective.data[14] = -1.0f;
	perspective.data[15] = 0.0f;

	return perspective;
}

Matrix
create_pndc(const Matrix *pclip) {
	Matrix pndc;

	pndc = matrix_from_coords(pclip->data[0] / pclip->data[3], pclip->data[1] / pclip->data[3], pclip->data[2] / pclip->data[3]);

	return pndc;
}

Matrix
create_pscreen(const Matrix *pndc, double height, double width) {
	Matrix pscreen;
	double xndc, yndc;

	xndc = (pndc->data[0] + 1) * width / 2;
	yndc = (pndc->data[1] + 1) * height / 2;

	pscreen = create_matrix(1, 2);
	pscreen.data[0] = xndc;
	pscreen.data[1] = yndc;

	return pscreen;
}
