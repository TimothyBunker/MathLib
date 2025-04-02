#include <stdio.h>
#include "mathlib.h"

int
main() {
	// Test Vector Operations
	Vector v1, v2, v3, v4, v5, proj, norm, cp, from_coords;
	double d, m;

	v1 = create_vector(3);
	v2 = create_vector(3);

	v1.data[0] = 1.0; v1.data[1] = 2.0; v1.data[2] = 3.0;
	v2.data[0] = 4.0; v2.data[1] = 5.0; v2.data[2] = 6.0;

	from_coords = vector_from_coords(3.0, 2.5, 1.9);
	norm = normalize(&from_coords);
	cp = cross(&v1, &v2);

	printf("\n----------\n");
	printf("\nVector v1:\n");
	print_vector(&v1);
	printf("\nVector v2:\n");
	print_vector(&v2);

	v3 = vector_add(&v1, &v2);
	printf("\nVector v1 + v2:\n");
	print_vector(&v3);

	d = dot(&v1, &v2);
	printf("\nDot Product v2 and v2:\n%f\n", d);

	v4 = zero_vector(&v1);
	printf("\nZero'd v1 vector:\n");
	print_vector(&v4);

	v5 = mult_vector(3.0, &v1);
	printf("\nVector 3 * v1:\n");
	print_vector(&v5);


	proj = proj_vector(&v1, &v2);
	printf("\nVector projection v1 onto v2:\n");
	print_vector(&proj);

	printf("\nVector from coords:\n");
	print_vector(&from_coords);

	printf("\nNormalized vector from_coords:\n");
	print_vector(&norm);

	printf("\ncross_product of v1 and v2:\n");
	print_vector(&cp);

	m = magnitude(&v1);
	printf("\nMagnitude v1:\n%f\n", m);
	printf("\n----------\n");

	free_vector(&v1);
	free_vector(&v2);
	free_vector(&v3);
	free_vector(&v4);
	free_vector(&v5);
	free_vector(&proj);
	free_vector(&cp);
	free_vector(&norm);
	free_vector(&from_coords);

	Matrix m1, m2, m3, m4, m5, m6, m7;
	m1 = create_matrix(2, 2);
	m2 = create_matrix(2, 2);
	m6 = create_matrix(3, 3);

	m1.data[0] = 1.0; m1.data[1] = 2.0; m1.data[2] = 3.0; m1.data[3] = 4.0;
	m2.data[0] = 5.0; m2.data[1] = 6.0; m2.data[2] = 7.0; m2.data[3] = 8.0;
	m6.data[0] = 1.0; m6.data[1] = 2.0; m6.data[2] = 3.0;
	m6.data[3] = 4.0; m6.data[4] = 5.0; m6.data[5] = 6.0;
	m6.data[6] = 7.0; m6.data[7] = 8.0; m6.data[8] = 9.0;

	printf("\nMatrix 1:\n");
	print_matrix(&m1);
	printf("\nMatrix 2:\n");
	print_matrix(&m2);

	m3 = matrix_add(&m1, &m2);
	printf("\nMatrix m1 + m2:\n");
	print_matrix(&m3);

	m4 = matrix_mult(&m1, &m2);
	printf("\nMatrix m1 x m2:\n");
	print_matrix(&m4);
	

	m5 = transpose(&m1);
	printf("\nMatrix before Transposition:\n");
	print_matrix(&m1);
	printf("\nTransposed m1:\n");
	print_matrix(&m5);

	m7 = get_sub_matrix(&m6, 2, 2);
	printf("\nCofactor of 3x3 Matrix m6:\n");
	print_matrix(&m7);
	
	double sub_det, det;
	Matrix inverse, projm, m8, m9, identity, trans;
	sub_det = get_sub_matrix_determinant(&m7);
	det = get_determinant(&m6);
	inverse = invert(&m2);
	printf("\nSub Determinant m7: %1f\n", sub_det);
	printf("\nDeterminant m7: %1f\n", det);
	printf("\nInverse m2:\n");
	print_matrix(&inverse);
	
	m8 = create_matrix(3,2);
	m8.data[0] = 1.0; m8.data[1] = 2.0;
 	m8.data[2] = 3.0; m8.data[3] = 4.0;
	m8.data[4] = 5.0; m8.data[5] = 6.0;

	m9 = create_matrix(3, 1);
	m9.data[0] = 1.0;
	m9.data[1] = 0.0;
	m9.data[2] = 0.0;
	
	projm = ortho_proj_matrix(&m8);
	printf("\nProjection m8 onto m9:\n");
	print_matrix(&projm);

	identity = create_identity_matrix(5);
	printf("\nIdentity matrix of size 5:\n");
	print_matrix(&identity);

	trans = create_translation_matrix(1, 2, 3);
	printf("\nTranslation matrix 1 2 3:\n");
	print_matrix(&trans);
		

	double sin_check;
	sin_check = mysin(35.0, 0);
	printf("\nsine 35.0 terms 10:\n%f\n", sin_check);

	Matrix scale, rot;
	scale = create_scaling_matrix(4, 5, 6);
	printf("\nScaling matrix 4 5 6:\n");
	print_matrix(&scale);

	rot = create_rotation_matrix(7, 8, 9, 35.0);
	printf("\nRotating Matrix 7 8 9, Theta 35 Degrees:\n");
	print_matrix(&rot);
	
	Vector vec_from_mat;
	vec_from_mat = row_vector(&m8, 2);
	printf("\nRow Vector from Matrix m8:\n");
	print_vector(&vec_from_mat);
	
	Matrix mat_from_vec, transform, p;

	p = homogenous_row_vector(3, 6, 9);
	mat_from_vec = matrix_from_vector(&vec_from_mat);
	printf("\nMatrix from vector, vector_from_mat:\n");
	print_matrix(&mat_from_vec);
	
	transform = transformation(&trans, &scale, &rot, &p);
	printf("\nTransform matrix from trans scale rot and 1 2 3 1 row vector:\n");
	print_matrix(&transform);

	printf("\n-----------\n");
	
	free_matrix(&p);
	free_matrix(&transform);
	free_matrix(&mat_from_vec);
	free_vector(&vec_from_mat);
	free_matrix(&m1);
	free_matrix(&m2);
	free_matrix(&m3);
	free_matrix(&m4);
	free_matrix(&m5);
	free_matrix(&m6);
	free_matrix(&m7);
	free_matrix(&m8);
	free_matrix(&m9);
	free_matrix(&projm);
	free_matrix(&inverse);
	free_matrix(&identity);
	free_matrix(&trans);
	free_matrix(&scale);
	free_matrix(&rot);

	return 0;
}
