#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathlib.h"

Vector
create_vector(size_t size) {
	Vector v;
	v.size = size;
	v.data = (double *)malloc(size * sizeof(double));
	return v;
}

void
free_vector(Vector *v){
	free(v->data);
	v->data = NULL;
	v->size = 0;
}

void
print_vector(const Vector *v) {
	size_t i;
	for(i = 0; i < v->size; ++i) {
		printf("%f ", v->data[i]);
	}
	printf("\n");
}

Vector
mult_vector(double n, const Vector *v) {
	size_t i;
	Vector result;

	result = create_vector(v->size);

	for (i = 0; i < v->size; ++i) {
		result.data[i] = v->data[i]*n; 
	}
	return result;
}

Vector 
vector_add(const Vector *v1, const Vector *v2) {
	Vector result;
	size_t i;
	if (v1->size != v2->size) {
		fprintf(stderr, "Error: Vector sizes do not match for addition.\n");
		exit(EXIT_FAILURE);
	}	
	
	result = create_vector(v1->size);
	for(i = 0; i < v1->size; ++i) {
		result.data[i] = v1->data[i] + v2->data[i];
	}
	return result;
}

Vector 
vector_sub(const Vector *v1, const Vector *v2) {
	Vector result;
	size_t i;
	if (v1->size != v2->size) {
		fprintf(stderr, "Error: Vector sizes do not match for addition.\n");
		exit(EXIT_FAILURE);
	}	
	
	result = create_vector(v1->size);
	for(i = 0; i < v1->size; ++i) {
		result.data[i] = v1->data[i] - v2->data[i];
	}
	return result;
}

Vector
zero_vector(const Vector *v) {
	size_t i;
	Vector result;

	result = create_vector(v->size);
	for (i= 0; i < v->size; ++i) {
		result.data[i] = 0.0;
	}
	return result;
}

double
dot(const Vector *v1, const Vector *v2) {
	size_t i;
	double result;

	
	if(v1->size != v2->size) {
		fprintf(stderr, "Error: Vectors do not match size for dot product.\n");
		exit(EXIT_FAILURE);
	}	

	result = 0.0;	
	for (i = 0; i < v1->size; ++i) {
		result = result + v1->data[i] * v2->data[i];
	}
	return result;
}

Vector
proj_vector(const Vector *v1, const Vector *v2) {
	double num, den;
	Vector result;

	num = dot(v1, v2);
	den = dot(v2, v2);
	
	result = mult_vector(num/den, v2);

	return result;
}

Vector
vector_from_coords(double x, double y, double z) {
	Vector result;

	result = create_vector(3.0);
	result.data[0] = x;
	result.data[1] = y;
	result.data[2] = z;

	return result; 
}

double
magnitude(const Vector *v){
	double result;
	size_t i;

	for (i = 0; i < v->size; ++i) {
		result += v->data[i] * v->data[i];
	}

	return sqrt_approx(result); 
}

Vector 
normalize(const Vector *v) {
	Vector result;
	size_t i;
	double m;

	result = create_vector(v->size);
	m = magnitude(v); 
	
	for (i = 0; i < v->size; ++i) {
		result.data[i] = v->data[i]/m;
	}

	return result;
}


Vector
cross(const Vector *v1, const Vector *v2) {
	Vector result;
	result = create_vector(3.0);

	result.data[0] = v1->data[1] * v2->data[2] - v1->data[2] * v2->data[1];
	result.data[1] = v1->data[2] * v2->data[0] - v1->data[0] * v2->data[2];
       	result.data[2] = v1->data[0] * v2->data[1] - v1->data[1] * v2->data[0];
	
	if (result.data[0] == -0.0f) 
		result.data[0] = 0.0f;
	if (result.data[1] == -0.0f) 
		result.data[1] = 0.0f;
	if (result.data[2] == -0.0f)
		result.data[2] = 0.0f;

	return result; 
}

Vector
row_vector(const Matrix *m, int row) {
	Vector result;
	size_t i;
	
	result = create_vector(m->cols);
	for (i = 0; i < m->cols; ++i) {
		result.data[i] = m->data[row * m->cols + i];  
	}

	return result;
}

