#include <stdio.h>
#include "mathlib.h"
#include <fcntl.h>
#include <linux/fb.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdlib.h>
#include "constants.h"
#include "render.h"
#include <string.h>

int
main() {
	Matrix view, eyePos, targetPos, upVec_temp, forward, right, upVec, perspective, p, p_transposed, t, s, r, trans, transfin;

	targetPos = create_matrix_vector(0.0, 0.0, 0.0, 1.0);
	upVec_temp = create_matrix_vector(0.0, 1.0, 0.0, 0.0);
	eyePos = create_matrix_vector(0.0, 0.0, 5.0, 1.0);
	forward = forward_matrix(&targetPos, &eyePos);
	right = right_matrix(&forward, &upVec_temp);
	upVec = up_matrix(&right, &forward);
	view = create_view_matrix(&right, &upVec, &forward, &eyePos);
	perspective = create_perspective_matrix(FOV, ASPECT, 0.1, 100.0);
	p = create_matrix_vector(6.0, 7.0, 8.0, 1.0);
	p_transposed = transpose(&p);
	t = create_translation_matrix(6.0, 7.0, 8.0);
	s = create_scaling_matrix(2.0, 2.0, 2.0);
	r = create_rotation_matrix(10.0, 1.0, 1.0, 40.0);
       	trans = transformation(&t, &s, &r);	
	transfin = matrix_mult(&trans, &p_transposed);

	free_matrix(&transfin);
	free_matrix(&trans);
	free_matrix(&p_transposed);
	free_matrix(&eyePos);
	free_matrix(&targetPos);
	free_matrix(&upVec);
	free_matrix(&forward);
	free_matrix(&right);
	free_matrix(&view);
	free_matrix(&p);
	free_matrix(&t);
	free_matrix(&s);
	free_matrix(&r);
	free_matrix(&perspective);
	free_matrix(&upVec_temp);

	return 0;

}
