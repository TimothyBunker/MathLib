#include <stdio.h>
#include <mathlib.h>
#include <constants.h>

int
main() {
	Matrix eye, target, up, forward, right, updated_up, view, perspective, trans, p, t, s, r; // C, T, U, p is a point
	Vector3 cubeVertices[NUM_VERTICES] = {
		{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}
	};
	Vector2 screenPoints[NUM_VERTICES];

	eye = create_matrix_vector(0.0, 0.0, 5.0, 1.0);
	target = create_matrix_vector(0.0, 0.0, 0.0, 1.0);
	up = create_matrix_vector(0.0, 1.0, 0.0, 0.0);
	up.data[3] = 0;
	
	p = create_matrix_vector(6.0, 7.0, 8.0, 1.0); // new point to go to
	t = create_translation_matrix(3.0, 6.0, 9.0);  // translation matrix
	s = create_scaling_matrix(2.0, 2.0, 2.0);      // scaling matrix
	r = create_rotation_matrix(1.0, 1.0, 1.0, 40); // rotation matrix

	forward = forward_matrix(&target, &eye);
	right = right_matrix(&forward, &up);
	updated_up = up_matrix(&right, &forward);
	view = create_view_matrix(&right, &updated_up, &forward, &eye);
	perspective = create_perspective_matrix(FOV, ASPECT, 0.1, 100.0);
       	trans = transformation(&t, &s, &r);	

	printf("\n-------------\n");
	
	printf("\nEye Matrix:\n");
	print_matrix(&eye);

	printf("\nTarget Matrix:\n");
	print_matrix(&target);

	printf("\nUp Matrix:\n");
	print_matrix(&up);

	printf("\nForward Matrix:\n");
	print_matrix(&forward);

	printf("\nRight Matrix:\n");
	print_matrix(&right);

	printf("\nUpdated Up Matrix:\n");
	print_matrix(&updated_up);

	printf("\nView Matrix:\n");
	print_matrix(&view);

	printf("\nTransformation Matrix:\n");
	print_matrix(&trans);

	printf("\nPerspective Matrix:\n");
	print_matrix(&perspective);
	
	Matrix pview, pclip, pndc, pscreen;

	pview = matrix_mult(&view, &trans);
	printf("\nPview Matrix View x Transformation:\n");
	print_matrix(&pview);

	pclip = matrix_mult(&perspective, &pview);
	printf("\nPClip Matrix Perspective x Pview:\n");
	print_matrix(&pclip);

	pndc = create_pndc(&pclip);
	printf("\nPndc Matrix from Pclip:\n");
	print_matrix(&pndc);

	pscreen = create_pscreen(&pndc, WIDTH, HEIGHT);
	printf("\nPscreen Matrix from pndc:\n");
	print_matrix(&pscreen);

	printf("\n-------------\n");

	free_matrix(&pscreen);
	free_matrix(&pclip);
	free_matrix(&pndc);
	free_matrix(&pview);
	free_matrix(&p);
	free_matrix(&t);
	free_matrix(&s);
	free_matrix(&r);
	free_matrix(&perspective);
	free_matrix(&trans);
	free_matrix(&view);
	free_matrix(&updated_up);
	free_matrix(&right);
	free_matrix(&forward);
	free_matrix(&eye);
	free_matrix(&target);
	free_matrix(&up);	
	return 0;
}
