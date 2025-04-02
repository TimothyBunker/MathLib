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
	long screensize;
	char *fb_ptr;
	// framebuffer init
	int fb_fd = open("/dev/fb0", O_RDWR);
	if (fb_fd == -1) {
		perror("error opening framebuffer device");
		return -1;
	}

	struct fb_var_screeninfo vinfo;
	struct fb_fix_screeninfo finfo;

	if (ioctl(fb_fd, FBIOGET_FSCREENINFO, &finfo)) {
		perror("Error reading fixed information");
		return -1;
	}

	if (ioctl(fb_fd, FBIOGET_VSCREENINFO, &vinfo)) {
		perror("Error reading variable information");
		return -1;
	}

	screensize = vinfo.yres_virtual * finfo.line_length;
	fb_ptr = (char *)mmap(0, screensize, PROT_READ | PROT_WRITE, MAP_SHARED, fb_fd, 0);

	if ((intptr_t)fb_ptr == -1) {
		perror("Error mapping framebuffer device to memory");
		return -1;
	}

	// painstakingly manually setup cube vertices like a boss
	cubeVertices[0] = create_matrix_vector(-1, -1, -1, 1);
	cubeVertices[1] = create_matrix_vector(1, -1, -1, 1);
	cubeVertices[2] = create_matrix_vector(1, 1, -1, 1);
	cubeVertices[3] = create_matrix_vector(-1, 1, 1, 1);
	cubeVertices[4] = create_matrix_vector(-1, -1, 1, 1);
	cubeVertices[5] = create_matrix_vector(1, -1, 1, 1);
	cubeVertices[6] = create_matrix_vector(1, 1, 1, 1);
	cubeVertices[7] = create_matrix_vector(-1, 1, 1, 1);

	for (size_t i = 0; i < NUM_VERTICES; ++i) {
		init_matrix(&screenPoints[i], NUM_VERTICES, 2);
	}
	Matrix view, eyePos, targetPos, upVec_temp, forward, right, upVec, perspective, p, p_transposed, t, s, r, trans;

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

	while(1) {
		// clear framebuffer (optional apparently lol)
		memset(fb_ptr, 0, screensize);
		transform_vertices(cubeVertices, &trans, &view, &perspective, screenPoints, vinfo, NUM_VERTICES);

		draw_cube(fb_ptr, vinfo, finfo, NUM_VERTICES);
		
		// delay
		usleep(32000); // about 60FPS
		

	}

	// clean
	munmap(fb_ptr, screensize);
	close(fb_fd);

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

	for (size_t i = 0; i < NUM_VERTICES; i++) {
		free_matrix(&screenPoints[i]);
		free_matrix(&cubeVertices[i]);
	}
	return 0;

}
