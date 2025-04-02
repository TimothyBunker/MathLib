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


Matrix cubeVertices[NUM_VERTICES];
Matrix screenPoints[NUM_VERTICES];

int edges[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
};

void
transform_vertices(const Matrix *vertices, const Matrix *trans, const Matrix *view, const Matrix *perspective, Matrix *screenPoints, struct fb_var_screeninfo vinfo, size_t num_vertices) {
	size_t i;
	Matrix t1, t2, t3, v_transposed; 
	double x, y, w;
	
	for(i = 0; i < num_vertices; ++i) {
		v_transposed = transpose(&vertices[i]);
		t1 = matrix_mult(trans, &v_transposed);
		t2 = matrix_mult(view, &t1);
		t3 = matrix_mult(perspective, &t2);
		
		w = t3.data[3];
		x = t3.data[0] / w;
		y = t3.data[1] / w;
		
		set_matrix_value(&screenPoints[i], 0, 0, (x + 1) * 0.5 * vinfo.xres);
		set_matrix_value(&screenPoints[i], 1, 0, (1 - y) * 0.5 * vinfo.yres);

		free_matrix(&v_transposed);
		free_matrix(&t1);
		free_matrix(&t2);
		free_matrix(&t3);
	}
}

void
draw_pixel(char *fb_ptr, struct fb_var_screeninfo vinfo, struct fb_fix_screeninfo finfo, double x, double y, int color) {
	long location;
	if (x >= 0 && x < vinfo.xres && y >= 0 && y < vinfo.yres) {
		location = (x + vinfo.xoffset) * (vinfo.bits_per_pixel / 8) + (y + vinfo.yoffset) * finfo.line_length;
		*((int *)(fb_ptr + location)) = color;
	}
}

void
draw_line(char *fb_ptr, struct fb_var_screeninfo vinfo, struct fb_fix_screeninfo finfo, int x0, int y0, int x1, int y1, int color) {
	int dx, dy, err, e2, sx, sy;

	dx = my_abs(x1 - x0);
       	sx = x0 < x1 ? 1 : -1;
	dy = my_abs(y1 - y0);
       	sy = y0 < y1 ? 1 : -1;
	err = dx + dy;

	while(1) {
		draw_pixel(fb_ptr, vinfo, finfo, x0, y0, color);
		if (x0 == x1 && y0 == y1) break;
		e2 = 2 * err;
		if (e2 >= dy) { err += dy; x0 += sx; }
		if (e2 >= dx) { err += dx; x0 += sy; }
	}
}

void
draw_cube(char *fb_ptr, struct fb_var_screeninfo vinfo, struct fb_fix_screeninfo finfo, size_t num_vertices) {
	size_t i;
	int x, y, x0, y0, x1, y1;
	for (i = 0; i < num_vertices; ++i) {
		x = (int)get_matrix_value(&screenPoints[i], 0, 0);
		y = (int)get_matrix_value(&screenPoints[i], 1, 0);
		draw_pixel(fb_ptr, vinfo, finfo, x, y, 0xFFFFFF);
	}

	for (i = 0; i < 12; ++i) {
		x0 = (int)get_matrix_value(&screenPoints[edges[i][0]], 0, 0);
		x1 = (int)get_matrix_value(&screenPoints[edges[i][0]], 1, 0);
		y0 = (int)get_matrix_value(&screenPoints[edges[i][1]], 0, 0);
		y1 = (int)get_matrix_value(&screenPoints[edges[i][1]], 1, 0);
		draw_line(fb_ptr, vinfo, finfo, x0, y0, x1, y1, 0xFFFFFF);
	}
}
