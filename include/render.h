// Rendering Header
#ifndef RENDER_H
#define RENDER_H

// framebuffer rendering
void draw_pixel(char *fb_ptr, struct fb_var_screeninfo vinfo, struct fb_fix_screeninfo finfo, double x, double y, int color);
void draw_line(char *fb_ptr, struct fb_var_screeninfo vinfo, struct fb_fix_screeninfo finfo, int x0, int y0, int x1, int y1, int color);
void draw_cube(char *fb_ptr, struct fb_var_screeninfo vinfo, struct fb_fix_screeninfo finfo, size_t num_vertices);
void transform_vertices(const Matrix *vertices, const Matrix *trans, const Matrix *view, const Matrix *perspective, Matrix *screenPoints, struct fb_var_screeninfo vinfo, size_t num_vertices);

#define NUM_VERTICES 8
extern Matrix cubeVertices[NUM_VERTICES];
extern Matrix screenPoints[NUM_VERTICES];


#endif // RENDER_H
