#include "stdafx.h"
#include "rend.h"
#include <math.h> 
#define	ARRAY(x,y)	(x+(y*display->xres))

float multipy_vectors(float x1, float y1, float z1, float x2, float y2, float z2) {
	return x1 * x2 + y1 * y2 + z1 * z2;
}

void normalize(float* x,  float* y) {
	float len = sqrt(*x * *x + *y * *y);
	*x /= len;
	*y /= len;
	return;
}

void Triangle::compute_cos(float* cos_ax, float* cos_ay, float* sin_ax, float* sin_ay, GzRender* map, GzLight* light){
	float w = 1.0;
	float image_x1, image_y1, image_z1, image_x2, image_y2, image_z2, image_x3, image_y3, image_z3;
	multiplyMatrixByVector(world_x1, world_y1, world_z1, map->camera.Xiw, &image_x1, &image_y1, &image_z1, &w);
	image_x1 /=w;
	image_y1 /=w;
	image_z1 /=w;
	multiplyMatrixByVector(world_x2, world_y2, world_z2, map->camera.Xiw, &image_x2, &image_y2, &image_z2, &w);
	image_x2 /=w;
	image_y2 /=w;
	image_z2 /=w;
	multiplyMatrixByVector(world_x3, world_y3, world_z3, map->camera.Xiw, &image_x3, &image_y3, &image_z3, &w);
	image_x3 /=w;
	image_y3 /=w;
	image_z3 /=w;

	Plane pln = init_plane(image_x1, image_y1, image_z1, image_x2, image_y2, image_z2, image_x3, image_y3, image_z3);
	GzCoord plane_norm_img;
	plane_norm_img[0] = pln.A;
	plane_norm_img[1] = pln.B;
	plane_norm_img[2] = pln.C;
	GzNormalize_vector(plane_norm_img + 0, plane_norm_img + 1, plane_norm_img + 2);

	if (-plane_norm_img[2] < 0.0) {
		plane_norm_img[0] *= -1.0;
		plane_norm_img[1] *= -1.0;
		plane_norm_img[2] *= -1.0;
	}

	float n_x, n_y, n_z;
	n_x = plane_norm_img[0];
	n_y = plane_norm_img[1];
	n_z = plane_norm_img[2];
	
	float sqrt_nx2_and_nz2 = sqrt(n_x * n_x + n_z * n_z);
	*cos_ay = -n_z / sqrt_nx2_and_nz2;
	*sin_ay = n_x / sqrt_nx2_and_nz2;
	*cos_ax = sqrt_nx2_and_nz2;
	*sin_ax = -n_y;
}

int rotX(float cos_ax, float sin_ax, float* x, float* y, float* z) {
	GzMatrix m = { 
		1.0,	0.0,	0.0,	0.0, 
		0.0,	cos_ax,	-sin_ax,	0.0, 
		0.0,	sin_ax,	cos_ax,	0.0, 
		0.0,	0.0,	0.0,	1.0 
	};
	float old_x = *x;
	float old_y = *y;
	float old_z = *z;
	float w = 1.0;
	multiplyMatrixByVector(old_x, old_y, old_z, m, x, y, z, &w);
	if (w != 0){
		*x /= w;
		*y /= w;
		*z /= w;
		return GZ_SUCCESS;
	}  else 
		return GZ_FAILURE;
}

int rotY(float cos_ay, float sin_ay, float* x, float* y, float* z) {
	GzMatrix m = { 
		cos_ay,	0.0,	sin_ay,	0.0, 
		0.0,	1.0,	0.0,	0.0, 
		-sin_ay,	0.0,	cos_ay,	0.0, 
		0.0,	0.0,	0.0,	1.0 
	};
	float old_x = *x;
	float old_y = *y;
	float old_z = *z;
	float w = 1.0;
	multiplyMatrixByVector(old_x, old_y, old_z, m, x, y, z, &w);
	if (w != 0){
		*x /= w;
		*y /= w;
		*z /= w;
		return GZ_SUCCESS;
	}  else 
		return GZ_FAILURE;
}

int is_in_display_range(GzDisplay* display, int valueX, int valueY) {
	if(valueX >= 0 && valueX < display->xres && valueY >= 0 && valueY < display->yres)
		return 1;
	else 
		return 0;
}

float distance_sqr(float x1, float y1, float x2, float y2) {
	return (x2- x1) * (x2- x1)  + (y2 - y1) * (y2 - y1);
}

float invers_sqr_dist_weight_fn(float x1, float y1, float x2, float y2) {
	return distance_sqr(x1, y1, x2, y2);
}
float sqr_dist_weight_fn(float x1, float y1, float x2, float y2, float radius_sqr) {
	return radius_sqr - distance_sqr(x1, y1, x2, y2);
}

float Triangle::GzPCFVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light, int filter_size_x, int filter_size_y) {
	float d = 1.0 / tan((map->camera.FOV / 2.0) * (PI / 180.0));
	float w = 1.0;
	GzDisplay* display = (GzDisplay*)map->display;
	///////////////////////////////////
	float image_x, image_y, image_z; // point coordinates in image space
	multiplyMatrixByVector(world_x, world_y, world_z,map->camera.Xiw, &image_x, &image_y, &image_z, &w);
	if (w <= 0)
		return 1.0;
	image_x /=w;
	image_y /=w;
	image_z /=w;
	///////////////////////////////////
	float screen_x, screen_y, screen_z; // point coordinates in screen space
	multiplyMatrixByVector(world_x, world_y, world_z, map->Ximage_from_world[3], &screen_x, &screen_y, &screen_z, &w);
	if (w <= 0)
		return 1.0;
	screen_x /=w; 
	int int_screen_x = (int)(screen_x + 0.5);
	float displacement_x = screen_x - int_screen_x;
	screen_y /=w; 
	int int_screen_y = (int)(screen_y + 0.5);
	float displacement_y = screen_y - int_screen_y;
	screen_z /=w;

	///////////////////////////////////
	if(!is_in_display_range(display, int_screen_x, int_screen_y) ||
		(INT_MAX == (display->fbuf[ARRAY(int_screen_x, int_screen_y)].z))) // point out of shadow map or map value is infinity
		return 1.0;
	float z_from_map; // in image space
	///if 1x1 filter (hard shadows)////
	if (filter_size_x == 1 && filter_size_y == 1) {
		// bilinear interpolation: = s t C + (1-s) t D + s (1-t) B + (1-s) (1-t) A 
		int A_x, A_y, B_x, B_y, C_x, C_y, D_x, D_y;
		A_x = (int)screen_x;
		A_y = (int)screen_y;
		C_x = (int)screen_x + 1;
		C_y = (int)screen_y + 1;
		B_x = (int)screen_x + 1;
		B_y = (int)screen_y;
		D_x = (int)screen_x;
		D_y = (int)screen_y + 1;
		// check screen bounds and finity of the surrounding pixels
		if (!is_in_display_range(display, A_x, A_y) || 
			!is_in_display_range(display, B_x, B_y) || 
			!is_in_display_range(display, C_x, C_y) || 
			!is_in_display_range(display, D_x, D_y) ||
			INT_MAX == (display->fbuf[ARRAY(A_x, A_y)].z) ||
			INT_MAX == (display->fbuf[ARRAY(B_x, B_y)].z) ||
			INT_MAX == (display->fbuf[ARRAY(C_x, C_y)].z) ||
			INT_MAX == (display->fbuf[ARRAY(D_x, D_y)].z) ){ // use only closest point (which is visible)
			int screen_z_from_map = display->fbuf[ARRAY(int_screen_x, int_screen_y)].z; 
			z_from_map = screen_z_from_map * d / (INT_MAX - screen_z_from_map);
		}
		else {	// else use bilinear interpolation in image space:
			float Ascreen_z = display->fbuf[ARRAY(A_x, A_y)].z;
			float Bscreen_z = display->fbuf[ARRAY(B_x, B_y)].z;
			float Cscreen_z = display->fbuf[ARRAY(C_x, C_y)].z;
			float Dscreen_z = display->fbuf[ARRAY(D_x, D_y)].z;
			float s = screen_x - A_x;
			float t = screen_y - A_y;
			float screen_z_from_map = s * t * Cscreen_z + (1.0 - s) * t * Dscreen_z + s * (1.0 - t) * Bscreen_z + (1.0 - s) * (1.0 - t) * Ascreen_z;
			z_from_map = screen_z_from_map * d / (INT_MAX - screen_z_from_map);
		}
		float diff = image_z - z_from_map; // compare z in image space 
		if (diff <= CONSTANT_BAIS)
			return 1.0;
		else 
			return 0.0;
	}
	///if orbitrary filter///
	int radius_x = (filter_size_x - 1)/2;
	int radius_y = (filter_size_y - 1)/2;
	int radius_sqr = (min(radius_x, radius_y) + 1) * (min(radius_x, radius_y) + 1);
	int max_x = min(int_screen_x + radius_x, display->xres - 1);
	int min_x = max(int_screen_x - radius_x, 0);
	int max_y = min(int_screen_y + radius_y, display->yres - 1);
	int min_y = max(int_screen_y - radius_y, 0);
	
	float visibility = 0.0;
	float norm = 0.0;
	// partial coverage
	max_x = displacement_x > 0 ? max_x + 1 : max_x;
	min_x = displacement_x < 0 ? min_x - 1 : min_x;
	max_y = displacement_y > 0 ? max_y + 1 : max_y; 
	min_y = displacement_y < 0 ? min_y - 1 : min_y;
	// screen boundaries
	max_x = max_x > (display->xres - 1) ? (display->xres - 1) : max_x;
	min_x = min_x < 0 ? 0 : min_x;
	max_y = max_y > (display->yres - 1) ? (display->yres - 1) : max_y; 
	min_y = min_y < 0 ? 0 : min_y;
	int no_clip_xr = max_x == (display->xres - 1) ?	 0 : 1;
	int no_clip_xl = min_x == 0 ? 0 : 1;
	int no_clip_yr = max_y == (display->yres - 1) ?	 0 : 1;
	int no_clip_yl = min_y == 0 ? 0 : 1;

	float image_x_rotated = image_x;
	float image_y_rotated = image_y;
	float image_z_rotated = image_z;

	float cos_ax, cos_ay, sin_ax, sin_ay;
	compute_cos(&cos_ax, &cos_ay, &sin_ax, &sin_ay, map, light);
	rotY(cos_ay, sin_ay, &image_x_rotated, &image_y_rotated, &image_z_rotated);
	rotX(cos_ax, sin_ax, &image_x_rotated, &image_y_rotated, &image_z_rotated);

	for (int i = min_y; i <= max_y; i++){
		float y_coaf = 1.0;
		if (i == 0 && no_clip_yl)
			y_coaf = displacement_y > 0 ? (1.0 - displacement_y) : displacement_y;
		if (i == max_y && no_clip_yr)
			y_coaf = displacement_y > 0 ? displacement_y : (1.0 - displacement_y);

		for (int j = min_x; j <= max_x; j++) {
			float x_coaf = 1.0;
			if (j == 0 && no_clip_xl)
				x_coaf = displacement_x > 0 ? (1.0 - displacement_x) : displacement_x;
			if (j  == max_x && no_clip_xr)
				x_coaf = displacement_x > 0 ? displacement_x : (1.0 - displacement_x);

			float weight = 1.0;
			weight = sqr_dist_weight_fn(j, i, screen_x, screen_y, radius_sqr);
			int screen_z_from_map = display->fbuf[ARRAY(j, i)].z; 
			if (INT_MAX == screen_z_from_map){ // map value is infinity
				visibility += x_coaf * y_coaf * weight;
				norm += x_coaf * y_coaf * weight;
				continue;
			}
			z_from_map = screen_z_from_map * d / (INT_MAX - screen_z_from_map); // in image space
			float x_from_map = (2.0 * (j - map->display->x_shift) / map->display->xres - 1.0) * (((float)screen_z_from_map) / (INT_MAX - (float)screen_z_from_map) + 1.0); // in image space
			float y_from_map = (1.0 - 2.0 * (i - map->display->y_shift) / map->display->yres) * (((float)screen_z_from_map) / (INT_MAX - (float)screen_z_from_map) + 1.0); // in image space
			rotY(cos_ay, sin_ay, &x_from_map, &y_from_map, &z_from_map);
			rotX(cos_ax, sin_ax, &x_from_map, &y_from_map, &z_from_map);

			float diff = image_z_rotated - z_from_map; // compare z in image space 
			if (diff <= CONSTANT_BAIS) // visible
				visibility += x_coaf * y_coaf * weight; // else add 0.0
			norm += x_coaf * y_coaf * weight;
		}
	}
	visibility = norm != 0.0 ? (visibility / norm) : 1;
	visibility = visibility > 1.0 ? 1.0 : visibility;
	return visibility;
}

float Triangle::GzSimpleVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light) {
	return GzPCFVisibilityFn(world_x, world_y, world_z, map, light, 1, 1);
}

float Triangle::GzPCFVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light) {
	return GzPCFVisibilityFn(world_x, world_y, world_z, map, light, FILTER_SIZE_X, FILTER_SIZE_Y);
}
