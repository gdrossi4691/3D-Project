#include "stdafx.h"
#include "rend.h"
#include <math.h> 
#define	ARRAY(x,y)	(x+(y*display->xres))

int is_in_display_range(GzDisplay* display, int valueX, int valueY) {
	if(valueX >= 0 && valueX < display->xres && valueY >= 0 && valueY < display->yres)
		return 1;
	else 
		return 0;
}

float GzPCFVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light, int filter_size_x, int filter_size_y) {
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
	screen_y /=w; 
	int int_screen_y = (int)(screen_y + 0.5);
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
		if (diff <= Z_DIFFERENCE_THRESHOLD)
			return 1.0;
		else 
			return 0.0;
	}
	///if orbitrary filter (hard shadows)///
	int radius_x = (filter_size_x - 1)/2;
	int radius_y = (filter_size_y - 1)/2;
	int max_x = min(int_screen_x + radius_x, display->xres - 1);
	int min_x = max(int_screen_x - radius_x, 0);
	int max_y = min(int_screen_y + radius_y, display->yres - 1);
	int min_y = max(int_screen_y - radius_y, 0);
	
	float visibility = 0.0;
	float counter = 0.0;

	for (int i = min_y; i <= max_y; i++)
		for (int j = min_x; j <= max_x; j++) {
			counter++;
			int screen_z_from_map = display->fbuf[ARRAY(j, i)].z; 
			if (INT_MAX == screen_z_from_map){ // map value is infinity
				visibility += 1.0;
				continue;
			}
			z_from_map = screen_z_from_map * d / (INT_MAX - screen_z_from_map); // in image space
			float diff = image_z - z_from_map; // compare z in image space 
			if (diff <= Z_DIFFERENCE_THRESHOLD) // visible
				visibility += 1.0; // else add 0.0
		}
	visibility = counter != 0 ? visibility / counter : 1;
	return visibility;
}

float GzSimpleVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light) {
	return GzPCFVisibilityFn(world_x, world_y, world_z, map, light, 1, 1);
}

float GzPCFVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light) {
	return GzPCFVisibilityFn(world_x, world_y, world_z, map, light, FILTER_SIZE_X, FILTER_SIZE_Y);
}
