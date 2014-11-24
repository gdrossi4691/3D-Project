#include "stdafx.h"
#include "rend.h"
#include <math.h> 

float GzPCFSoftShadowVisibilityFn(float world_x, float world_y, float world_z, GzRender* map, GzLight* light, float cos_a) {
	float d = tan((map->camera.FOV / 2.0) * PI / 180.0);
	float w = 1.0;
	GzDisplay* display = (GzDisplay*)map->display;
	///////////////////////////////////
	float image_x, image_y, image_z;
	multiplyMatrixByVector(world_x, world_y, world_z, (*map).camera.Xiw, &image_x, &image_y, &image_z, &w);
	if (w <= 0)
		return 1.0;
	image_x /= w;
	image_y /= w;
	image_z /= w;
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
	float light_size = light->size;
	float area = (light_size * image_z) / (d + image_z); 
	int area_radius_in_screen = (int)(area * 0.5 * display->xres / 2.0 + 0.5); // half of the area !
	int max_x = min(int_screen_x + area_radius_in_screen, display->xres - 1);
	int max_y = min(int_screen_y + area_radius_in_screen, display->yres - 1);
	int min_x = max(int_screen_x - area_radius_in_screen, 0);
	int min_y = max(int_screen_y - area_radius_in_screen, 0);
	float avg = 0.0;
	float counter = 0.0; // counts number of samples in avg

	for (int i = min_y; i <= max_y; i++)
		for (int j = min_x; j <= max_x; j++) {
			int screen_z_from_map = display->fbuf[ARRAY(j, i)].z; 
			if (INT_MAX == screen_z_from_map) // map value is infinity
				continue;
			float z_from_map = screen_z_from_map * d / (INT_MAX - screen_z_from_map); // in image space
			if (z_from_map < image_z) {
				avg += z_from_map;
				counter++;
			}
		}
	if (counter == 0 || avg == 0)
		return GzPCFVisibilityFn(world_x, world_y, world_z, map, light, 1, 1, cos_a);
	
	avg /= counter;
	float filter = light_size * (d + image_z - avg) / avg;	
	int filter_size = (int)(filter * 0.5 * display->xres + 0.5);
	if (filter_size % 2 == 0)
		filter_size++;
	if (filter_size > FILTER_SIZE_LIMIT)
		filter_size = FILTER_SIZE_LIMIT;
	return GzPCFVisibilityFn(world_x, world_y, world_z, map, light, filter_size, filter_size, cos_a);
}
