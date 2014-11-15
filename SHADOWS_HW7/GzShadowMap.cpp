#include "stdafx.h"
#include "rend.h"
#include "math.h"

// light position and box positions are in model space...
int GzNewShadowMapCamera(GzRender** map, GzLight* light, GzBoundingBox* bbox) {


	// get all the things in camera initialized
	*map = new GzRender();
	(*map)->numlights = 0;
	(*map)->spec = 0;
	(*map)->tex_fun = NULL;
	(*map)->shift_x = 0.0;
	(*map)->shift_y = 0.0;

	// temprorary stub:
	(*map)->matlevel = 0;
	(*map)->interp_mode = GZ_FLAT;
	GzNewDisplay(&(*map)->display, 512, 512);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			(*map)->Xsp[i][j] = 0.0;
	(*map)->Xsp[0][0] = ((double)(*map)->display->xres / 2.0);
	(*map)->Xsp[1][1] = (-(double)(*map)->display->yres / 2.0);
	(*map)->Xsp[2][2] = INT_MAX;
	(*map)->Xsp[0][3] = ((double)(*map)->display->xres / 2.0);
	(*map)->Xsp[1][3] = ((double)(*map)->display->yres / 2.0);
	(*map)->Xsp[3][3] = 1.0;

	// TODO:
	// calculate the FOV - method 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// to be continue if current code has error... do it properly i guess...

	// calculate the FOV - method 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// find the longest path on the plane surface
	
	float diagono_distance_on_the_plane = sqrt((bbox->Xmax- bbox->Xmin)*(bbox->Xmax- bbox->Xmin) + (bbox->Ymax- bbox->Ymax)*(bbox->Ymax- bbox->Ymin) + (bbox->Ymax- bbox->Zmax)-(bbox->Ymax- bbox->Zmin));
	 
	// find the shortest path from the camera to the point
	float x_short1 =  abs(light->position[0]-bbox->Xmax);
	float x_short2 =  abs(light->position[0]-bbox->Xmin);

	float y_short1 =  abs(light->position[1]-bbox->Ymax);
	float y_short2 =  abs(light->position[1]-bbox->Ymin);

	float z_short1 =  abs(light->position[2]-bbox->Zmax);
	float z_short2 =  abs(light->position[2]-bbox->Zmin);

	float x_MIN = x_short1;
	if(x_short1>x_short2){
		x_MIN= x_short2;
	}

	float y_MIN = y_short1;
	if(y_short1>y_short2){
		y_MIN= y_short2;
	}

	float z_MIN = z_short1;
	if(z_short1>z_short2){
		z_MIN= z_short2;
	}


	// shortest_distance_from_camera_to_plane can be replace by distance_from_camera_to_center_cube
	float shortest_distance_from_camera_to_plane = sqrt((x_MIN*x_MIN)+ (y_MIN*y_MIN) + (z_MIN*z_MIN));

	GzCoord look_at_point = {(bbox->Xmax-bbox->Xmin)/2, (bbox->Ymax-bbox->Ymin)/2, (bbox->Zmax-bbox->Zmin)/2};
	float distance_from_camera_to_center_cube = sqrt(  ((*light).position[0]-look_at_point[0])*((*light).position[0]-look_at_point[0])+ (((*light).position[1]-look_at_point[1])*((*light).position[1]-look_at_point[1])) + ((*light).position[2]-look_at_point[2])*  ((*light).position[2]-look_at_point[2]));

	float FOV_in_radius = (diagono_distance_on_the_plane/distance_from_camera_to_center_cube);


	(*map)->camera.FOV =  FOV_in_radius;
	(*map)->camera.position[0] = (*light).position[0];
	(*map)->camera.position[1] = (*light).position[1];
	(*map)->camera.position[2] = (*light).position[2];
	(*map)->camera.lookat[0] = look_at_point[0];
	(*map)->camera.lookat[1] = look_at_point[1];
	(*map)->camera.lookat[2] = look_at_point[2];
	(*map)->camera.worldup[0] = 0;
	(*map)->camera.worldup[1] = 1;
	(*map)->camera.worldup[2] = 0;
	GzBeginRender(*map);
	return GZ_SUCCESS;
}

int GzNewPerspectiveShadowMapCamera(GzRender** map, GzLight* light, GzBoundingBox bbox){
	*map = new GzRender();
	(*map)->numlights = 0;
	(*map)->spec = 0;
	(*map)->tex_fun = NULL;
	(*map)->shift_x = 0.0;
	(*map)->shift_y = 0.0;
	// temprorary stub:
	(*map)->matlevel = 1;
	(*map)->interp_mode = GZ_FLAT;
	GzNewDisplay(&(*map)->display, 512, 512);
	// TODO:

	return GZ_SUCCESS;
}

int GzDeleteShadowMapCamera(GzRender* map) {
	GzFreeDisplay(map->display);
	delete map;
	return GZ_SUCCESS;
}
