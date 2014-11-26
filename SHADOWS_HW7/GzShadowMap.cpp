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
	(*map)->matlevel = 0;
	(*map)->interp_mode = GZ_NONE;
	GzNewDisplay(&(*map)->display, SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			(*map)->Xsp[i][j] = 0.0;
	(*map)->Xsp[0][0] = ((double)(*map)->display->xres / 2.0);
	(*map)->Xsp[1][1] = (-(double)(*map)->display->yres / 2.0);
	(*map)->Xsp[2][2] = INT_MAX;
	(*map)->Xsp[0][3] = ((double)(*map)->display->xres / 2.0);
	(*map)->Xsp[1][3] = ((double)(*map)->display->yres / 2.0);
	(*map)->Xsp[3][3] = 1.0;

	// FIXME: hardcode below:
	GzCoord look_at_point = {0, 0, 0};
	(*map)->camera.FOV = 83; 
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
	(*map)->matlevel = 0;
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
