#include "stdafx.h"
#include "rend.h"

int GzNewShadowMapCamera(GzRender** map, GzLight* light, GzBoundingBox* bbox) {
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
