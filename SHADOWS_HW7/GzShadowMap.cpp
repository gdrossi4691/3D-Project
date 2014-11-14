#include "stdafx.h"
#include "rend.h"

int GzNewShadowMapCamera(GzRender* map, GzLight* light, GzBoundingBox* bbox) {
	map = new GzRender();
	// TODO:

	return GZ_SUCCESS;
}

int GzNewPerspectiveShadowMapCamera(GzRender* map, GzLight* light, GzBoundingBox bbox){
	map = new GzRender();
	// TODO:

	return GZ_SUCCESS;
}

int GzDeleteShadowMapCamera(GzRender* map) {
	delete map;
	return GZ_SUCCESS;
}
