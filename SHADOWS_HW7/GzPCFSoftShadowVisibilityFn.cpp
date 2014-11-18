#include "stdafx.h"
#include "rend.h"
#include <math.h> 

float GzPCFSoftShadowVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light) {
	// In world space
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,(*map).camera.Xiw,&screenX,&screenY,&screenZ,&screenW);
	
	if (screenW < 0)
		return 1.0;

	if (screenZ < 0)
		return 1.0;

	screenZ /= screenW;
	
	float distance_from_light_to_obj = sqrt(((*light).position[0]-x)*((*light).position[0]-x) + ((*light).position[1]-y)*((*light).position[1]-y) + ((*light).position[2]-z)*((*light).position[2]-z));

	float d = tan((map->camera.FOV / 2.0) * PI / 180.0);
	float area = ((0.0000001)*screenZ)/(d+screenZ); 
	multiplyMatrixByVector(area, 0, 0,(*map).camera.Xpi,&screenX,&screenY,&screenZ,&screenW);
	float area_in_perspective = screenX;
	multiplyMatrixByVector(area_in_perspective, 0, 0,(*map).Xsp,&screenX,&screenY,&screenZ,&screenW);
	float area_in_screen = screenX;
	if(area_in_screen<=(*map).display->xres){
		
	}
		
	
	
	
	// 1.    Get the z-average
	// 1.1   Compute the area in image space 
	// 
	//1. Get the target points into image space 
	/*
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,map->Ximage_im[3],&screenX,&screenY,&screenZ,&screenW);
	if (screenW < 0)
		return 1.0;

	if (screenZ < 0)
		return 1.0;

	screenZ /= screenW;
	//Propably error : should be in perspective space 
	//2. calculate the size of the z-buff area
	float d = 1/(tan((*map).camera.FOV/2));
	float area_in_image = (*light).size*screenZ/(d+screenZ);
	multiplyMatrixByVector(area_in_image, 0, 0,map->Xwi,&screenX,&screenY,&screenZ,&screenW);
	float area_in_world = screenX;
	multiplyMatrixByVector(area_in_world, 0, 0,map->Ximage_from_world[3],&screenX,&screenY,&screenZ,&screenW);
	float area_in_screen = screenX;
	int t =0;
	
	//Transforming this pixel value to screen space
	/*
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,map->Ximage_from_world[3],&screenX,&screenY,&screenZ,&screenW);
	if (screenW < 0)
		return 1.0;
	screenX /= screenW;
	screenY /= screenW;
	screenZ /= screenW;
	
	// get the x,y,z point in the shadow map
	(*light).size;
	int t =0;
	*/

	return 1.0;
}
