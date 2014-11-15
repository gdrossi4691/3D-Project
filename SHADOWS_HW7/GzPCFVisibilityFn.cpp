#include "stdafx.h"
#include "rend.h"
#define	ARRAY(x,y)	(x+(y*display->xres))

float GzPCFVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light) {
	// TODO:
	//Transforming this pixel value to screen space
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,map->Ximage[map->matlevel],&screenX,&screenY,&screenZ,&screenW);

	//Get the Z values for neighbouring pixels using display->fbuf and compare it with transformed pixel Z-values
	float visibility;
	GzDisplay* display = (GzDisplay*)map->display;
	//Checking for screen space bounds
	if(screenX >= 0 && screenX < display->xres && screenY >= 0 && screenY < display->yres){
		GzDepth neighbourZ = display->fbuf[ARRAY((int)screenX,(int)screenY)].z;
		if(screenZ <= neighbourZ){
			visibility = 1.0;
		}
		else{
			visibility = 0.0;
		}
	}
	else if((screenX < 0) || (screenX > display->xres) || (screenY < 0) || (screenY > display->yres)){
		visibility = 1.0;
	}
	return visibility;
}
