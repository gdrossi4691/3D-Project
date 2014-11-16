#include "stdafx.h"
#include "rend.h"
#define	ARRAY(x,y)	(x+(y*display->xres))

float GzPCFVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light) {
	//Transforming this pixel value to screen space
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,map->Ximage[map->matlevel],&screenX,&screenY,&screenZ,&screenW);
	screenX /= screenW;
	screenY /= screenW;
	screenZ /= screenW;

	float visibility = 1.0;
	if (screenW < 0)
		return 1.0;
	//Get the Z values for neighbouring pixels using display->fbuf and compare it with transformed pixel Z-values
	GzDisplay* display = (GzDisplay*)map->display;
	int fbX = (int)(screenX + 0.5);
	int fbY = (int)(screenY + 0.5);
	int fbZ = (int)(screenZ + 0.5);
	//Checking for screen space bounds
	if(fbX >= 0 && fbX < display->xres && fbY >= 0 && fbY < display->yres){
		GzDepth neighbourZ = display->fbuf[ARRAY(fbX,fbY)].z;
		float diff = screenZ - (float)neighbourZ;
		if(diff <= 1000000){
			return 1.0;
		}
		else{
			return 0.0;
		}
	}
	else {
		return 1.0;
	}
	return visibility;
}
