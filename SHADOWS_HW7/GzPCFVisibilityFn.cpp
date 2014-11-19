#include "stdafx.h"
#include "rend.h"
#define	ARRAY(x,y)	(x+(y*display->xres))
#define	Z_DIFFERENCE_THRESHOLD  100000

int is_in_display_range(GzDisplay* display, int valueX, int valueY) {
	if(valueX >= 0 && valueX < display->xres && valueY >= 0 && valueY < display->yres)
		return 1;
	else 
		return 0;
}


float GzSimpleVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light) {
	//Transforming this pixel value to screen space
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,map->Ximage_from_world[3],&screenX,&screenY,&screenZ,&screenW);
	if (screenW < 0)
		return 1.0;
	screenX /= screenW;
	screenY /= screenW;
	screenZ /= screenW;

	GzDisplay* display = (GzDisplay*)map->display;

	// bilinear interpolation:
	int A_x, A_y, B_x, B_y, C_x, C_y, D_x, D_y;
	A_x = (int)screenX;
	A_y = (int)screenY;
	C_x = (int)screenX + 1;
	C_y = (int)screenY + 1;
	B_x = (int)screenX + 1;
	B_y = (int)screenY;
	D_x = (int)screenX;
	D_y = (int)screenY + 1;
	
	// Color(p) = s t C + (1-s) t D + s (1-t) B + (1-s) (1-t) A 
	float s = screenX - A_x;
	float t = screenY - A_y;
	
	//Checking for screen space bounds
	if(is_in_display_range(display, A_x, A_y) && is_in_display_range(display, B_x, B_y) && is_in_display_range(display, C_x, C_y) && is_in_display_range(display, D_x, D_y)){
		int interpolated_z_from_map	= s * t * display->fbuf[ARRAY(C_x, C_y)].z + (1 - s) * t * display->fbuf[ARRAY(D_x, D_y)].z +
			s * (1 - t) *  display->fbuf[ARRAY(B_x, B_y)].z + (1 - s) * (1 - t) *  display->fbuf[ARRAY(A_x, A_y)].z;
		int screen_z = (int)(screenZ + 0.5);
		int diff = screen_z - interpolated_z_from_map;
		if(diff <= Z_DIFFERENCE_THRESHOLD)
			return 1.0;
		else
			return 0.0;
	}
	int scr_x = (int)(screenX + 0.5); 
	int scr_y = (int)(screenY + 0.5);
	if(is_in_display_range(display, scr_x, scr_y)){
		int interpolated_z_from_map	= display->fbuf[ARRAY(scr_x, scr_y)].z;
		int screen_z = (int)(screenZ + 0.5);
		int diff = screen_z - interpolated_z_from_map;
		if(diff <= Z_DIFFERENCE_THRESHOLD)
			return 1.0;
		else
			return 0.0;
	} 

	return 1.0;
}

float GzPCFVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light) {
	//Transforming this pixel value to screen space
	float screenX, screenY, screenZ, screenW;
	multiplyMatrixByVector(x, y, z,map->Ximage_from_world[3],&screenX,&screenY,&screenZ,&screenW);
	screenX /= screenW;
	screenY /= screenW;
	screenZ /= screenW;

	float visibility = 1.0;
	float filter[5][5];
	if (screenW < 0)
		return 1.0;
	//Get the Z values for neighbouring pixels using display->fbuf and compare it with transformed pixel Z-values
	GzDisplay* display = (GzDisplay*)map->display;
	int fbX = (int)(screenX + 0.5);
	int fbY = (int)(screenY + 0.5);
	int fbZ = (int)(screenZ + 0.5);
	//Checking for screen space bounds
	for(int i = fbX-2;i <= fbX+2; i++){
		for(int j = fbY-2; j <= fbY+2; j++){
			if((i<0 ||i>=display->xres) ||(j<0 || j>=display->yres)){
				return 1.0;
			}
				if(fbX >= 0 && fbX < display->xres && fbY >= 0 && fbY < display->yres){
					GzDepth neighbourZ = display->fbuf[ARRAY(i,j)].z;
					float diff = fbZ - (float)neighbourZ;
					if(diff <= Z_DIFFERENCE_THRESHOLD){
						filter[i-(fbX-2)][j-(fbY-2)] = 1.0;
					}
					else{
						filter[i-(fbX-2)][j-(fbY-2)] = 0.0;
					}
				}
			else {
				filter[i-(fbX-2)][j-(fbY-2)] = 1.0;
			}
		}
	}
	filter[2][2] = visibility;
	visibility = (filter[0][0] + filter[0][1] + filter[0][2]+ filter[0][3]+ filter[0][4] + filter[1][0]+ filter[1][1] + filter[1][2]+ filter[1][3]+ filter[1][4]+ filter[2][0]+ filter[2][1]+ filter[2][2]+ filter[2][3]+ filter[2][4]+ filter[3][0]+ filter[3][1]+ filter[3][2]+ filter[3][3]+ filter[3][4]+ filter[4][0]+ filter[4][1]+ filter[4][2]+ filter[4][3]+ filter[4][4])/25;
	return visibility;
}
