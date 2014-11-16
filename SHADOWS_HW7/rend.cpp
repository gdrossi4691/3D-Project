/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"


/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}


int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	//degree = -degree;
	GzMatrix m = { 
		1.0,	0.0,						0.0,						0.0, 
		0.0,	cos(degree * PI / 180.0 ),	-sin(degree * PI / 180.0 ),	0.0, 
		0.0,	sin(degree * PI / 180.0 ),	cos(degree * PI / 180.0 ),	0.0, 
		0.0,	0.0,						0.0,						1.0 
	};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat[i][j] = m[i][j];
	
	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	//degree = -degree;
	GzMatrix m = { 
		cos(degree * PI / 180.0 ),	0.0,	sin(degree * PI / 180.0 ),	0.0, 
		0.0,						1.0,	0.0,						0.0, 
		-sin(degree * PI / 180.0 ),	0.0,	cos(degree * PI / 180.0 ),	0.0, 
		0.0,						0.0,	0.0,						1.0 
	};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat[i][j] = m[i][j];
	
	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	//degree = -degree;
	GzMatrix m = { 
		cos(degree * PI / 180.0 ),	-sin(degree * PI / 180.0 ),	0.0,	0.0, 
		sin(degree * PI / 180.0 ),	cos(degree * PI / 180.0 ),	0.0,	0.0, 
		0.0,						0.0,						1.0,	0.0, 
		0.0,						0.0,						0.0,	1.0 
	};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat[i][j] = m[i][j];
	
	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	GzMatrix m = { 
		1.0,	0.0,	0.0,	translate[0], 
		0.0,	1.0,	0.0,	translate[1], 
		0.0,	0.0,	1.0,	translate[2], 
		0.0,	0.0,	0.0,	1.0 
	};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat[i][j] = m[i][j];
	
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	GzMatrix m = { 
		scale[0],	0.0,		0.0,		0.0, 
		0.0,		scale[1],	0.0,		0.0, 
		0.0,		0.0,		scale[2],	0.0, 
		0.0,		0.0,		0.0,		1.0 
	};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat[i][j] = m[i][j];
	
	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions


// positioin, worldup, lookup are input parameters.
void GzInitCameraXiw(GzRender *render){
	// Z from position and lookat
	GzCoord z;
	z[0] = (render->camera.lookat[0] - render->camera.position[0]);
	z[1] = (render->camera.lookat[1] - render->camera.position[1]);
	z[2] = (render->camera.lookat[2] - render->camera.position[2]);
	
	float normZ = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
	z[0] = z[0] / normZ;
	z[1] = z[1] / normZ;
	z[2] = z[2] / normZ;

	// Y from worldup and z
	float up_dot_Z = render->camera.worldup[0] * z[0] + render->camera.worldup[1] * z[1] + render->camera.worldup[2] * z[2];
	GzCoord y;
	y[0] =  render->camera.worldup[0] - up_dot_Z * z[0];
	y[1] =  render->camera.worldup[1] - up_dot_Z * z[1];
	y[2] =  render->camera.worldup[2] - up_dot_Z * z[2];

	float normY = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	y[0] = y[0] / normY;
	y[1] = y[1] / normY;
	y[2] = y[2] / normY;

	// X from z and y
	// i (Yy Zz - Yz Zy) + 
	// j (Yz Zx - Yx Zz) + 
	// k (Yx Zy - Yy Zx) 
	GzCoord x;
	x[0] = y[1] * z[2] - y[2] * z[1];
	x[1] = y[2] * z[0] - y[0] * z[2];
	x[2] = y[0] * z[1] - y[1] * z[0];

	// Xiw
	GzMatrix	iw = { 
		x[0],	x[1],	x[2],	-(x[0] * render->camera.position[0] + x[1] * render->camera.position[1] + x[2] * render->camera.position[2]), 
		y[0],	y[1],	y[2],	-(y[0] * render->camera.position[0] + y[1] * render->camera.position[1] + y[2] * render->camera.position[2]), 
		z[0],	z[1],	z[2],	-(z[0] * render->camera.position[0] + z[1] * render->camera.position[1] + z[2] * render->camera.position[2]), 
		0.0,	0.0,	0.0,	1.0 
	}; 

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			render->camera.Xiw[i][j] = iw[i][j];
	// Xwi
	GzMatrix	wi = { 
		x[0],	y[0],	z[0], render->camera.position[0],
		x[1],	y[1],	z[1], render->camera.position[1],
		x[2],	y[2],	z[2], render->camera.position[2],
		0.0,	0.0,	0.0,	1.0 
	}; 

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			render->Xwi[i][j] = wi[i][j];

	return;
}


int GzUpdateRender(GzRender *render, GzDisplay	*display)
{
	render->display = display;
	return GZ_SUCCESS;
}

int GzNewRender(GzRender **render, GzDisplay	*display)
{
/*  
- malloc a renderer struct 
- setup Xsp and anything only done once 
- save the pointer to display 
- init default camera 
*/ 
	// malloc a renderer struct 
	GzRender* newRender = new GzRender();

	// save the pointer to display 
	newRender->display = display;
	
	// setup Xsp and anything only done once 
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			newRender->Xsp[i][j] = 0.0;
	newRender->Xsp[0][0] = ((double)display->xres / 2.0);
	newRender->Xsp[1][1] = (-(double)display->yres / 2.0);
	newRender->Xsp[2][2] = INT_MAX;
	newRender->Xsp[0][3] = ((double)display->xres / 2.0);
	newRender->Xsp[1][3] = ((double)display->yres / 2.0);
	newRender->Xsp[3][3] = 1.0;

	// init default camera 
	newRender->camera.FOV = DEFAULT_FOV;
	newRender->matlevel = 0;

	newRender->camera.position[0] = DEFAULT_IM_X;
	newRender->camera.position[1] = DEFAULT_IM_Y;
	newRender->camera.position[2] = DEFAULT_IM_Z;

	newRender->camera.lookat[0] = 0.0;
	newRender->camera.lookat[1] = 0.0;
	newRender->camera.lookat[2] = 0.0;

	newRender->camera.worldup[0] = 0.0;
	newRender->camera.worldup[1] = 1.0;
	newRender->camera.worldup[2] = 0.0;

	*render = newRender;
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	for (int i = 0; i < render->numlights; i++)
		GzDeleteShadowMapCamera((GzRender*)render->lights_shadow_maps[i]);

	delete render;
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render)
{
/*  
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 

	// - setup for start of each frame - init frame buffer color,alpha,z
	GzInitDisplay(render->display);

	// -Xpi
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			render->camera.Xpi[i][j] = 0.0;
	float d_inv = tan((render->camera.FOV / 2.0) * PI / 180.0);
	render->camera.Xpi[0][0] = 1.0;
	render->camera.Xpi[1][1] = 1.0;
	render->camera.Xpi[3][3] = 1.0;
	render->camera.Xpi[2][2] = d_inv;
	render->camera.Xpi[3][2] = d_inv;

	// -Xiw
	GzInitCameraXiw(render);

	// push
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			render->Ximage[0][i][j] = (i != j) ? 0.0 : 1.0;
			render->Ximage_im[0][i][j] = (i != j) ? 0.0 : 1.0;
			render->Ximage_im[1][i][j] = (i != j) ? 0.0 : 1.0;
			render->Ximage_im[2][i][j] = (i != j) ? 0.0 : 1.0;
			render->Xnorm[0][i][j] = (i != j) ? 0.0 : 1.0;
			render->Xnorm[1][i][j] = (i != j) ? 0.0 : 1.0;
			render->Xnorm[2][i][j] = (i != j) ? 0.0 : 1.0;
		}
	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render, render->camera.Xpi);
	GzPushMatrix(render, render->camera.Xiw);

	// init world to screen Xform:
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++){
			render->Ximage_from_world[0][i][j] = render->Ximage[0][i][j];
			render->Ximage_from_world[1][i][j] = render->Ximage[1][i][j];
			render->Ximage_from_world[2][i][j] = render->Ximage[2][i][j];
			render->Ximage_from_world[3][i][j] = render->Ximage[3][i][j];
		}
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	render->camera = *camera;
	return GZ_SUCCESS;	
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (render->matlevel == 99)
		return GZ_FAILURE;
	render->matlevel++;
	
	GzMatrix matrix_norm;
	if (render->matlevel > 2) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				matrix_norm[i][j] = matrix[i][j];
			}
		}

		matrix_norm[0][3] = 0;
		matrix_norm[1][3] = 0;
		matrix_norm[2][3] = 0;
	
		matrix_norm[3][3] = 1;
	
		matrix_norm[3][0] = 0;
		matrix_norm[3][1] = 0;
		matrix_norm[3][2] = 0;
	}
	float scale = sqrt(matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1] + matrix[0][2] * matrix[0][2]);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix_norm[i][j] = matrix[i][j] / scale;
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			float sum = 0;
			float sum_im = 0;
			float sum_norm = 0;
			for (int k = 0; k < 4; k++) {
				sum +=  render->Ximage[render->matlevel - 1][i][k] * matrix[k][j];
				sum_im +=  render->Ximage_im[render->matlevel - 1][i][k] * matrix[k][j];
				sum_norm +=  render->Xnorm[render->matlevel - 1][i][k] * matrix_norm[k][j];
			}
			render->Ximage[render->matlevel][i][j] = sum; 	
			if (render->matlevel > 2) {
				render->Xnorm[render->matlevel][i][j] = sum_norm;
				render->Ximage_im[render->matlevel][i][j] = sum_im; 	
			}
		}
	}
	
	for (int i = 0; i < render->numlights; i++)
		GzPushMatrix(render->lights_shadow_maps[i], matrix);
	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (render->matlevel == 0)
		return GZ_FAILURE;
	render->matlevel--;

	
	for (int i = 0; i < render->numlights; i++)
		GzPopMatrix(render->lights_shadow_maps[i]);
	return GZ_SUCCESS;
} 


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for (int i = 0; i < numAttributes; i++){
		GzToken token = nameList[i];

		// GZ_RGB_COLOR
		if (token == GZ_RGB_COLOR) {
			GzPointer pointer = valueList[i];

			render->flatcolor[0] = *((float*)(pointer));
			render->flatcolor[1] = *((float*)((pointer)) + 1);
			render->flatcolor[2] = *((float*)((pointer)) + 2);

		}
		if (token == GZ_INTERPOLATE) {
			GzPointer pointer = valueList[i];
			int interpolation_mod = *((int*)(pointer));
			render->interp_mode = interpolation_mod;
		}
		if (token == GZ_DIRECTIONAL_LIGHT) {
			GzPointer pointer = valueList[i];
			GzLight light = *((GzLight*)(pointer));
			render->lights[render->numlights] = light;
			// update with other init procedure for perspective map 
			GzNewShadowMapCamera(&render->lights_shadow_maps[render->numlights], &light, &(render->bbox));
			render->numlights++;
		}
		if (token == GZ_AMBIENT_LIGHT) {
			GzPointer pointer = valueList[i];
			GzLight ambient_light = *((GzLight*)(pointer));
			render->ambientlight = ambient_light;
		}
		if (token == GZ_AMBIENT_COEFFICIENT) {
			GzPointer pointer = valueList[i];
			render->Ka[0] = *((float*)(pointer));
			render->Ka[1] = *((float*)((pointer)) + 1);
			render->Ka[2] = *((float*)((pointer)) + 2);
		}
		if (token == GZ_DIFFUSE_COEFFICIENT) {
			GzPointer pointer = valueList[i];
			render->Kd[0] = *((float*)(pointer));
			render->Kd[1] = *((float*)((pointer)) + 1);
			render->Kd[2] = *((float*)((pointer)) + 2);
		}
		if (token == GZ_SPECULAR_COEFFICIENT) {
			GzPointer pointer = valueList[i];
			render->Ks[0] = *((float*)(pointer));
			render->Ks[1] = *((float*)((pointer)) + 1);
			render->Ks[2] = *((float*)((pointer)) + 2);
		}
		if (token == GZ_DISTRIBUTION_COEFFICIENT) {
			GzPointer pointer = valueList[i];
			float dist = *((float*)(pointer));
			render->spec = dist;
		}
		if (token == GZ_TEXTURE_MAP) {
			GzPointer pointer = valueList[i];
			render->tex_fun = *((GzTexture)(pointer));
		}
		if (token == GZ_AASHIFTX) {
			float shift_x = *((float*)valueList[i]);
			render->display->x_shift = shift_x;
		}
		if (token == GZ_AASHIFTY) {
			float shift_y = *((float*)valueList[i]);
			render->display->y_shift = shift_y;
		}
		if (token == GZ_BBOX) {
			GzBoundingBox bbox = *((GzBoundingBox*)valueList[i]);
			render->bbox = bbox;
		}
	}
	
	return GZ_SUCCESS;
}

void multiplyMatrixByVector(float pX, float pY, float pZ, GzMatrix matrix, float* x, float* y, float* z, float* w) {
	*x = pX * matrix[0][0] + pY * matrix[0][1] + pZ * matrix[0][2] + 1 * matrix[0][3];
	*y = pX * matrix[1][0] + pY * matrix[1][1] + pZ * matrix[1][2] + 1 * matrix[1][3];
	*z = pX * matrix[2][0] + pY * matrix[2][1] + pZ * matrix[2][2] + 1 * matrix[2][3];
	*w = pX * matrix[3][0] + pY * matrix[3][1] + pZ * matrix[3][2] + 1 * matrix[3][3];
	return;
}

void retrieveData( GzPointer *valueList, int i, float* p1X, float* p1Y, float* p1Z, float* p2X, float* p2Y, float* p2Z, float* p3X, float* p3Y, float* p3Z) {
	GzCoord* p1 = (GzCoord*)(valueList[i]);
	GzCoord* p2 = ((GzCoord*)valueList[i] + 1);
	GzCoord* p3 = ((GzCoord*)valueList[i] + 2);
			
			
	*p1X = (*p1)[0];
	*p1Y = (*p1)[1];
	*p1Z = (*p1)[2];

	*p2X = (*p2)[0];
	*p2Y = (*p2)[1];
	*p2Z = (*p2)[2];

	*p3X = (*p3)[0];
	*p3Y = (*p3)[1];
	*p3Z = (*p3)[2];	
	return;
}
void retrieveData( GzPointer *valueList, int i, float* p1X, float* p1Y,  float* p2X, float* p2Y,  float* p3X, float* p3Y) {
	GzTextureIndex* p1 = ( GzTextureIndex*)(valueList[i]);
	GzTextureIndex* p2 = ((GzTextureIndex*)valueList[i] + 1);
	GzTextureIndex* p3 = ((GzTextureIndex*)valueList[i] + 2);
			
	*p1X = (*p1)[0];
	*p1Y = (*p1)[1];

	*p2X = (*p2)[0];
	*p2Y = (*p2)[1];

	*p3X = (*p3)[0];
	*p3Y = (*p3)[1];	
	return;
}


void GzNormalize_vector(float* x, float* y, float* z) {
	float norm = sqrt(*x * *x + *y * *y + *z * *z);
	*x /= norm;
	*y /= norm;
	*z /= norm;
	return;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts using matrix on top of stack 
- Clip - just discard any triangle with any vert(s) behind view plane 
       - optional: test for triangles with all three verts off-screen (trivial frustum cull)
- invoke triangle rasterizer  
*/ 
	Triangle tr;
	bool is_position_set = false;
	bool is_normal_set = false;
	bool is_texture_set = false;
	float x1, y1, z1, x2, y2, z2, x3, y3, z3;
	float x1_im, y1_im, z1_im, x2_im, y2_im, z2_im, x3_im, y3_im, z3_im;
	float nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3;
	float p1U, p1V, p2U, p2V, p3U, p3V;

	for (int i = 0; i < numParts; i++){
		GzToken token = nameList[i];

		// GZ_POSITION
		if (token == GZ_POSITION) {			
			float p1X, p1Y, p1Z, p2X, p2Y, p2Z, p3X, p3Y, p3Z;
			retrieveData(valueList, i, &p1X, &p1Y, &p1Z, &p2X, &p2Y, &p2Z, &p3X, &p3Y, &p3Z);
			
			float w1, w2, w3;
			multiplyMatrixByVector(p1X, p1Y, p1Z, render->Ximage[render->matlevel], &x1, &y1, &z1, &w1);
			multiplyMatrixByVector(p2X, p2Y, p2Z, render->Ximage[render->matlevel], &x2, &y2, &z2, &w2);
			multiplyMatrixByVector(p3X, p3Y, p3Z, render->Ximage[render->matlevel], &x3, &y3, &z3, &w3);

			if ((z1 < 0 || z2 < 0 || z3 < 0) ||	(w1 == 0 || w2 == 0 || w3 == 0))
				continue;
			
			x1 /= w1;
			y1 /= w1;
			z1 /= w1;
			x2 /= w2;
			y2 /= w2;
			z2 /= w2;
			x3 /= w3;
			y3 /= w3;
			z3 /= w3;

			if ((x1 < 0 && x2 < 0 && x3 < 0) || (x1 > render->display->xres && x2 > render->display->xres && x3 > render->display->xres) ||
				(y1 < 0 && y2 < 0 && y3 < 0) || (y1 > render->display->yres && y2 > render->display->yres && y3 > render->display->yres)
				)
					continue;	

			multiplyMatrixByVector(p1X, p1Y, p1Z, render->Ximage_im[render->matlevel], &x1_im, &y1_im, &z1_im, &w1);
			multiplyMatrixByVector(p2X, p2Y, p2Z, render->Ximage_im[render->matlevel], &x2_im, &y2_im, &z2_im, &w2);
			multiplyMatrixByVector(p3X, p3Y, p3Z, render->Ximage_im[render->matlevel], &x3_im, &y3_im, &z3_im, &w3);
			
			x1_im /= w1;
			y1_im /= w1;
			z1_im /= w1;
			x2_im /= w2;
			y2_im /= w2;
			z2_im /= w2;
			x3_im /= w3;
			y3_im /= w3;
			z3_im /= w3;

			is_position_set = true;
		}

		
		// GZ_NORMAL
		if (token == GZ_NORMAL) {
			float p1X, p1Y, p1Z, p2X, p2Y, p2Z, p3X, p3Y, p3Z;
			retrieveData(valueList, i, &p1X, &p1Y, &p1Z, &p2X, &p2Y, &p2Z, &p3X, &p3Y, &p3Z);
			
			GzNormalize_vector(&p1X, &p1Y, &p1Z);
			GzNormalize_vector(&p2X, &p2Y, &p2Z);
			GzNormalize_vector(&p3X, &p3Y, &p3Z);

			float nw1, nw2, nw3;
			multiplyMatrixByVector(p1X, p1Y, p1Z, render->Xnorm[render->matlevel], &nx1, &ny1, &nz1, &nw1);
			multiplyMatrixByVector(p2X, p2Y, p2Z, render->Xnorm[render->matlevel], &nx2, &ny2, &nz2, &nw2);
			multiplyMatrixByVector(p3X, p3Y, p3Z, render->Xnorm[render->matlevel], &nx3, &ny3, &nz3, &nw3);

			nx1 /= nw1;
			ny1 /= nw1;
			nz1 /= nw1;
			nx2 /= nw2;
			ny2 /= nw2;
			nz2 /= nw2;
			nx3 /= nw3;
			ny3 /= nw3;
			nz3 /= nw3;

			is_normal_set = true;
		}


		// GZ_TEXTURE_INDEX
		if (token == GZ_TEXTURE_INDEX) {			
			retrieveData(valueList, i, &p1U, &p1V, &p2U, &p2V, &p3U, &p3V);
			is_texture_set = true;
		}

		// GZ_NULL_TOKEN
		if (token == GZ_NULL_TOKEN) {
			// do nothing;
		}

		if (is_normal_set && is_position_set && is_texture_set ) {
			tr.init_triangle(x1 + render->display->x_shift, y1 + render->display->y_shift, z1, 
							x2 + render->display->x_shift, y2 + render->display->y_shift, z2, 
							x3 + render->display->x_shift, y3 + render->display->y_shift, z3, 
							x1_im, y1_im, z1_im, x2_im, y2_im, z2_im, x3_im, y3_im, z3_im, 
							nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, 
							p1U, p1V, p2U, p2V, p3U, p3V, render);
			tr.rasterize(render);	
			is_normal_set = false;
			is_position_set = false;
		}
	}
	return GZ_SUCCESS;
}
