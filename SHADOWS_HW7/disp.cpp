/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"
#include <iostream>
#include <fstream>

using namespace std;


int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
 -- NOTE: this function is optional and not part of the API, but you may want to use it within the display function.
*/
	char* frmbf = new char[3 * (width * height)];
	*framebuffer = frmbf;
	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	GzDisplay* disp = new GzDisplay;
	disp->xres = xRes;
	disp->yres = yRes;
	disp->fbuf = new GzPixel[xRes * yRes];
	*display = disp;
	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* clean up, free memory */
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;
	delete[] display->fbuf;
	delete display;
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* pass back values for a display */
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;
	*xRes = display->xres;
	*yRes = display->yres;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;
	for (int i = 0; i < display->yres; i++) {
		for (int j = 0; j < display->xres; j++){
			display->fbuf[ARRAY(j,i)].alpha = 0;
			display->fbuf[ARRAY(j,i)].z = INT_MAX;
			display->fbuf[ARRAY(j,i)].red = 0;
			display->fbuf[ARRAY(j,i)].green = 3000;
			display->fbuf[ARRAY(j,i)].blue = 1000;	
		}
	}
	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;

	if (i >= display->xres)
		return GZ_SUCCESS;
	if (j >= display->yres)
		return GZ_SUCCESS;
	if (i < 0)
		return GZ_SUCCESS;
	if (j < 0)
		return GZ_SUCCESS;

	display->fbuf[ARRAY(i,j)].alpha = a;
	GzDepth curZ = display->fbuf[ARRAY(i,j)].z;
	if ((z < curZ) && (z > 0)){
		display->fbuf[ARRAY(i,j)].z = z;
		display->fbuf[ARRAY(i,j)].red = r;
		display->fbuf[ARRAY(i,j)].green = g;
		display->fbuf[ARRAY(i,j)].blue = b;
	}
	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;

	*a = display->fbuf[ARRAY(i,j)].alpha;
	*z = display->fbuf[ARRAY(i,j)].z;
	*r = display->fbuf[ARRAY(i,j)].red;
	*g = display->fbuf[ARRAY(i,j)].green;
	*b = display->fbuf[ARRAY(i,j)].blue;
	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
	/* write pixels to ppm file -- "P6 %d %d 255\r" */
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;
	fprintf(outfile, "P6 %d %d 255\r", display->xres, display->yres);
	for (int i = 0; i < display->yres; i++) {
		for (int j = 0; j < display->xres; j++){
			GzIntensity ib = display->fbuf[ARRAY(j,i)].blue ;
			GzIntensity ig = display->fbuf[ARRAY(j,i)].green;
			GzIntensity ir = display->fbuf[ARRAY(j,i)].red  ;

			GzIntensity cl_b = (ib > 4095) ? 4095 : ib;
			GzIntensity cl_g = (ig > 4095) ? 4095 : ig;
			GzIntensity cl_r = (ir > 4095) ? 4095 : ir;

			unsigned char b = (unsigned char) (cl_b >> 4);
			unsigned char g = (unsigned char) (cl_g >> 4);
			unsigned char r = (unsigned char) (cl_r >> 4);

			fwrite(&r, sizeof(unsigned char), 1, outfile);
			fwrite(&g, sizeof(unsigned char), 1, outfile);
			fwrite(&b, sizeof(unsigned char), 1, outfile);

		}
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{
	/* write pixels to framebuffer: 
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
		- Not red, green, and blue !!!
	*/
	if ((display == NULL)  || (display->fbuf == NULL))
		return GZ_FAILURE;

	for (int i = 0; i < display->yres; i++) {
		for (int j = 0; j < display->xres; j++){
			GzIntensity ib = display->fbuf[ARRAY(j,i)].blue ;
			GzIntensity ig = display->fbuf[ARRAY(j,i)].green;
			GzIntensity ir = display->fbuf[ARRAY(j,i)].red  ;

			GzIntensity cl_b = (ib > 4095) ? 4095 : ib;
			GzIntensity cl_g = (ig > 4095) ? 4095 : ig;
			GzIntensity cl_r = (ir > 4095) ? 4095 : ir;

			unsigned char b = (unsigned char) (cl_b >> 4);
			unsigned char g = (unsigned char) (cl_g >> 4);
			unsigned char r = (unsigned char) (cl_r >> 4);

			framebuffer[ARRAY(j,i) * 3 + 0] = b ;
			framebuffer[ARRAY(j,i) * 3 + 1] = g ;
			framebuffer[ARRAY(j,i) * 3 + 2] = r ;

		}
	}

	return GZ_SUCCESS;
}