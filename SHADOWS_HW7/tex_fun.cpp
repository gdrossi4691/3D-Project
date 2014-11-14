/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include	<cmath>
#include	<cstdio>
#include	<limits>

GzColor	*image=NULL;
int xs, ys;
int reset = 1;
int proc_reset = 1;

int interpolate(float u, float v, GzColor color) {
	if (u > 1){
		u = 1;
//		return GZ_FAILURE;
	}
	if (u < 0){
		u = 0;
//		return GZ_FAILURE;
	}
	if (v > 1){
		v = 1;
//		return GZ_FAILURE;
	}
	if (v < 0){
		v = 0;
//		return GZ_FAILURE;
	}
	float x = u * (xs - 1);
	float y = v * (ys - 1);

	int A_x, A_y, B_x, B_y, C_x, C_y, D_x, D_y;
	A_x = (int)x;
	A_y = (int)y;
	C_x = (int)x + 1;
	C_y = (int)y + 1;
	B_x = (int)x + 1;
	B_y = (int)y;
	D_x = (int)x;
	D_y = (int)y + 1;
	
	// Color(p) = s t C + (1-s) t D + s (1-t) B + (1-s) (1-t) A 
	float s = x - A_x;
	float t = y - A_y;
	color[RED]	= s * t * image[C_y * xs + C_x][RED]	+ (1 - s) * t * image[D_y * xs + D_x][RED]	+ s * (1 - t) *  image[B_y * xs + B_x][RED]	 + (1 - s) * (1 - t) *  image[A_y * xs + A_x][RED]	;
	color[GREEN]= s * t * image[C_y * xs + C_x][GREEN]	+ (1 - s) * t * image[D_y * xs + D_x][GREEN]+ s * (1 - t) *  image[B_y * xs + B_x][GREEN]+ (1 - s) * (1 - t) *  image[A_y * xs + A_x][GREEN]	;
	color[BLUE] = s * t * image[C_y * xs + C_x][BLUE]	+ (1 - s) * t * image[D_y * xs + D_x][BLUE]	+ s * (1 - t) *  image[B_y * xs + B_x][BLUE] + (1 - s) * (1 - t) *  image[A_y * xs + A_x][BLUE]	;
	
	return GZ_SUCCESS;
}


/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */
	return interpolate(u, v, color);
}

#define C_REAL  0.4 
#define C_IMG  -0.25
#define POLYNOM_POWER  2
#define MAX_LEVEL  4

double start_u;
double start_v;

void to_power(double* u, double* v, int power) {
	float x = *u;
	float x_i = *v;

	float y, y_i;
	for (int i = 0; i < power; i++) {
		y  = x * *u - x_i * *v;
		y_i = x * *v + x_i * *u; 
		x = y;
		x_i = y_i;
	}
	*u = x;
	*v = x_i;
}

void julia(double * u, double* v) {
	static int level = 0;
	if (level == MAX_LEVEL) {
		level = 0;
		*u = start_u;
		*v = start_v;
		return;
	} else {
		level++;
		julia(u, v);
		double u2 = *u;
		double v2 = *v;
		to_power(&u2, &v2, POLYNOM_POWER);
		u2 = u2 + C_REAL;
		v2 = v2 + C_IMG;
		if ((u2 == u2) && fabs(u2) != std::numeric_limits<double>::infinity( ))
			*u = u2;
		if ((v2 == v2) && fabs(v2) != std::numeric_limits<double>::infinity( ))
			*v = v2;
		if (u2 != u2)
			return;
		if (v2 != v2){
			return;
		}
	}
}


void to_color(double length, GzColor color) {
	static float colors[11] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
	static float colors_r[11] = {1.0, 55.00/255.0, 255.0/255.0, 0.000/255.0, 112.0/255.0, 153.0/255.0, 255.0/255.0, 0.000/255.0, 0.000/255.0, 166.0/255.0, 255.0/255.0};
	static float colors_g[11] = {0.0, 255.0/255.0, 255.0/255.0, 112.0/255.0, 48.00/255.0, 255.0/255.0, 192.0/255.0, 112.0/255.0, 51.00/255.0, 166.0/255.0, 0.000/255.0};
	static float colors_b[11] = {0.0, 55.00/255.0, 0.000/255.0, 192.0/255.0, 160.0/255.0, 153.0/255.0, 1.000/255.0, 192.0/255.0, 102.0/255.0, 166.0/255.0, 255.0/255.0};


	if (length > 2)
		length = 2;
	
	int i = 0;
	while (length > colors[i] && i < 11) 
		i++;
	i--;


	color[RED]	= colors_r[i] * (colors[i + 1] - length) / 0.2 + colors_r[i + 1] * (length - colors[i]) / 0.2;
	color[GREEN]= colors_g[i] * (colors[i + 1] - length) / 0.2 + colors_g[i + 1] * (length - colors[i]) / 0.2;
	color[BLUE] = colors_b[i] * (colors[i + 1] - length) / 0.2 + colors_b[i + 1] * (length - colors[i]) / 0.2;
}

#define K 10

double turbulence(double u, double v) {
	double len = sqrt(u*u + v*v);
	double sum = 0;
	for (int i = 0; i < K; i++) {
		double two = 2;
		for (int j = 0; j < i; j++)
			two *=2;
		int rand_ = rand(); 
		int two_x = (int)(two * len + 1);
		double noise = rand_ % two_x;
		noise /= two;
		sum += noise;
	}
	return sum;
}

/* Procedural texture function */
int ptex_fun(float u_start, float v_start, GzColor color)
{
	/*
	if (proc_reset) {          // create texture
		proc_reset = 0;
		xs = 200;
		ys = 200;
		image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
		if (image == NULL) {
			fprintf (stderr, "malloc for texture image failed\n");
			exit(-1);
		}
		for (int i = 0; i < xs*ys; i++) {	// create array of GzColor values 
			start_u = ((int)(i - ((int)i / (int)ys) * ys)) / (double)xs;
			start_v = ((int)((int)i / (int)ys)) / (double)ys;
			double u,v;
			julia(&u, &v);
			double length = sqrt(u * u + v * v) / 100.0;
			to_color(length, color);

			image[i][RED] =	color[RED]	;
			image[i][GREEN] = color[GREEN];
			image[i][BLUE]  = color[BLUE]	;
		}
	}

	return interpolate(u_start, v_start, color);
	*/
	start_u = (u_start > 0.5) ? u_start : (1.0 - u_start );
	start_v = (v_start > 0.5) ? v_start : (1.0 - v_start );
	double u,v;
	julia(&u, &v);
	double length = sqrt(u * u + v * v);
	to_color(length, color);
	return GZ_SUCCESS; 
	
}




/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

