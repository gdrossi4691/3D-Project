#include "disp.h" /* include your own disp.h file (e.g. hw1)*/

/* Camera defaults */
#define	DEFAULT_FOV	35.0/* world coords for image plane origin */
#define	DEFAULT_IM_Z	(-10.0)  
#define	DEFAULT_IM_Y	(5.0)    
#define	DEFAULT_IM_X	(-10.0)/* default look-at point = 0,0,0 */

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100	/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10	/* how many lights allowed */

#define PI 3.141592653589793


#pragma once
typedef struct {
	float Xmin;
	float Xmax;
	float Ymin;
	float Ymax;
	float Zmin;
	float Zmax;
} GzBoundingBox;

#ifndef GZRENDER
#define GZRENDER
typedef struct GzRender {			/* define a renderer */
  GzDisplay	*display;
  GzCamera		camera;
  short		matlevel;  /* top of stack -> current xform */
  GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
  GzMatrix		Ximage_im[MATLEVELS];	/* stack of xforms (Xsm) */
  GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
  GzMatrix		Xsp;		 /* NDC to screen (pers-to-screen) */
  GzColor		flatcolor;  /* color for flat shaded triangles */
  int			interp_mode;
  int			numlights;
  GzLight		lights[MAX_LIGHTS];
//  struct GzRender*	lights_shadow_maps[MAX_LIGHTS];
  struct GzRender*	lights_shadow_maps[MAX_LIGHTS];
  GzLight		ambientlight;
  GzColor		Ka, Kd, Ks;
  float		    spec;		/* specular power */
  GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
  float			shift_x;
  float			shift_y;
  GzBoundingBox bbox;
}  GzRender;
#endif


short	ctoi(float color);
void multiplyMatrixByVector(float pX, float pY, float pZ, GzMatrix matrix, float* x, float* y, float* z, float* w); 
void GzInitCameraXiw(GzRender *render);
// Function declaration
// HW2
int GzNewRender(GzRender **render, GzDisplay *display);
int GzUpdateRender(GzRender *render, GzDisplay	*display);
int GzFreeRender(GzRender *render);
int GzBeginRender(GzRender	*render);
int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList);
int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList);

// HW3
int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix	matrix);
int GzPopMatrix(GzRender *render);

// HW5
int GzFreeTexture();
int tex_fun(float u, float v, GzColor color);
int ptex_fun(float u, float v, GzColor color);

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);


// HW7

// shadow map initialization applicable for both visibility functions (step I)
int GzNewShadowMapCamera(GzRender** map, GzLight* light, GzBoundingBox* bbox);
// shadow map initialization applicable only for perspective shadow mapping (step II)
int GzNewPerspectiveShadowMapCamera(GzRender** map, GzLight* light, GzBoundingBox bbox);
// free memory
int GzDeleteShadowMapCamera(GzRender* map);

// Percentage Closer Filter. Filter has a fixed-size kernel. x,y,z are coordinates of the point in world.
float GzPCFVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light);
// Percentage Closer Filter with with dynamic kernel size. Penumbra size depends on light area size. x,y,z are coordinates of the point in world.
float GzPCFSoftShadowVisibilityFn(float x, float y, float z, GzRender* map, GzLight* light);



void GzNormalize_vector(float* x, float* y, float* z);
#pragma once
class Triangle {
private:
	class Edge
        {
         public:
                  double x1;
                  double y1;
                  double z1;

				  double x2;
                  double y2;
                  double z2;

				  double nx1;
                  double ny1;
                  double nz1;

				  double nx2;
                  double ny2;
                  double nz2;

				  bool isL;
				  bool isR;
				  bool isT;
				  bool isB;

				  double A;
				  double B;
				  double C;
		public :
			void init_edge(bool is_bottom, bool is_top, bool is_left, bool is_right, 
					double x1, double y1, double z1, double x2, double y2, double z2,
					double nx1, double ny1, double nz1, double nx2, double ny2, double nz2){
				this->x1 = x1;
				this->y1 = y1;
				this->z1 = z1;
				this->x2 = x2;
				this->y2 = y2;
				this->z2 = z2;

				this->nx1 = nx1;
				this->ny1 = ny1;
				this->nz1 = nz1;
				this->nx2 = nx2;
				this->ny2 = ny2;
				this->nz2 = nz2;

				this->isB = is_bottom;
				this->isT = is_top;
				this->isL = is_left;
				this->isR = is_right;

				// line equation:
				this->A =	this->y2 - this->y1;
				this->B = -(this->x2 - this->x1);
				this->C = (	this->x2 - this->x1) - (this->y2 - this->y1);
				
				return;
			}
			int sign(int x, int y) {
				double res = -(((y2 - y1) * (x - x1) - (x2 - x1) * (y - y1)));
				if (res == 0 ) {
					if (isL)
						return 1;
					if (isR)
						return -1;
					if (isT)
						return 1;
					if (isB)
						return -1;
					return GZ_FAILURE;
				}
				if (res > 0) {
					return 1;
				}
				if (res < 0) {
					return -1;
				}
				return GZ_FAILURE;
			}
        };

	class BoundingBox {
		public :
			int minX;
			int minY;
			int maxX;
			int maxY;
		};

	class Plane {
		public :
			double A;
			double B;
			double C;
			double D;
		};

private :
	Edge edges[3];
	BoundingBox b_box;
	
	Plane plane;
	Plane plane_im;

	Plane color_plane_r;
	Plane color_plane_g;
	Plane color_plane_b;

	Plane color_plane_r_no_txt;
	Plane color_plane_g_no_txt;
	Plane color_plane_b_no_txt;

	Plane normal_plane_x;
	Plane normal_plane_y;
	Plane normal_plane_z;

	Plane text_plane_u;
	Plane text_plane_v;

	bool is_flat;
private : 
	void init_edges(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float n1X, float n1Y, float n1Z, float n2X, float n2Y, float n2Z, float n3X, float n3Y, float n3Z);

	void init_orientation(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, 
		double nx1, double ny1, double nz1, double nx2, double ny2, double nz2, double nx3, double ny3, double nz3);

	bool is_flat_triangle(double x1, double y1, double x2, double y2, double x3, double y3);

	BoundingBox init_bounding_box(double p1X, double p1Y, double p2X, double p2Y, double p3X, double p3Y);

	Plane init_plane(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);

	int calculate_sign(int x, int y /*Edge edges[3]*/);

	int int_interpolate(int x, int y, Plane plane);

	float float_interpolate(double x, double y, Plane plane);

	float multipy_vectors(float x1, float y1, float z1, float x2, float y2, float z2);

	void compute_color(double x_im, double y_im, double z_im, GzColor* color, GzRender *render, bool is_flat);

	void compute_gouraud_color(double x_im, double y_im, double z_im, GzColor* color, GzRender *render);

public :
	void init_triangle(float p1X, float p1Y, float p1Z, float p2X, float p2Y, float p2Z, float p3X, float p3Y, float p3Z,    
		float p1X_im, float p1Y_im, float p1Z_im, float p2X_im, float p2Y_im, float p2Z_im, float p3X_im, float p3Y_im, float p3Z_im,    
		float n1X, float n1Y, float n1Z, float n2X, float n2Y, float n2Z, float n3X, float n3Y, float n3Z, 
		float p1U, float p1V, float p2U, float p2V, float p3U, float p3V, GzRender* render);
	
	void rasterize(GzRender *render);
};

