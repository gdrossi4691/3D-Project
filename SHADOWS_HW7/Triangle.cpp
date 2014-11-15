#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"disp.h"


// private 
void Triangle::init_edges(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float n1X, float n1Y, float n1Z, float n2X, float n2Y, float n2Z, float n3X, float n3Y, float n3Z) {	
	// sort edges by Y:
	is_flat = is_flat_triangle(x1, y1, x2, y2, x3, y3); // special case - check if all three points on one line.
	if (is_flat)
		return;
	if ( (y1 <= y2) && (y2 <= y3))
		init_orientation(x1, y1, z1, x2, y2, z2, x3, y3, z3, n1X, n1Y, n1Z, n2X, n2Y, n2Z, n3X, n3Y, n3Z);
	if ( (y1 <= y3) && (y3 <= y2)) 
		init_orientation(x1, y1, z1, x3, y3, z3, x2, y2, z2, n1X, n1Y, n1Z, n3X, n3Y, n3Z, n2X, n2Y, n2Z);
	if ( (y2 <= y1) && (y1 <= y3)) 
		init_orientation(x2, y2, z2, x1, y1, z1, x3, y3, z3, n2X, n2Y, n2Z, n1X, n1Y, n1Z, n3X, n3Y, n3Z);
	if ( (y2 <= y3) && (y3 <= y1)) 
		init_orientation(x2, y2, z2, x3, y3, z3, x1, y1, z1, n2X, n2Y, n2Z, n3X, n3Y, n3Z, n1X, n1Y, n1Z);
	if ( (y3 <= y1) && (y1 <= y2)) 
		init_orientation(x3, y3, z3, x1, y1, z1, x2, y2, z2, n3X, n3Y, n3Z, n1X, n1Y, n1Z, n2X, n2Y, n2Z);
	if ( (y3 <= y2) && (y2 <= y1)) 
		init_orientation(x3, y3, z3, x2, y2, z2, x1, y1, z1, n3X, n3Y, n3Z, n2X, n2Y, n2Z, n1X, n1Y, n1Z);
		
	return;
}

void Triangle::init_orientation(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double nx1, double ny1, double nz1, double nx2, double ny2, double nz2, double nx3, double ny3, double nz3){
	// CW orientation of the edges (needs flipping! (x(-1)).
	if (y1 == y2) {// (1-2) is a horizontal edge, and (1-2) is Top, and (3-1) is Left and (2-3) is Right => order is (1-2) (2-3) (3-1) 
		if (x2 < x1) { // switch x1 and x2, otherwise it is inconsistent (counterclockwise)
			int tmp = x2;
			x2 = x1;
			x1 = tmp;
		}
		// (1-2)
		edges[0].init_edge(false, true, false, false, x1, y1, z1, x2, y2, z2, nx1, ny1, nz1, nx2, ny2, nz2);
		// (2-3)
		edges[1].init_edge(false, false, false, true,  x2, y2, z2, x3, y3, z3, nx2, ny2, nz2, nx3, ny3, nz3);
		// (3-1)
		edges[2].init_edge(false, false, true, false, x3, y3, z3, x1, y1, z1, nx3, ny3, nz3, nx1, ny1, nz1);

		return;
	}
	if (y2 == y3) {// (2-3) is a horizontal edge, and (2-3) is Bottom, and (3-1) is Left and (1-2) is Right => order is (1-2) (2-3) (3-1) 
		if (x3 > x2) { // switch x3 and x2, otherwise it is inconsistent (counterclockwise)
			int tmp = x3;
			x3 = x2;
			x2 = tmp;
		}
		// (1-2)
		edges[0].init_edge(false, false, false, true, x1, y1, z1, x2, y2, z2, nx1, ny1, nz1, nx2, ny2, nz2);
		// (2-3)
		edges[1].init_edge(true, false, false, false,  x2, y2, z2, x3, y3, z3, nx2, ny2, nz2, nx3, ny3, nz3);
		// (3-1)
		edges[2].init_edge(false, false, true, false, x3, y3, z3, x1, y1, z1, nx3, ny3, nz3, nx1, ny1, nz1);

		return;
	}	
	double midX = x1 + ((x3 - x1) / (y3 - y1)) * (y2 - y1); //mid_x(y2, x1, y1, x3, y3);// expect y2 != y1

	if (midX > x2){ // (2-1) and (3-2) are Lefft edges, and (1-3) is Right => order is (1-3) (3-2) (2-1)
		// (1-3)
		edges[0].init_edge(false, false, false, true, x1, y1, z1, x3, y3, z3, nx1, ny1, nz1, nx3, ny3, nz3);
		// (3-2)
		edges[1].init_edge(false, false, true, false, x3, y3, z3, x2, y2, z2, nx3, ny3, nz3, nx2, ny2, nz2);
		// (2-1)
		edges[2].init_edge(false, false, true, false, x2, y2, z2, x1, y1, z1, nx2, ny2, nz2, nx1, ny1, nz1);

		return;
	}
	if (midX < x2){ // (1-2) and (2-3) are Right edges, and (3-1) is Left => order is (1-2) (2-3) (3-1)
		// (1-2)
		edges[0].init_edge(false, false, false, true, x1, y1, z1, x2, y2, z2, nx1, ny1, nz1, nx2, ny2, nz2);
		// (2-3)
		edges[1].init_edge(false, false, false, true, x2, y2, z2, x3, y3, z3, nx2, ny2, nz2, nx3, ny3, nz3);
		// (3-1)
		edges[2].init_edge(false, false, true, false, x3, y3, z3, x1, y1, z1, nx3, ny3, nz3, nx1, ny1, nz1);
			
		return;
	}
}

bool Triangle::is_flat_triangle(double x1, double y1, double x2, double y2, double x3, double y3){
	if ( (y1 == y2) && (y2 == y3)) // special case - check if all three pints on one horizontal line.
		return true;
	else
		if ( (x3 - x1) * (y2 - y1) == (y3 - y1) * (x2 - x1) )
			return true;
		else 
			return false;
	return false;
}

Triangle::BoundingBox Triangle::init_bounding_box(double p1X, double p1Y, double p2X, double p2Y, double p3X, double p3Y) {
	BoundingBox box;
	box.minX = min(min(p1X, p2X), p3X) - 1;
	box.maxX = max(max(p1X, p2X), p3X) + 1;
	box.minY = min(min(p1Y, p2Y), p3Y) - 1;
	box.maxY = max(max(p1Y, p2Y), p3Y) + 1;
	return box;
}

Triangle::Plane Triangle::init_plane(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
	Plane pl;
	pl.A = (y2 - y1) * (z3 - z1) - (y3 -y1) * (z2 - z1);
	pl.B = -(x2 - x1) * (z3 - z1) + (x3 -x1) * (z2 - z1);
	pl.C = (x2 - x1) * (y3 - y1) - (x3 -x1) * (y2 - y1);
	pl.D = -x1 * pl.A - y1 * pl.B - z1 * pl.C;
	//double z = -(pl.A * 252.0 + pl.B * 151.0 + pl.D) / pl.C;
	return pl;
}

int Triangle::calculate_sign(int x, int y /*Edge edges[3]*/){
	int edg1 = edges[0].sign(x, y);
	int edg2 = edges[1].sign(x, y);
	int edg3 = edges[2].sign(x, y);

	if (((edg1 > 0) && (edg2 > 0) && (edg3 > 0)) ) {
		return 1;
	}
	else {
		return -1;
	}
}

int Triangle::int_interpolate(int x, int y, Plane plane) {
	if (plane.C != 0.0) {
		double z = -(plane.A * x + plane.B * y + plane.D) / plane.C;
		//z = (z >= 0) ? z + 0.5 : z - 0.5;
		int int_z = (z > INT_MAX) ? (INT_MAX - 1) : static_cast<int>(z);
		return int_z;
	} else { // C == 0
		return INT_MAX;
	}
	return GZ_FAILURE;
}

float Triangle::float_interpolate(double x, double y, Plane plane) {
	if (plane.C != 0.0) {
		double z = -(plane.A * x + plane.B * y + plane.D) / plane.C;
		float fl_z = static_cast<float>(z);
		return fl_z ;
	} else { // C == 0
		return INT_MAX;
	}
	return GZ_FAILURE;
}

float Triangle::multipy_vectors(float x1, float y1, float z1, float x2, float y2, float z2) {
	return x1 * x2 + y1 * y2 + z1 * z2;
}

void Triangle::compute_color(double x_im, double y_im, double z_im, GzColor* color, GzRender *render, bool is_flat) {
	float nx = float_interpolate(x_im, y_im, normal_plane_x);
	float ny = float_interpolate(x_im, y_im, normal_plane_y);
	float nz = float_interpolate(x_im, y_im, normal_plane_z);

	float u = float_interpolate(x_im, y_im, text_plane_u);
	float v = float_interpolate(x_im, y_im, text_plane_v);
	if (is_flat) {
		float nx_f = (float)plane_im.A;
		float ny_f = (float)plane_im.B;
		float nz_f = (float)plane_im.C;
		float N_NF_dot_product = multipy_vectors(nx, ny, nz, nx, ny, nz);
		nx = (N_NF_dot_product < 0) ? -nx_f : nx_f;
		ny = (N_NF_dot_product < 0) ? -ny_f : ny_f;
		nz = (N_NF_dot_product < 0) ? -nz_f : nz_f;
	}
	GzNormalize_vector(&nx, &ny, &nz);
		
	// implementing eye vector directed to the focal point
	float d = 1.0 / tan((render->camera.FOV / 2.0) * (PI / 180.0));
	// FIXME: don't know why "ex = -x_im / 2.0;" works beter then "ex = -x_im;" ???
	float ex = -x_im / 2.0;		//0;
	float ey = -y_im / 2.0;		//0;
	float ez = -d -z_im;		//-1;
	GzNormalize_vector(&ex, &ey, &ez);

	GzColor amb_color;
	if (render->tex_fun == NULL){
		amb_color[0] = render->Ka[0];
		amb_color[1] = render->Ka[1];
		amb_color[2] = render->Ka[2];
	} else {
		int status = render->tex_fun(u, v, amb_color);
		if (status) 
			exit(GZ_FAILURE); 
	}
	(*color)[0] = amb_color[0] * render->ambientlight.color[0];
	(*color)[1] = amb_color[1] * render->ambientlight.color[1];
	(*color)[2] = amb_color[2] * render->ambientlight.color[2];

	// TODO: compute world coordinates of curent point
	float x = 0.0; 
	float y = 0.0;
	float z = 0.0;
	float w = 0.0;
	multiplyMatrixByVector(x_im, y_im, z_im, render->Xwi, &x, &y, &z, &w);
	
	for (int i = 0; i < render->numlights; i++) {
		float visibility = GzPCFVisibilityFn(x, y, z, render->lights_shadow_maps[i], &render->lights[i]);

		float lx = render->lights[i].position[0] - x_im;
		float ly = render->lights[i].position[1] - y_im;
		float lz = render->lights[i].position[2] - z_im;
		GzNormalize_vector(&lx, &ly, &lz);

		float N_L_dot_product = multipy_vectors(nx, ny, nz, lx, ly, lz);
		float N_E_dot_product = multipy_vectors(nx, ny, nz, ex, ey, ez);
		if ((N_L_dot_product  < 0 && N_E_dot_product > 0) || (N_L_dot_product > 0 && N_E_dot_product < 0))
			continue;
		if (N_L_dot_product  < 0 && N_E_dot_product < 0){
			nx *= -1;
			ny *= -1;
			nz *= -1;
			N_L_dot_product = multipy_vectors(nx, ny, nz, lx, ly, lz);
			N_E_dot_product = multipy_vectors(nx, ny, nz, ex, ey, ez);
		}

		float rx = 2 * nx * N_L_dot_product - lx;
		float ry = 2 * ny * N_L_dot_product - ly;
		float rz = 2 * nz * N_L_dot_product - lz;
		GzNormalize_vector(&rx, &ry, &rz);

		float R_E_dot_product = pow(multipy_vectors(rx, ry, rz, ex, ey, ez), render->spec);
		R_E_dot_product = R_E_dot_product < 0 ? 0 : R_E_dot_product;

		if (render->tex_fun == NULL) {
			amb_color[0] = render->Kd[0];
			amb_color[1] = render->Kd[1];
			amb_color[2] = render->Kd[2];
		}

		(*color)[0] += visibility * (render->Ks[0] * render->lights[i].color[0] * R_E_dot_product + amb_color[0] * render->lights[i].color[0] * N_L_dot_product);
		(*color)[1] += visibility * (render->Ks[1] * render->lights[i].color[1] * R_E_dot_product + amb_color[1] * render->lights[i].color[1] * N_L_dot_product);
		(*color)[2] += visibility * (render->Ks[2] * render->lights[i].color[2] * R_E_dot_product + amb_color[2] * render->lights[i].color[2] * N_L_dot_product);
	}

	(*color)[0] = ((*color)[0] > 1) ? 1 : (*color)[0];
	(*color)[1] = ((*color)[1] > 1) ? 1 : (*color)[1];
	(*color)[2] = ((*color)[2] > 1) ? 1 : (*color)[2];
	return;
}

void Triangle::compute_gouraud_color(double x_im, double y_im, double z_im, GzColor* color, GzRender *render) {
	float nx = float_interpolate(x_im, y_im, normal_plane_x);
	float ny = float_interpolate(x_im, y_im, normal_plane_y);
	float nz = float_interpolate(x_im, y_im, normal_plane_z);
	GzNormalize_vector(&nx, &ny, &nz);
		
	float u = float_interpolate(x_im, y_im, text_plane_u);
	float v = float_interpolate(x_im, y_im, text_plane_v);
		
	// implementing eye vector directed to the focal point
	float d = 1.0 / tan((render->camera.FOV / 2.0) * (PI / 180.0));
	// FIXME: don't know why "ex = -x_im / 2.0;" works beter then "ex = -x_im;" ???
	float ex = -x_im / 2.0;		//0;
	float ey = -y_im / 2.0;		//0;
	float ez = -d -z_im;		//-1;
	GzNormalize_vector(&ex, &ey, &ez);

	(*color)[0] = render->ambientlight.color[0];
	(*color)[1] = render->ambientlight.color[1];
	(*color)[2] = render->ambientlight.color[2];

	for (int i = 0; i < render->numlights; i++) {
		float lx = render->lights[i].position[0];
		float ly = render->lights[i].position[1];
		float lz = render->lights[i].position[2];
		GzNormalize_vector(&lx, &ly, &lz);

		float N_L_dot_product = multipy_vectors(nx, ny, nz, lx, ly, lz);
		float N_E_dot_product = multipy_vectors(nx, ny, nz, ex, ey, ez);
		if ((N_L_dot_product  < 0 && N_E_dot_product > 0) || (N_L_dot_product > 0 && N_E_dot_product < 0))
			continue;
		if (N_L_dot_product  < 0 && N_E_dot_product < 0){
			nx *= -1;
			ny *= -1;
			nz *= -1;
			N_L_dot_product = multipy_vectors(nx, ny, nz, lx, ly, lz);
			N_E_dot_product = multipy_vectors(nx, ny, nz, ex, ey, ez);
		}

		float rx = 2 * nx * N_L_dot_product - lx;
		float ry = 2 * ny * N_L_dot_product - ly;
		float rz = 2 * nz * N_L_dot_product - lz;
		GzNormalize_vector(&rx, &ry, &rz);

		float R_E_dot_product = pow(multipy_vectors(rx, ry, rz, ex, ey, ez), render->spec);
		R_E_dot_product = R_E_dot_product < 0 ? 0 : R_E_dot_product;
			
		(*color)[0] += render->lights[i].color[0] * R_E_dot_product + render->lights[i].color[0] * N_L_dot_product;
		(*color)[1] += render->lights[i].color[1] * R_E_dot_product + render->lights[i].color[1] * N_L_dot_product;
		(*color)[2] += render->lights[i].color[2] * R_E_dot_product + render->lights[i].color[2] * N_L_dot_product;
	}
	return;
}

// public
void Triangle::init_triangle(float p1X, float p1Y, float p1Z, float p2X, float p2Y, float p2Z, float p3X, float p3Y, float p3Z,    
	float p1X_im, float p1Y_im, float p1Z_im, float p2X_im, float p2Y_im, float p2Z_im, float p3X_im, float p3Y_im, float p3Z_im,    
	float n1X, float n1Y, float n1Z, float n2X, float n2Y, float n2Z, float n3X, float n3Y, float n3Z, 
	float p1U, float p1V, float p2U, float p2V, float p3U, float p3V, GzRender* render) {

	b_box = init_bounding_box(p1X, p1Y, p2X, p2Y, p3X, p3Y);
	plane = init_plane(p1X, p1Y, p1Z, p2X, p2Y, p2Z, p3X, p3Y, p3Z);
	plane_im = init_plane(p1X_im, p1Y_im, p1Z_im, p2X_im, p2Y_im, p2Z_im, p3X_im, p3Y_im, p3Z_im);
	init_edges(p1X, p1Y, p1Z, p2X, p2Y, p2Z, p3X, p3Y, p3Z, 
		p1X_im, p1Y_im, p1Z_im, p2X_im, p2Y_im, p2Z_im, p3X_im, p3Y_im, p3Z_im);

	GzNormalize_vector(&n1X, &n1Y, &n1Z);
	GzNormalize_vector(&n2X, &n2Y, &n2Z);
	GzNormalize_vector(&n3X, &n3Y, &n3Z);
		
	normal_plane_x = init_plane(p1X_im, p1Y_im, n1X, p2X_im, p2Y_im, n2X, p3X_im, p3Y_im, n3X);
	normal_plane_y = init_plane(p1X_im, p1Y_im, n1Y, p2X_im, p2Y_im, n2Y, p3X_im, p3Y_im, n3Y);
	normal_plane_z = init_plane(p1X_im, p1Y_im, n1Z, p2X_im, p2Y_im, n2Z, p3X_im, p3Y_im, n3Z);


	text_plane_u = init_plane(p1X_im, p1Y_im, p1U, p2X_im, p2Y_im, p2U, p3X_im, p3Y_im, p3U);
	text_plane_v = init_plane(p1X_im, p1Y_im, p1V, p2X_im, p2Y_im, p2V, p3X_im, p3Y_im, p3V);

	GzColor color[3];
	compute_gouraud_color(p1X_im, p1Y_im, p1Z_im, color + 0, render);
	compute_gouraud_color(p2X_im, p2Y_im, p2Z_im, color + 1, render);
	compute_gouraud_color(p3X_im, p3Y_im, p3Z_im, color + 2, render);
		
	color_plane_r = init_plane(p1X_im, p1Y_im, color[0][0], p2X_im, p2Y_im, color[1][0], p3X_im, p3Y_im, color[2][0]);
	color_plane_g = init_plane(p1X_im, p1Y_im, color[0][1], p2X_im, p2Y_im, color[1][1], p3X_im, p3Y_im, color[2][1]);
	color_plane_b = init_plane(p1X_im, p1Y_im, color[0][2], p2X_im, p2Y_im, color[1][2], p3X_im, p3Y_im, color[2][2]);

	compute_color(p1X_im, p1Y_im, p1Z_im, color + 0, render, false);
	compute_color(p2X_im, p2Y_im, p2Z_im, color + 1, render, false);
	compute_color(p3X_im, p3Y_im, p3Z_im, color + 2, render, false);
	color_plane_r_no_txt = init_plane(p1X_im, p1Y_im, color[0][0], p2X_im, p2Y_im, color[1][0], p3X_im, p3Y_im, color[2][0]);
	color_plane_g_no_txt = init_plane(p1X_im, p1Y_im, color[0][1], p2X_im, p2Y_im, color[1][1], p3X_im, p3Y_im, color[2][1]);
	color_plane_b_no_txt = init_plane(p1X_im, p1Y_im, color[0][2], p2X_im, p2Y_im, color[1][2], p3X_im, p3Y_im, color[2][2]);
}
	
void Triangle::rasterize(GzRender *render) {
	if (is_flat)
		return;
	for (int i = b_box.minY; i <= b_box.maxY; i++) {
		for (int j = b_box.minX; j <= b_box.maxX; j++) {
			int sign = calculate_sign(j, i);

			if (sign == 1) {
				GzDepth z = int_interpolate(j, i, plane);
				float d = 1.0 / tan((render->camera.FOV / 2.0) * (PI / 180.0));
				float z_im = z * d / (INT_MAX - z);
				float x_im = (2.0 * (j - render->display->x_shift) / render->display->xres - 1.0) * (((float)z) / (INT_MAX - (float)z) + 1.0);
				float y_im = (1.0 - 2.0 * (i - render->display->y_shift) / render->display->yres) * (((float)z) / (INT_MAX - (float)z) + 1.0);

				GzIntensity r,g,b;
				if (render->interp_mode == GZ_FLAT) {
					GzColor cl;
					compute_color(x_im, y_im, z_im, &cl, render, true);
					r = ctoi(cl[0]);
					g = ctoi(cl[1]);
					b = ctoi(cl[2]);
				}
				if (render->interp_mode == GZ_COLOR) {
					float u = float_interpolate(x_im, y_im, text_plane_u);
					float v = float_interpolate(x_im, y_im, text_plane_v);
					GzColor color; 
					if (render->tex_fun != NULL){
						int status = render->tex_fun(u, v, color);
						if (status) 
							exit(GZ_FAILURE); 
						color[0] = color[0] * float_interpolate(x_im, y_im, color_plane_r);
						color[1] = color[1] * float_interpolate(x_im, y_im, color_plane_g);
						color[2] = color[2] * float_interpolate(x_im, y_im, color_plane_b);

					} else {
						color[0] = float_interpolate(x_im, y_im, color_plane_r_no_txt);
						color[1] = float_interpolate(x_im, y_im, color_plane_g_no_txt);
						color[2] = float_interpolate(x_im, y_im, color_plane_b_no_txt);
					}
					color[0] = (color[0] > 1) ? 1 : color[0];
					color[1] = (color[1] > 1) ? 1 : color[1];
					color[2] = (color[2] > 1) ? 1 : color[2];

					color[0] = (color[0] < 0) ? 0 : color[0];
					color[1] = (color[1] < 0) ? 0 : color[1];
					color[2] = (color[2] < 0) ? 0 : color[2];


					r = ctoi(color[0]);
					g = ctoi(color[1]);
					b = ctoi(color[2]);
				}
				if (render->interp_mode == GZ_NORMALS) {
					GzColor cl;
					compute_color(x_im, y_im, z_im, &cl, render, false);
					r = ctoi(cl[0]);
					g = ctoi(cl[1]);
					b = ctoi(cl[2]);
				}
					
				GzPutDisplay(render->display, j, i, r, g, b, 1, z);
			}				
		}
	}
	return;
}

