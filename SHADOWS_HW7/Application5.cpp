// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/
#include <vector>
#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "newmodel_.asc"
//#define INFILE  "board.obj"
#define OUTFILE "output.ppm"

#define IMAGE_SIZE  512
#define MAX_NUMBER_OF_TRIANGLES 150000

#define AA_ENABLED
//#undef AA_ENABLED


//#define OBJ_ENABLED

#define NUMBER_OF_LIGHTS 3 // no more then 3!

float   AAFilter[AAKERNEL_SIZE][3] 	= /* each sample is defined by Xshift, Yshift, weight*/
		{  -0.52, 0.38, 0.128,                  0.41, 0.56, 0.119,                     0.27, 0.08, 0.294,
		-0.17, -0.29, 0.249,                    0.58, -0.55, 0.104,                   -0.31, -0.71, 0.106    };
float zero = 0.0;

extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */

void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 
 	m_nWidth = IMAGE_SIZE;		// frame buffer and display width
	m_nHeight = IMAGE_SIZE;    // frame buffer and display height

	status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

	status |= GzNewDisplay(&m_pDisplay, m_nWidth, m_nHeight);
	for (int i = 0; i < AAKERNEL_SIZE; i++) {
		status |= GzNewDisplay(&(m_pDisplays[i]), m_nWidth, m_nHeight);
	}	
	
	status |= GzGetDisplayParams(m_pDisplay, &xRes, &yRes); 

	status |= GzNewRender(&m_pRender, m_pDisplay); 

/* Translation matrix */
GzMatrix	scale = 
{ 
	3.25,	0.0,	0.0,	0.0, 
	0.0,	3.25,	0.0,	-3.25, 
	0.0,	0.0,	3.25,	3.5, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateX = 
{ 
	1.0,	0.0,	0.0,	0.0, 
	0.0,	.7071,	.7071,	0.0, 
	0.0,	-.7071,	.7071,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateY = 
{ 
	.866,	0.0,	-0.5,	0.0, 
	0.0,	1.0,	0.0,	0.0, 
	0.5,	0.0,	.866,	0.0, 
	0.0,	0.0,	0.0,	1.0 
};

GzMatrix rotateX2=
{
	1,0,0,0,
	0,.939,-0.342,0,
	0,.3420,.939,0,
	0,0,0,1
};
GzMatrix rotateZ2=
{
		.939,-0.342,0,0,
	.3420,.939,0,0,
	0,0,1,0,

	0,0,0,1
};
GzMatrix rotateY2=
{
	0,0,-1,0,
	0,1,0,0,
	1,0,0,0,
	0,0,0,1
};
GzMatrix rotateY3=
{
	0.866,0,0.5,0,
	0,1,0,0,
	-0.5,0,0.866,0,
	0,0,0,1
};
GzMatrix Translate2=
{ 1,0,0,1,
0,1,0,0,
0,0,1,0,
0,0,0,1
};

#if 1 	
    camera.position[X] = 0.0;    
  	camera.position[Y] = 10.0;
  	camera.position[Z] = -15.0;

  	camera.lookat[X] = 0.0;
  	camera.lookat[Y] = 0.0;
  	camera.lookat[Z] = 0.0;

  	camera.worldup[X] = 0.0;
  	camera.worldup[Y] = 1.0;
  	camera.worldup[Z] = 0.0;

	camera.FOV = 50.0;  
	status |= GzPutCamera(m_pRender, &camera); 
	
#endif 

	/* Start Renderer */
	status |= GzBeginRender(m_pRender);

	float w;
	/* Light */
	GzLight	light1 = { {10.0, 8.0, 10.0}, {0.0, 0.0, 0.0}, {0.5, 0.5, 0.9}, 5};
	multiplyMatrixByVector(light1.position[0], light1.position[1], light1.position[2], m_pRender->Ximage_im[m_pRender->matlevel], 
		&(light1.position_im[0]), &(light1.position_im[1]), &(light1.position_im[2]), &w);
	light1.position_im[0] /= w;
	light1.position_im[1] /= w;
	light1.position_im[2] /= w;

	GzLight	light2 = { {-10.0, 8.0, 10.0}, {0.0, 0.0, 0.0}, {0.9, 0.2, 0.3}, 5};
	multiplyMatrixByVector(light2.position[0], light2.position[1], light2.position[2], m_pRender->Ximage_im[m_pRender->matlevel], 
		&(light2.position_im[0]), &(light2.position_im[1]), &(light2.position_im[2]), &w);
	light2.position_im[0] /= w;
	light2.position_im[1] /= w;
	light2.position_im[2] /= w;
	
	GzLight	light3 = { {0.0, 8.0, 15.0}, {0.0, 0.0, 0.0}, {0.2, 0.7, 0.3}, 5};
	multiplyMatrixByVector(light3.position[0], light3.position[1], light3.position[2], m_pRender->Ximage_im[m_pRender->matlevel],
		&(light3.position_im[0]), &(light3.position_im[1]), &(light3.position_im[2]), &w);
	light3.position_im[0] /= w;
	light3.position_im[1] /= w;
	light3.position_im[2] /= w;
	
	GzLight	ambientlight = { {0, 0, 0}, {0, 0, 0}, {0.3, 0.3, 0.3}, 1};

	GzBoundingBox	bbox = { -10.0, 10.0, -10.0, 10.0, -10.0, 10.0};

	/* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
	GzColor ambientCoefficient = { 0.2, 0.2, 0.2 };
	GzColor diffuseCoefficient = {0.7, 0.7, 0.7};

/* 
  renderer is ready for frame --- define lights and shader at start of frame 
*/

        /*
         * Tokens associated with light parameters
         */
		nameListLights[0] = GZ_BBOX;
        valueListLights[0] = (GzPointer)&bbox;
        status |= GzPutAttribute(m_pRender, 1, nameListLights, valueListLights);

        nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[0] = (GzPointer)&light1;
		nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[1] = (GzPointer)&light2;
        nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[2] = (GzPointer)&light3;
        status |= GzPutAttribute(m_pRender, NUMBER_OF_LIGHTS, nameListLights, valueListLights);
		
		nameListLights[0] = GZ_AMBIENT_LIGHT;
        valueListLights[0] = (GzPointer)&ambientlight;
        status |= GzPutAttribute(m_pRender, 1, nameListLights, valueListLights);

		nameListLights[0] = GZ_AASHIFTX;
		valueListLights[0] = (GzPointer)&(zero);
		status |= GzPutAttribute(m_pRender,  1, nameListLights, valueListLights);
		nameListLights[0] = GZ_AASHIFTY;
		valueListLights[0] = (GzPointer)&(zero);
		status |= GzPutAttribute(m_pRender,  1, nameListLights, valueListLights);

		for (int i = 0; i < AAKERNEL_SIZE; i++) {
			status |= GzUpdateRender(m_pRender, m_pDisplays[i]);
			
			nameListLights[0] = GZ_AASHIFTX;
			valueListLights[0] = (GzPointer)&(AAFilter[i][0]);
			status |= GzPutAttribute(m_pRender,  1, nameListLights, valueListLights);
	
			nameListLights[0] = GZ_AASHIFTY;
			valueListLights[0] = (GzPointer)&(AAFilter[i][1]);
			status |= GzPutAttribute(m_pRender,  1, nameListLights, valueListLights);
		}

        /*
         * Tokens associated with shading 
         */
        nameListShader[0]  = GZ_DIFFUSE_COEFFICIENT;
        valueListShader[0] = (GzPointer)diffuseCoefficient;

	/* 
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode  
	*/
		nameListShader[1]  = GZ_INTERPOLATE;
		interpStyle = GZ_NORMALS;         /* Phong shading */
        valueListShader[1] = (GzPointer)&interpStyle;

        nameListShader[2]  = GZ_AMBIENT_COEFFICIENT;
        valueListShader[2] = (GzPointer)ambientCoefficient;
        nameListShader[3]  = GZ_SPECULAR_COEFFICIENT;
        valueListShader[3] = (GzPointer)specularCoefficient;
        nameListShader[4]  = GZ_DISTRIBUTION_COEFFICIENT;
        specpower = 32;
        valueListShader[4] = (GzPointer)&specpower;

        nameListShader[5]  = GZ_TEXTURE_MAP;
#if 0   /* set up null texture function or valid pointer */
        valueListShader[5] = (GzPointer)0;
#else
        valueListShader[5] = (GzPointer)(NULL);	/* or use ptex_fun */
#endif
        status |= GzPutAttribute(m_pRender, 6, nameListShader, valueListShader);

		// model space transformations:
		/*
		status |= GzPushMatrix(m_pRender, scale);
		status |= GzPushMatrix(m_pRender, rotateY); 
		status |= GzPushMatrix(m_pRender, rotateX); 
		status|=GzPushMatrix(m_pRender,rotateX2);
		status|=GzPushMatrix(m_pRender,rotateZ2);
		status|=GzPushMatrix(m_pRender,rotateY2);
		status|=GzPushMatrix(m_pRender,rotateY2);
		status|=GzPushMatrix(m_pRender,rotateY3);
		status|=GzPushMatrix(m_pRender,Translate2);
		*/
		if (status) exit(GZ_FAILURE); 

		if (status) 
			return(GZ_FAILURE); 
		else 
			return(GZ_SUCCESS); 
}

int Application5::LoadObjModel(Model** pModel){
	char strLine[255]		= {0};
	char ch					= 0;
	char c=0;
	int i=0,j=0,k=0;
	GzCoord* NewVertex = new GzCoord[MAX_NUMBER_OF_TRIANGLES];
	GzCoord* NewNormal = new GzCoord[MAX_NUMBER_OF_TRIANGLES];
	GzTextureIndex* Newtext = new GzTextureIndex[MAX_NUMBER_OF_TRIANGLES];
	int number_of_triangles=0;
	int v[3];
	int t[3];
	int vn[3];
	GzCoord Vect1,Vect2;
	GzCoord Final;
	int vnpresent=0;

	FILE *output;
	output= fopen(INFILE, "r");
	if(!output) {	
		return -1;
	}

	int obj_count = -1;
   	while(!feof(output)){
		float x = 0.0f, y = 0.0f, z = 0.0f;
		ch = fgetc(output);
		pModel[number_of_triangles] = new Model();
		switch(ch) {
		case 'o':fgets(strLine, 100, output);
				obj_count++;
				break;
		case 'v':c = fgetc(output);
                 if(c==' ') {
                    fscanf(output, "%f %f %f", &NewVertex[i][0], &NewVertex[i][1], &NewVertex[i][2]);
					i++;
				 } else if(c=='t'){
					fscanf(output, "%f %f", &Newtext[j][0], &Newtext[j][1]);
					j++;
				 } else if(c=='n'){
					fscanf(output, "%f %f %f", &NewNormal[k][0], &NewNormal[k][1],&NewNormal[k][2]);
					k++;
				 }
				 c = fgetc(output);
				 break;
		case 'f':c = fgetc(output);
			    if(c==' ') {
					fscanf(output, "%d/%d/%d %d/%d/%d %d/%d/%d",&v[0],&t[0],&vn[0]
														,&v[1],&t[1],&vn[1]
														,&v[2],&t[2],&vn[2]);
					for(int m=0;m<3;m++) {
						(pModel[number_of_triangles])->side[m][0]=NewVertex[v[m]-1][0]/1.0;
						(pModel[number_of_triangles])->side[m][1]=NewVertex[v[m]-1][1]/1.0;
						(pModel[number_of_triangles])->side[m][2]=NewVertex[v[m]-1][2]/1.0;
						(pModel[number_of_triangles])->normal[m][0]=NewNormal[vn[m]-1][0];
						(pModel[number_of_triangles])->normal[m][1]=NewNormal[vn[m]-1][1];
						(pModel[number_of_triangles])->normal[m][2]=NewNormal[vn[m]-1][2];
						(pModel[number_of_triangles])->text[m][0]=Newtext[t[m]-1][0];
						(pModel[number_of_triangles])->text[m][1]=Newtext[t[m]-1][1];
					}													  
					number_of_triangles++;
				}
				c = fgetc(output);
				break;
		default:fgets(strLine, 100, output);
		}
	}
	fclose(output);

	delete NewNormal;
	delete NewVertex;
	delete Newtext;
	return number_of_triangles;
}


int Application5::LoadModel(GzCoord* vertexLists, GzCoord* normalLists, GzTextureIndex* uvLists){
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */ 
	char		dummy[256]; 

	// I/O File open
	FILE *infile;
	if( (infile  = fopen( INFILE , "r" )) == NULL )	{
         AfxMessageBox( "The input file was not opened\n" );
		 return GZ_FAILURE;
	}
	/* 
	* Walk through the list of triangles, set color 
	* and render each triangle 
	*/ 
	int number_of_triangles = 0;
	while( fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
			&(vertexList[0][0]), &(vertexList[0][1]),  
			&(vertexList[0][2]), 
			&(normalList[0][0]), &(normalList[0][1]), 	
			&(normalList[0][2]), 
			&(uvList[0][0]), &(uvList[0][1]) ); 
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
			&(vertexList[1][0]), &(vertexList[1][1]), 	
			&(vertexList[1][2]), 
			&(normalList[1][0]), &(normalList[1][1]), 	
			&(normalList[1][2]), 
			&(uvList[1][0]), &(uvList[1][1]) ); 
			fscanf(infile, "%f %f %f %f %f %f %f %f", 
			&(vertexList[2][0]), &(vertexList[2][1]), 	
			&(vertexList[2][2]), 
			&(normalList[2][0]), &(normalList[2][1]), 	
			&(normalList[2][2]), 
			&(uvList[2][0]), &(uvList[2][1]) ); 
			
			//FIXME: deep copy
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					vertexLists[number_of_triangles + i][j] = vertexList[i][j];
					normalLists[number_of_triangles + i][j] = normalList[i][j]; 
				}	
				uvLists[number_of_triangles + i][0] = uvLists[i][0];
				uvLists[number_of_triangles + i][1] = uvLists[i][1];
			}
			number_of_triangles = number_of_triangles + 3;
	}	
	/* 
	 * Close file
	 */ 
	if( fclose( infile ) )
      AfxMessageBox( "The input file was not closed\n" );
	return (number_of_triangles / 3);
}

int Application5::Render() 
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */ 

	#ifndef OBJ_ENABLED
	GzCoord			vertexLists[MAX_NUMBER_OF_TRIANGLES / 10];	/* vertex position coordinates */ 
	GzCoord			normalLists[MAX_NUMBER_OF_TRIANGLES / 10];	/* vertex normals */ 
	GzTextureIndex  uvLists[MAX_NUMBER_OF_TRIANGLES / 10];		/* vertex texture map indices */ 
	#endif
	int			status; 


	/* Initialize Display */
	status |= GzInitDisplay(m_pDisplay); 
	for (int i = 0; i < AAKERNEL_SIZE; i++)
		status |= GzInitDisplay(m_pDisplays[i]); 
	for (int i = 0; i < m_pRender->numlights; i++) {
		GzRender* map = m_pRender->lights_shadow_maps[i];
		status |= GzInitDisplay(map->display);
	}
	/* 
	* Tokens associated with triangle vertex values 
	*/ 
	nameListTriangle[0] = GZ_POSITION; 
	nameListTriangle[1] = GZ_NORMAL; 
	nameListTriangle[2] = GZ_TEXTURE_INDEX;  


	// PASS I
	#ifdef OBJ_ENABLED
	Model* triangles[MAX_NUMBER_OF_TRIANGLES];
	int number_of_triangles = LoadObjModel(triangles);
	#endif
	#ifndef OBJ_ENABLED
	int number_of_triangles = LoadModel(vertexLists, normalLists, uvLists);
	#endif
	
	for (int k = 0; k < number_of_triangles; k++) {
		#ifdef OBJ_ENABLED
		valueListTriangle[0] = (GzPointer)(triangles[k]->side);   
		valueListTriangle[1] = (GzPointer)(triangles[k]->normal); 
		valueListTriangle[2] = (GzPointer)(triangles[k]->text);  
		#endif
		#ifndef OBJ_ENABLED
		valueListTriangle[0] = (GzPointer)vertexLists[3*k]; 
		valueListTriangle[1] = (GzPointer)normalLists[3*k]; 
		valueListTriangle[2] = (GzPointer)uvLists[3*k]; 
		#endif
		// shadow map rendering
		for (int i = 0; i < m_pRender->numlights; i++) {
			GzRender* map = m_pRender->lights_shadow_maps[i];
			GzPutTriangle(map, 3, nameListTriangle, valueListTriangle); 
		}
	}
	// PASS II
	for (int k = 0; k < number_of_triangles; k++) {
		#ifdef OBJ_ENABLED
		valueListTriangle[0] = (GzPointer)(triangles[k]->side);   
		valueListTriangle[1] = (GzPointer)(triangles[k]->normal); 
		valueListTriangle[2] = (GzPointer)(triangles[k]->text);  
		#endif
		#ifndef OBJ_ENABLED
		valueListTriangle[0] = (GzPointer)vertexLists[3 * k]; 
		valueListTriangle[1] = (GzPointer)normalLists[3 * k]; 
		valueListTriangle[2] = (GzPointer)uvLists[3 * k]; 
		#endif
		status |= GzUpdateRender(m_pRender, m_pDisplay);
		GzPutTriangle(m_pRender, 3, nameListTriangle, valueListTriangle); 
		#ifdef AA_ENABLED
		for (int i = 0; i < AAKERNEL_SIZE; i++) {
			status |= GzUpdateRender(m_pRender, m_pDisplays[i]);
			GzPutTriangle(m_pRender, 3, nameListTriangle, valueListTriangle); 
		}
		#endif
	}

	#ifdef OBJ_ENABLED
	for (int i = 0; i < number_of_triangles; i++)
		delete triangles[i];
	#endif

	#ifdef AA_ENABLED
	for (int l = 0; l < m_pDisplay->xres * m_pDisplay->yres; l++) {
		m_pDisplay->fbuf[l].red   = 0.0;
		m_pDisplay->fbuf[l].green = 0.0;
		m_pDisplay->fbuf[l].blue  = 0.0;
		m_pDisplay->fbuf[l].z     = 0.0;
		for (int i = 0; i < AAKERNEL_SIZE; i++) {
				m_pDisplay->fbuf[l].red   += AAFilter[i][2]  * m_pDisplays[i]->fbuf[l].red  ;
				m_pDisplay->fbuf[l].green += AAFilter[i][2]  * m_pDisplays[i]->fbuf[l].green;
				m_pDisplay->fbuf[l].blue  += AAFilter[i][2]  * m_pDisplays[i]->fbuf[l].blue ;
				m_pDisplay->fbuf[l].z     += AAFilter[i][2]  * m_pDisplays[i]->fbuf[l].z    ;
		}
	}
	#endif
	
	FILE *outfile;
	if( (outfile  = fopen( OUTFILE , "wb" )) == NULL ) {
         AfxMessageBox( "The output file was not opened\n" );
		 return GZ_FAILURE;
	}
	
	GzFlushDisplay2File(outfile, m_pDisplay); 	// write out or update display to file
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, m_pDisplay);	// write out or update display to frame buffer
	
	if( fclose( outfile ) )
	    AfxMessageBox( "The output file was not closed\n" );


	// THE END!!!
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	status |= GzFreeRender(m_pRender); 
	status |= GzFreeDisplay(m_pDisplay);
	for (int i = 0; i < AAKERNEL_SIZE; i++) {
		status |= GzFreeDisplay(m_pDisplays[i]);
	}

	status |= GzFreeTexture();
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}



