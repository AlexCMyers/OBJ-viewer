//
//		          Programming Assignment #2
//
//			      Alex Myers
//				  9-23-19
//		assn1.cpp contains code for most directly opengl related functionality
//      including the main function, display functiion, and functions for
//		recording mouse clicks and key strokes
/***************************************************************************/

                                                 
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#ifdef __APPLE__
#include <OpenGl/glu.h>
#else
#include <GL/glu.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>  
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <algorithm>

#define WIDTH 600
#define HEIGHT 600

/* structs to hold points that define each shape */
typedef struct{
	double x, y, z, w;
} vert;

typedef struct{
	double i, j, k;
} normal;

typedef struct{
	double u, v, w;
} texture;

typedef struct{
	long vertNum, textNum, normNum;
}vertSet;

typedef struct{
	double r, g, b;
}color_t;

typedef struct{
	std::vector<vertSet> verts;
	normal normV;
	color_t c;
}face;

typedef struct{
	normal normV;
	double d;
}plane;

typedef struct{
	normal v;
	color_t La, Ld, Ls;
}light;

typedef struct{
	color_t Ka, Kd, Ks;
	int specCoef;
}surface;

/* Globals to keep track of recent points and the setting */
int x_last,y_last;
char currSetting;

bool render;
bool perspective;
bool color;
bool smooth;
bool bump;
vert center;
vert eyeLoc;
surface material;
double r, g, b;



std::vector<light> myLights;
std::vector<vert> myVerts;
std::vector<vert> perspectiveVerts;
std::vector<normal> myNormals;
std::vector<texture> myTextures;
std::vector<face> myFaces;
std::vector<normal> myVertNormals;


/* function prototypes */

void parse(std::ifstream& infile);
void normalizeVerts();
void drawFaces();
void drawPerspectiveFaces();
void drawColorFaces();
void drawColorPerspectiveFaces();
void drawColorQuad(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4);
void drawColorTriangle(int x1, int y1, int x2, int y2, int x3, int y3);

void updatePerspectiveVerts();
void storeVertNormals();
void drawALine(int x1, int y1, int x2, int y2);

void transform(char setting, char direction);
void translate(char direction);
void rotate(char direction);
void scale (char direction);

void vertMult(double matrix[4][4]);

void matrixMult(double mat1[4][4], double mat2[4][4], double result[4][4]);

void write_pixel(int x, int y, double intensity);

bool faceSort(face& left, face& right);

normal calcNorm(face triangle);
normal calcSmoothNorm(int x, int y, double z, face currTri);
normal crossProduct(vert a, vert b, vert c);
normal unitVec(normal v);
normal bumpAdjust(normal oldNorm, face tri, double x, double y, double z);
double dotProd(normal a, normal b);
plane calcPlane(normal n, vert v);
bool insideTriangle(int x, int y, face tri);
bool insideTriangle(int x, int y, vert v1, vert v2, vert v3);
double getZ(int x, int y, plane p);
void renderImage();
void rayTrace(int x, int y, double z, face currTri);

double dist(double x, double y, double z, vert v);
bool sameLine(vert v1, vert v2, vert v3);
normal weightVec(double w1, normal v1, double w2, normal v2);
double magnitude(normal vec);




/* Given by Zordan; used to initialize window*/
/***************************************************************************/

void init_window()
                 /* Clear the image area, and set up the coordinate system */
{

        					       /* Clear the window */
    glClearColor(0.0,0.0,0.0,0.0);
	glShadeModel(GL_SMOOTH);
    glOrtho(0,WIDTH,0,HEIGHT,-1.0,1.0);
}

/* Function that turns pixel at (x,y) "on" with the specified intensity */
/***************************************************************************/

void write_pixel(int x, int y, double intensity)
                                         
{		

		glColor3f (intensity * r, intensity * g, intensity * b);  
        glBegin(GL_POINTS);
		glVertex3i(300 + x, 300 + y, 0);
        glEnd();	
}

/* Function that continuously runs and updates the window */
//***************************************************************************/

void display ( void )  
{
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen 

	//determines if enough points have been selected to draw a new shape
	// based on the current shape setting
    
	
	//draw functions
	if(render){
			renderImage();
	}else if(!perspective){
		if(color){
			drawColorFaces();
		}else{
			drawFaces();
		}
	}else{
		if(color){
			drawColorPerspectiveFaces();
		}else{
			drawPerspectiveFaces();
		}
	}


	glutSwapBuffers(); 
                                    // Draw Frame Buffer 
}


void renderImage(){

	bool hit = false;
	double maxZ;
	face closestTri;

	//for each pixel
	for(int x = -300; x < 300; x++){
		for(int y = -300; y < 300; y++){
			hit = false;
			//for each triangle
			for(int t = 0; t < (int)myFaces.size(); t++){
				//if that pixel is inside that triangle
				if(insideTriangle(x, y, myFaces.at(t))){
					//if this is first relevant triangle found
					if(!hit){
						hit = true;
						//set maxZ to keep track if any closer triangles are found in the future
						if(!perspective){
							maxZ = getZ(x, y, calcPlane(calcNorm(myFaces.at(t)), 
												myVerts.at(myFaces.at(t).verts.at(0).vertNum - 1)));
						}else{
							
							maxZ = getZ(x, y, calcPlane(calcNorm(myFaces.at(t)), 
												perspectiveVerts.at(myFaces.at(t).verts.at(0).vertNum - 1)));
						}
						//set closest triangle to the current triangle
						closestTri = myFaces.at(t);
					}else{
						//if not the first triangle, get z value and compare
						double testZ;
						if(!perspective){
							testZ = getZ(x, y, calcPlane(calcNorm(myFaces.at(t)), 
												myVerts.at(myFaces.at(t).verts.at(0).vertNum - 1)));
						}else{
							testZ = getZ(x, y, calcPlane(calcNorm(myFaces.at(t)), 
												perspectiveVerts.at(myFaces.at(t).verts.at(0).vertNum - 1)));
						}
						if(testZ > maxZ){
							//update values
							maxZ = testZ;
							closestTri = myFaces.at(t);
						}
					}
				}//end if inside triangle
			}//end t loop
			//if we hit any of the triangles, ray trace it to get color
			if(hit){
				r = closestTri.c.r;
				g = closestTri.c.g;
				b = closestTri.c.b;
				rayTrace(x, y, maxZ, closestTri);
			}else{
				//otherwise color black (background color)
				r = g = b = 0;
				write_pixel(x, y, 1.0);
			}
		}
	}
}
normal bumpAdjust(normal oldNorm, face tri, double x, double y, double z){
	vert triVerts[3];
	triVerts[0] = myVerts.at(tri.verts.at(0).vertNum - 1);
	triVerts[1] = myVerts.at(tri.verts.at(1).vertNum - 1);
	triVerts[2] = myVerts.at(tri.verts.at(2).vertNum - 1);

	double d0 = dist(x, y, z, triVerts[0]);
	double d1 = dist(x, y, z, triVerts[1]);
	double d2 = dist(x, y, z, triVerts[2]);
	srand((unsigned int)(d0 + d1 + d2));

	int fac1 = rand()%2 == 1 ? 1 : -1;
	int fac2 = rand()%2 == 1 ? 1 : -1;
	int fac3 = rand()%2 == 1 ? 1 : -1;

	normal bumpNorm;
	bumpNorm.i = oldNorm.i + (fac1 * 1.0/((double)(rand() % 6 + 20)));
	bumpNorm.j = oldNorm.j + (fac2 * 1.0/((double)(rand() % 6 + 20)));
	bumpNorm.k = oldNorm.k + (fac3 * 1.0/((double)(rand() % 6 + 20)));

	return bumpNorm;

}
/***************************************************************************/
void rayTrace(int x, int y, double z, face currTri){
	double lightR = 0.0;
	double lightG = 0.0;
	double lightB = 0.0;

	normal eyeVec = {eyeLoc.x - x, eyeLoc.y - y, eyeLoc.z - z};
	normal unitNorm;
	if(!smooth){
		unitNorm = unitVec(calcNorm(currTri));
	}else{
		unitNorm = unitVec(calcSmoothNorm(x, y, z, currTri));
	}

	if(bump){
		unitNorm = unitVec(bumpAdjust(unitNorm, currTri, x, y, z));
	}
	//for each light source, find vector from point x, y, z that goes to that light
	//comput dot product of the faces normal vector and that light vector (use unit vecs).
	for(auto l : myLights){
		normal lightV = {l.v.i - x, l.v.j - y, l.v.k - z};
		double cosAngleLN = dotProd(unitVec(lightV), unitNorm);
		//if the angle is > 90 degrees, then the light doesnt hit this surface
		//and only ambient light is accounted for
		double angleLN = acos(cosAngleLN);
		//if light hits point, add phong model for direct and spec to RGB
		if(angleLN <= 1.5708){
			//direct light
			//std::cout << lightR << " " << lightG << " " << lightB << std::endl;
			lightR += l.Ld.r * material.Kd.r * dotProd(unitVec(lightV), unitNorm);
			lightG += l.Ld.g * material.Kd.g * dotProd(unitVec(lightV), unitNorm);
			lightB += l.Ld.b * material.Kd.b * dotProd(unitVec(lightV), unitNorm);
			//std::cout << unitVec(lightV).i << " " << unitVec(lightV).j << " " << unitVec(lightV).k << std::endl;
			//std::cout << unitNorm.i << " " << unitNorm.j << " " << unitNorm.k << std::endl;
		//	std::cout << lightR << " " << lightG << " " << lightB << std::endl;
			//finding reflection vector
			//r = lightV - 2(dotProd(lightV, faceNorm) x n);
			// and cosAngleLN is that dot product already
			normal refl = {lightV.i - (2 * cosAngleLN * unitNorm.i), 
							lightV.j - (2 * cosAngleLN * unitNorm.j),
							lightV.k - (2 * cosAngleLN * unitNorm.k)};
			//spec light
			lightR += l.Ls.r * material.Ks.r * pow(dotProd(unitVec(eyeVec), unitVec(refl)), material.specCoef);
			lightG += l.Ls.g * material.Ks.g * pow(dotProd(unitVec(eyeVec), unitVec(refl)), material.specCoef);
			lightB += l.Ls.b * material.Ks.b * pow(dotProd(unitVec(eyeVec), unitVec(refl)), material.specCoef);

		}
		lightR += l.La.r * material.Ka.r;
		lightG += l.La.g * material.Ka.g;
		lightB += l.La.b * material.Ka.b;
	}
	r = lightR;
	g = lightG;
	b = lightB;

	write_pixel(x, y, 1.0);
}

/* Given by Dr. Zordan; altered to record points */
/***************************************************************************/
void mouse(int button, int state, int x, int y)
{
	
}
 
/* Alters the current Setting based on key input */
/***************************************************************************/
void keyboard ( unsigned char key, int x, int y )  
{

	switch ( key ) {
		case 'v':
			//perspective = !perspective;
			break;
		case 'r':
		case 't':
		case 'e':
			currSetting = key;
			break;
		case 'a':
		case 'w':
		case 's':
		case 'd':
			transform(currSetting, key);
			break;
		case 'z':
			//color = !color;
			//if(!color){
			//	r = b = g = 1.0;
			//}
			break;
		case 'x':
			smooth = !smooth;
			break;
		case 'b':
			bump = !bump;
			break;
		case 27:              // When Escape Is Pressed...
			exit ( 0 );   // Exit The Program
			break;        
	    case '1':            
		    printf("New screen\n");
			break;
		default:       
			break;
	}
}
/***************************************************************************/
void parse(std::ifstream& infile){
	std::string currLine;
	std::string signal;
	srand(time(NULL));
	//while there are still lines in the file
	while(getline(infile, currLine)){
		std::istringstream line(currLine);
		line >> signal;
		//vertices
		if(signal == "v" && currLine != ""){
			double x = 0.0, y = 0.0, z = 0.0, w = 0.0;
			line >> x >> y >> z;
			if(line >> w){
				myVerts.push_back({x, y, z, w});
				perspectiveVerts.push_back({x, y, z, w});
			}else{
				myVerts.push_back({x, y, z, 1.0});
				perspectiveVerts.push_back({x, y, z, 1.0});
			}
		}
		//textures
		else if(signal == "vt" && currLine != ""){
			double u = 0.0, v = 0.0, w = 0.0;
			line >> u;
			if(line >> v){
				if(line >> w){
					myTextures.push_back({u, v, w});
				}else{
					myTextures.push_back({u, v, 0.0});
				}
			}else{
				myTextures.push_back({u, 0.0, 0.0});
			}
		}
		//normals
		else if(signal == "vn" && currLine != ""){
			double i = 0.0, j = 0.0, k = 0.0;
			line >> i >> j >> k;
			myNormals.push_back({i, j, k});

		}
		//faces
		else if(signal == "f" && currLine != ""){
			face newFace;
			vertSet newVertSet;
			//note: max length for long is 18 digits
			//vertCombo can hold 3 longs and 2 '/' chars
			// 1 2 3 4
			// 1/3 2/2 
			// f 1//4 5//6 7//4
			// 355737
			char vertCombo[(18 * 3) + 2];
			//while still verts for the face (should be 3 or 4)
			while(line.getline(vertCombo, (18 * 3) + 2,' ')){
				std::string sVertCombo(vertCombo);
				if(sVertCombo != "" && sVertCombo != " "){
					std::stringstream ssVertCombo(sVertCombo);
					char nextVal[18];
					ssVertCombo.getline(nextVal, 18, '/');
					std::string sVertNum(nextVal);
					newVertSet.vertNum = std::stol(sVertNum, nullptr);
					if(ssVertCombo.getline(nextVal, 18, '/')){
						if(std::string(nextVal) != ""){
							std::string sTextNum(nextVal);
							newVertSet.textNum = std::stol(nextVal, nullptr);
						}else{
							newVertSet.textNum = 0;
						}
						if(ssVertCombo.getline(nextVal, 18, ' ')){
							std::string sNormNum(nextVal);
							newVertSet.normNum = std::stol(sNormNum, nullptr);
						}
						else{
							newVertSet.normNum = 0;
						}
					}else{
						newVertSet.normNum = 0;
						newVertSet.textNum = 0;
					}
					
					newFace.verts.push_back(newVertSet);
				}
			}
			newFace.c.r = double((rand() % 1000))/(double)1000;
			newFace.c.g = double((rand() % 1000))/(double)1000;
			newFace.c.b = double((rand() % 1000))/(double)1000;
			//at this point, the new face may not be a triangle
			//need to divide into triangles
			std::vector<face> newTris;
			for(int t = 2; t < (int)newFace.verts.size(); t++ ){
				face newTri;
				newTri.verts.push_back(newFace.verts.at(0));
				newTri.verts.push_back(newFace.verts.at(t-1));
				newTri.verts.push_back(newFace.verts.at(t));
				newTri.c = newFace.c;
				newTris.push_back(newTri);
			}

			normal normVec = calcNorm(newTris.at(0));

			//loading triangles into the face vector
			for(int t = 0; t < (int)newTris.size(); t++){
				newTris.at(t).normV = normVec;
				myFaces.push_back(newTris.at(t));
			}
			
		}
		
	}
}
/***************************************************************************/
void storeVertNormals(){
	//for each vert, check if it's part of each face
	//if it is, add the components of the normal of that face to the running totals
	//and then store the average in vertNormals
	for(int i = 0; i < (int)myVerts.size(); i++){
		normal newNormal;
		double total = 0;
		//check each face
		for(auto t: myFaces){
			//check each vertNum for each face
			for(auto v: t.verts){
				//if the vertNum is the current vert
				if(v.vertNum == i+1){
					newNormal.i += t.normV.i;
					newNormal.j += t.normV.j;
					newNormal.k += t.normV.k;
					total++;
				}
			}
		}//done checking all faces
		newNormal.i /= total;
		newNormal.j /= total;
		newNormal.k /= total;
		myVertNormals.push_back(newNormal);
	}
}
/***************************************************************************/
//utilizes the cross product function to get the normal of the plane
normal calcNorm(face triangle){
	vert v1, v2, v3;
	normal fNorm;

	if(!perspective){
		v1 = myVerts.at(triangle.verts.at(0).vertNum - 1);
		v2 = myVerts.at(triangle.verts.at(1).vertNum - 1);
		v3 = myVerts.at(triangle.verts.at(2).vertNum - 1);
	}else{
		v1 = perspectiveVerts.at(triangle.verts.at(0).vertNum - 1);
		v2 = perspectiveVerts.at(triangle.verts.at(1).vertNum - 1);
		v3 = perspectiveVerts.at(triangle.verts.at(2).vertNum - 1);
	}
	fNorm = crossProduct(v1,v2,v3);
	return fNorm;
}


normal calcSmoothNorm(int x, int y, double z, face currTri){
	vert currVert = {(double)x, (double)y, z, 0};
	vert triVerts[3];
	triVerts[0] = myVerts.at(currTri.verts.at(0).vertNum - 1);
	triVerts[1] = myVerts.at(currTri.verts.at(1).vertNum - 1);
	triVerts[2] = myVerts.at(currTri.verts.at(2).vertNum - 1);

	normal currNormals[3];
	currNormals[0] = myVertNormals.at(currTri.verts.at(0).vertNum - 1);
	currNormals[1] = myVertNormals.at(currTri.verts.at(1).vertNum - 1);
	currNormals[2] = myVertNormals.at(currTri.verts.at(2).vertNum - 1);

	double area012 = magnitude(crossProduct(triVerts[0], triVerts[1], triVerts[2]));
	double area0P2 = magnitude(crossProduct(triVerts[0], currVert, triVerts[2]));
	double areaP12 = magnitude(crossProduct(currVert, triVerts[1], triVerts[2]));

	double propV0 = areaP12/area012;
	double propV1 = area0P2/area012;
	double propV2 = 1 - propV0 - propV1;

	normal newNorm;
	newNorm.i = (propV0 * currNormals[0].i) + 
				(propV1 * currNormals[1].i) + 
				(propV2 * currNormals[2].i);
	newNorm.j = (propV0 * currNormals[0].j) + 
				(propV1 * currNormals[1].j) + 
				(propV2 * currNormals[2].j);
	newNorm.k = (propV0 * currNormals[0].k) + 
				(propV1 * currNormals[1].k) + 
				(propV2 * currNormals[2].k);

	return newNorm;

}
double magnitude(normal vec){
	return pow(pow(vec.i, 2.0) + pow(vec.j, 2.0) + pow(vec.k, 2.0),.5);
}
bool sameLine(vert v1, vert v2, vert v3){
	normal checkNorm = crossProduct(v1, v2, v3);
	if((int)checkNorm.i == 0 && (int)checkNorm.j == 0 && (int)checkNorm.k == 0){
		return true;
	}
	else{
		return false;
	}
}
normal weightVec(double w1, normal v1, double w2, normal v2){
	return{(w1 * v1.i) + (w2 * v2.i), (w1 * v1.j) + (w2 * v2.j), (w1 * v1.k) + (w2 * v2.k)};
}

double dist(double x, double y, double z, vert v){
	return pow(pow(x - v.x, 2.0) + pow(y - v.y, 2.0) + pow(z - v.z, 2.0), .5);
}
/***************************************************************************/

//finds cross product of 2 vectors defined by vertices of a triangle
normal crossProduct(vert a, vert b, vert c){
	double ABi = b.x - a.x;
	double ABj = b.y - a.y;
	double ABk = b.z - a.z;
	double ACi = c.x - a.x;
	double ACj = c.y - a.y;
	double ACk = c.z - a.z;

	normal result;
	//cross product formula
	result.i = ((ABj * ACk) - (ACj * ABk));
	result.j = ((ACi * ABk) - (ABi * ACk));
	result.k = ((ABi * ACj) - (ABj * ACi));

	return result;
}
/***************************************************************************/
/* uses a normal vector to a plane and a point on that plane to complete plane equation */
plane calcPlane(normal n, vert v){
	plane p;
	p.normV = n;
	//need to find d, recall plane equation where ijk are defined by the normal vector:
	//ix + jy + kz + d = 0
	//ix + jy + kz = -d

	double planeD = (n.i * v.x) + (n.j * v.y) + (n.k * v.z);
	planeD *= -1;

	p.d = planeD;
	return p;
	
}
/***************************************************************************/
/* given x and y, finds z value for point on a plane */
double getZ(int x, int y, plane p){
	//recall plane equation where ijk are defined by the normal vector:
	//ix + jy + kz + d = 0
	//need to use x and y to find z
	//-kz = ix + jy + d
	//z = (ix + jy + d)/-k

	double z = ((p.normV.i * x) + (p.normV.j * y) + p.d)/(-1 * p.normV.k);

	return z;
}
/***************************************************************************/
bool insideTriangle(int x, int y, face tri){
	double x1, x2, x3, y1, y2, y3;
	if(!perspective){
		x1 = myVerts.at(tri.verts.at(0).vertNum - 1).x;
		x2 = myVerts.at(tri.verts.at(1).vertNum - 1).x;
		x3 = myVerts.at(tri.verts.at(2).vertNum - 1).x;
		y1 = myVerts.at(tri.verts.at(0).vertNum - 1).y;
		y2 = myVerts.at(tri.verts.at(1).vertNum - 1).y;
		y3 = myVerts.at(tri.verts.at(2).vertNum - 1).y;
	}else{
		x1 = perspectiveVerts.at(tri.verts.at(0).vertNum - 1).x;
		x2 = perspectiveVerts.at(tri.verts.at(1).vertNum - 1).x;
		x3 = perspectiveVerts.at(tri.verts.at(2).vertNum - 1).x;
		y1 = perspectiveVerts.at(tri.verts.at(0).vertNum - 1).y;
		y2 = perspectiveVerts.at(tri.verts.at(1).vertNum - 1).y;
		y3 = perspectiveVerts.at(tri.verts.at(2).vertNum - 1).y;
	}

			
	double denomenator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
   	double a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denomenator;
   	double b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denomenator;
   	double c = 1 - a - b;
	if((a>=0)&&(a<=1)&&(b>=0)&&(b<=1)&&(c>=0)&&(c<=1)){
		return true;
	}else{
		return false;
	}
}
/***************************************************************************/
bool insideTriangle(int x, int y, vert v1, vert v2, vert v3){
	double x1, x2, x3, y1, y2, y3;

	x1 = v1.x;
	x2 = v2.x;
	x3 = v3.x;
	y1 = v1.y;
	y2 = v2.y;
	y3 = v3.y;
			
	double denomenator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
   	double a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denomenator;
   	double b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denomenator;
   	double c = 1 - a - b;
	if((a>=0)&&(a<=1)&&(b>=0)&&(b<=1)&&(c>=0)&&(c<=1)){
		return true;
	}else{
		return false;
	}
}
normal unitVec(normal v){
	double length = pow(((v.i * v.i) + (v.j * v.j) + (v.k * v.k)), .5);
	normal unit = {v.i/length, v.j/length, v.k/length};
	return unit;
}
/***************************************************************************/
double dotProd(normal a, normal b){
	return((a.i * b.i) + (a.j * b.j) + (a.k * b.k));
}
/***************************************************************************/
void normalizeVerts(){
	//need to find maxes and mins to center the object at 0,0,0
	double xMax, xMin;
	xMax = xMin = myVerts.at(0).x;

	double yMax, yMin;
	yMax = yMin = myVerts.at(0).y;

	double zMax, zMin;
	zMax = zMin = myVerts.at(0).z;
	for(auto v: myVerts){
		//checking mins and maxes for centering purposes
		if(v.x < xMin){xMin = v.x;}
		if(v.x > xMax){xMax = v.x;}
		if(v.y < yMin){yMin = v.y;}
		if(v.y > yMax){yMax = v.y;}
		if(v.z < zMin){zMin = v.z;}
		if(v.z > zMax){zMax = v.z;}
	}
	double xDiff = (xMin + xMax)/2.0;
	double yDiff = (yMin + yMax)/2.0;
	double zDiff = (zMin + zMax)/2.0;
	for(auto &v: myVerts){
		v.x -= xDiff;
		v.y -= yDiff;
		v.z -= zDiff;
	}
	//window is 600 x 600. Want to set it up so everything fits inside a 300 x 300 x 300 cube.
	//which means scaling points so that the absolute max value is 150 (to fit in window nicely)
	//recall: scaling is just multiplying points. If all three coordinates are scaled
	//for every point, then the entire object will scale equally in all directions

	//getting new maxes now that object has been centered
	xMax = xMin = myVerts.at(0).x;
	yMax = yMin = myVerts.at(0).y;
	zMax = zMin = myVerts.at(0).z;
	for(auto v: myVerts){
		if(v.x > xMax){xMax = v.x;}
		if(v.y > yMax){yMax = v.y;}
		if(v.z > zMax){zMax = v.z;}
	}
	double maxXY = (xMax > yMax) ? xMax : yMax;
	double maxVal = (maxXY > zMax) ? maxXY : zMax;
	double factor = 150/maxVal;

	//scaling all verts
	for(auto& v: myVerts){
		v.x *= factor;
		v.y *= factor;
		v.z *= factor;

		//std::cout << "v " << v.x << " " << v.y << " " << v.z << std::endl;
	}
	center.x = 0.0;
	center.y = 0.0;
	center.z = 0.0;
	updatePerspectiveVerts();

}
/***************************************************************************/
void updatePerspectiveVerts(){
	//after initial normalization, the z values range from -150 to 150
	//let's assume a vanishing point of -800 in the z direction.
	//everything along the plane z = 0 will keep it's x and y values.
	//-800 - z/-800 should then be the scale factor for x and y, because where z = 0,
	//that scale factor is 1.  Where z = 100, that scale factor is > 1, where z = -800,
	//that factor is 0, and so on. 
	for(int i = 0; i < (int)myVerts.size(); i++){
		double scaleFactor = (-800.0 - myVerts.at(i).z)/-800.0;
		perspectiveVerts.at(i).x = myVerts.at(i).x * scaleFactor;
		perspectiveVerts.at(i).y = myVerts.at(i).y * scaleFactor;
	}
}
/***************************************************************************/

bool faceSort(face& left, face& right){
	double leftAvg = 0.0, rightAvg = 0.0;
	for(auto& v: left.verts){
		leftAvg += myVerts.at(v.vertNum - 1).z;
	}
	for(auto& v: right.verts){
		rightAvg += myVerts.at(v.vertNum - 1).z;
	}
	leftAvg /= left.verts.size();
	rightAvg /= right.verts.size();

	return leftAvg < rightAvg;
}
void transform(char setting, char direction){
	switch(setting){
		case 't':
			translate(direction);
			break;
		case 'r':
			rotate(direction);
			break;
		case 'e':
			scale(direction);
			break;
		default:
			break;
	}
	updatePerspectiveVerts();
}

/***************************************************************************/
void rotate(char direction){
	double rotationMatrix[4][4];
	double originMatrix[4][4];
	double resetMatrix[4][4];
	double resultMatrix[4][4];
	double finalMatrix[4][4];
	double degree;
	char axis;
	switch(direction){
		case 'a':
			axis = 'y';
			degree = 10.0 * M_PI/180;
			break;
		case 'd':
			axis = 'y';
			degree = -1 * 10.0 * M_PI/180;
			break;
		case 'w':
			axis = 'x';
			degree = 10.0 * M_PI/180;
			break;
		case 's':
			axis = 'x';
			degree = -1 * 10.0 * M_PI/180;
			break;
		default:
			axis = '?';
			break;
	}

	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			rotationMatrix[i][j] = (i == j)? 1.0 : 0.0;
		}
	}
	if(axis == 'y'){
		rotationMatrix[0][0] = cos(degree);
		rotationMatrix[0][2] = -sin(degree);
		rotationMatrix[2][0] = sin(degree);
		rotationMatrix[2][2] = cos(degree);
	}else if(axis == 'x'){
		rotationMatrix[1][1] = cos(degree);
		rotationMatrix[2][1] = -sin(degree);
		rotationMatrix[1][2] = sin(degree);
		rotationMatrix[2][2] = cos(degree);
	}

	for(int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			originMatrix[i][j] = (i == j)? 1.0 : 0.0;
			resetMatrix[i][j] = (i == j)? 1.0 : 0.0;
		}
	}
	originMatrix[0][3] = center.x;
	originMatrix[1][3] = center.y;
	resetMatrix[0][3] = -1 * center.x;
	resetMatrix[1][3] = -1 * center.y;

	matrixMult(originMatrix, rotationMatrix, resultMatrix);
	
	matrixMult(resultMatrix, resetMatrix, finalMatrix);

	vertMult(finalMatrix);
}
/***************************************************************************/
void translate(char direction){
	double translationMatrix[4][4];
	int dx = 0, dy = 0, dz = 0;
	switch(direction){
		case 'w':
			dy = 10.0;
			center.y += 10.0;
			break;
		case 's':
			dy = -10.0;
			center.y -= 10.0;
			break;
		case 'a':
			dx = -10.0;
			center.x -= 10.0;
			break;
		case 'd':
			dx = 10.0;
			center.x += 10.0;
			break;
		default:
			break;

	}
	for(int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			translationMatrix[i][j] = (i == j)? 1.0 : 0.0;
		}
	}
	translationMatrix[0][3] = dx;
	translationMatrix[1][3] = dy;
	translationMatrix[2][3] = dz;
	translationMatrix[3][3] = 1;

	vertMult(translationMatrix);


}
/***************************************************************************/
void scale(char direction){
	double scaleMatrix[4][4];
	double originMatrix[4][4];
	double resetMatrix[4][4];
	double resultMatrix[4][4];
	double finalMatrix[4][4];
	double scaleFactor = 1.0;
	switch(direction){
		case 'a':
		case 's':
			scaleFactor = .9;
			break;
		case 'w':
		case 'd':
			scaleFactor = 1.1;
			break;
		default:
			break;

	}
	for(int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			scaleMatrix[i][j] = (i == j)? scaleFactor : 0.0;
		}
	}
	scaleMatrix[3][3] = 1.0;

	for(int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			originMatrix[i][j] = (i == j)? 1.0 : 0.0;
			resetMatrix[i][j] = (i == j)? 1.0 : 0.0;
		}
	}
	originMatrix[0][3] = center.x;
	originMatrix[1][3] = center.y;
	resetMatrix[0][3] = -1 * center.x;
	resetMatrix[1][3] = -1 * center.y;

	matrixMult(originMatrix, scaleMatrix, resultMatrix);
	matrixMult(resultMatrix, resetMatrix, finalMatrix);

	vertMult(finalMatrix);

}
/***************************************************************************/
void matrixMult(double mat1[4][4], double mat2[4][4], double result[4][4]){
	double dotProd = 0.0;
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			dotProd = 0.0;
			for(int k = 0; k < 4; k++){
				dotProd += mat1[i][k] * mat2[k][j];
			}
			result[i][j] = dotProd;
		}
	}
}
/***************************************************************************/
void vertMult(double matrix[4][4]){
	double newX, newY, newZ;
	for(auto& v: myVerts){
		newX = (matrix[0][0] * v.x) + (matrix[0][1] * v.y) +
			   (matrix[0][2] * v.z) + (matrix[0][3] * 1);
		newY = (matrix[1][0] * v.x) + (matrix[1][1] * v.y) +
			   (matrix[1][2] * v.z) + (matrix[1][3] * 1);
		newZ = (matrix[2][0] * v.x) + (matrix[2][1] * v.y) +
			   (matrix[2][2] * v.z) + (matrix[2][3] * 1);
		v.x = newX;
		v.y = newY;
		v.z = newZ;
	}

}
/***************************************************************************/
int main (int argc, char *argv[])
{
/* This main function sets up the main loop of the program and continues the
   loop until the end of the data is reached.  Then the window can be closed
   using the escape key.						  */

	if(argc != 2){
		std::cout << "Error, expected use is <executable> <filename.obj>" << std::endl;
		return 1;
	}
	std::ifstream input(argv[1]);
	if(!input.is_open()){
		std::cout << "Error opening " << argv[1] << std::endl << "Terminating program" << std::endl;
		return 2;
	}
	perspective = false;
	render = true;
	color = false;
	bump = true;
	smooth = true;
	currSetting = 't';
	r = 1.0;
	g = 1.0;
	b = 1.0;

	eyeLoc.x = 0;
	eyeLoc.y = 0;
	eyeLoc.z = 800;
	//light: location, then 3 colors for ambient, direct, and specular
	light l1 = {{-300, 200, -200}, {.04, .04, .2}, {.7, .7, .9}, {.5, .5, .85}};
	light l2 = {{300, -300, 200}, {.2, .2, .2}, {.5, .5, .5}, {.65, .65, .65}};
	light l3 = {{100, 0, 500}, {.01, .2, .01}, {.4, .9, .6}, {.45, .8, .2}};

	myLights.push_back(l1);
	myLights.push_back(l2);
	myLights.push_back(l3);

	material.Ka = {.25, .25, .25};
	material.Ks = {.4, .4, .4};
	material.Kd = {.9, .9, .9};
	material.specCoef = 150;

	parse(input);
	normalizeVerts();
	storeVertNormals();


	glutInit            ( &argc, argv ); 
    glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
	glutInitWindowSize  ( 600,600 ); 
	glutCreateWindow    ( "Computer Graphics" ); 
	glutDisplayFunc     ( display );  
	glutIdleFunc	    ( display );
	glutMouseFunc       ( mouse );
	glutKeyboardFunc    ( keyboard );
        					      
    init_window();				             //create_window
						       		
	glutMainLoop        ( );                 // Initialize The Main Loop
}

/***************************************************************************/
void drawFaces(){
	r = g = b = 1.0;
	for(auto f: myFaces){
		//if the face is a triangle
		if(f.verts.size() == 3){
			drawALine(int(myVerts.at(f.verts.at(0).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(0).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(1).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(1).vertNum - 1).y));
			drawALine(int(myVerts.at(f.verts.at(1).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(1).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(2).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(2).vertNum - 1).y));
			drawALine(int(myVerts.at(f.verts.at(2).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(2).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(0).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(0).vertNum - 1).y));
		}
		//if it's a quad (note: should only be tris or quads)
		else if(f.verts.size() == 4){
			drawALine(int(myVerts.at(f.verts.at(0).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(0).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(1).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(1).vertNum - 1).y));
			drawALine(int(myVerts.at(f.verts.at(1).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(1).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(2).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(2).vertNum - 1).y));
			drawALine(int(myVerts.at(f.verts.at(2).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(2).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(3).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(3).vertNum - 1).y));
			drawALine(int(myVerts.at(f.verts.at(3).vertNum - 1).x), 
					  int(myVerts.at(f.verts.at(3).vertNum - 1).y),
					  int(myVerts.at(f.verts.at(0).vertNum - 1).x),
					  int(myVerts.at(f.verts.at(0).vertNum - 1).y));
		}
	}
}
/***************************************************************************/
void drawPerspectiveFaces(){
	r = g = b = 1.0;
	for(auto f: myFaces){
		//if the face is a triangle
		if(f.verts.size() == 3){
			drawALine(int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).y));
			drawALine(int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).y));
			drawALine(int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).y));
		}
		//if it's a quad (note: should only be tris or quads)
		else if(f.verts.size() == 4){
			drawALine(int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).y));
			drawALine(int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).y));
			drawALine(int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(3).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(3).vertNum - 1).y));
			drawALine(int(perspectiveVerts.at(f.verts.at(3).vertNum - 1).x), 
					  int(perspectiveVerts.at(f.verts.at(3).vertNum - 1).y),
					  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).x),
					  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).y));
		}
	}
}
/***************************************************************************/
void drawColorFaces(){
	std::sort(myFaces.begin(), myFaces.end(), faceSort);

	for(auto& f: myFaces){
		r = f.c.r;
		b = f.c.b;
		g = f.c.g;
		if(f.verts.size() == 3){
			drawColorTriangle(int(myVerts.at(f.verts.at(0).vertNum - 1).x),
							  int(myVerts.at(f.verts.at(0).vertNum - 1).y),
							  int(myVerts.at(f.verts.at(1).vertNum - 1).x),
							  int(myVerts.at(f.verts.at(1).vertNum - 1).y),
							  int(myVerts.at(f.verts.at(2).vertNum - 1).x),
							  int(myVerts.at(f.verts.at(2).vertNum - 1).y));
		}else{
			drawColorQuad(int(myVerts.at(f.verts.at(0).vertNum - 1).x),
						  int(myVerts.at(f.verts.at(0).vertNum - 1).y),
						  int(myVerts.at(f.verts.at(1).vertNum - 1).x),
						  int(myVerts.at(f.verts.at(1).vertNum - 1).y),
						  int(myVerts.at(f.verts.at(2).vertNum - 1).x),
						  int(myVerts.at(f.verts.at(2).vertNum - 1).y),
						  int(myVerts.at(f.verts.at(3).vertNum - 1).x),
						  int(myVerts.at(f.verts.at(3).vertNum - 1).y));
		}
	}
}
/***************************************************************************/
void drawColorPerspectiveFaces(){
	std::sort(myFaces.begin(), myFaces.end(), faceSort);

	for(auto& f: myFaces){
		r = f.c.r;
		b = f.c.b;
		g = f.c.g;
		if(f.verts.size() == 3){
			drawColorTriangle(int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).x),
							  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).y),
							  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).x),
							  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).y),
							  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).x),
							  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).y));
			
		}else{
			drawColorQuad(int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).x),
						  int(perspectiveVerts.at(f.verts.at(0).vertNum - 1).y),
						  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).x),
						  int(perspectiveVerts.at(f.verts.at(1).vertNum - 1).y),
						  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).x),
						  int(perspectiveVerts.at(f.verts.at(2).vertNum - 1).y),
						  int(perspectiveVerts.at(f.verts.at(3).vertNum - 1).x),
						  int(perspectiveVerts.at(f.verts.at(3).vertNum - 1).y));
		}
	}
}
/***************************************************************************/
void drawColorQuad(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4){
	drawColorTriangle(x1, y1, x2, y2, x3, y3);
	drawColorTriangle(x3, y3, x4, y4, x1, y1);
}
/***************************************************************************/
void drawColorTriangle(int x1, int y1, int x2, int y2, int x3, int y3){
   int maxX = (x1 > x2) ? x1 : x2;
   int maxY = (y1 > y2) ? y1 : y2;
   int minX = (x1 < x2) ? x1 : x2;
   int minY = (y1 < y2) ? y1 : y2;

   maxX = (maxX > x3) ? maxX : x3;
   minX = (minX < x3) ? minX : x3;
   maxY = (maxY > y3) ? maxY : y3; 
   minY = (minY < y3) ? minY : y3;
   
   for(int x = minX; x <= maxX; x++){
	   for(int y = minY; y <= maxY; y++){
		   	double denomenator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
   			double a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denomenator;
   			double b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denomenator;
   			double c = 1 - a - b;
			if((a>=0)&&(a<=1)&&(b>=0)&&(b<=1)&&(c>=0)&&(c<=1)){
				write_pixel(x, y, 1.0);
			}

	   }
   }
   
   
}

/* Drawing line from (x1, y1) to (x2, y2) using DDA */
/***************************************************************************/
void drawALine(int x1, int y1, int x2, int y2){
	
	//calculating equation of the line
	float slope = (float)(y2-y1)/(float)(x2-x1);
	float b = (float)y1 - (slope * x1);

	bool steep;
	
	int xChange, yChange;
	
	//whether line is steep or shallow alters the algorithm
	if(abs(slope) > 1){
		steep = true;
	}else{
		steep = false;
	}

	//determining if will increment or decrement each variable
	if(x2 > x1){
		xChange = 1;
	}else{
		xChange = -1;
	}

	if(y2 > y1){
		yChange = 1;
	}else{
		yChange = -1;
	}
		

	double change;
	int currX = x1;
	int currY = y1;
	double intensity = 1.0;
	double temp = 0;

	//vertical line
	if(x1 == x2){
		while(currY != y2){
			write_pixel(currX, currY, intensity);
			currY += yChange;
		}
	}
	//horizontal line
	else if(y1 == y2){
		while(currX != x2){
			write_pixel(currX, currY, intensity);
			currX += xChange;
		}
	}
	else if(!steep){
		write_pixel(currX,currY,intensity);
		while(currX != x2){
			//determining next pixel
			currX += xChange;
			currY = (int)((slope * currX) + b);
			//taking care of the jaggies
			change = modf((slope * currX) + b, &temp);
			write_pixel(currX, currY, 1 - change);
			write_pixel(currX, currY + 1, change);

		}
	}else{
		write_pixel(currX,currY, intensity);
		while(currY != y2){
			//determining next pixel
			currY += yChange;
			currX = (int)((currY - b)/slope);
			//taking care of jaggies
			change = modf((currY -b)/slope, &temp);
			write_pixel(currX, currY, 1 - change);
			write_pixel(currX + 1, currY, change);

		}
	}
	
	
}

