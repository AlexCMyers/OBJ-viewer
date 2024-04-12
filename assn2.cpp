//
//		          Programming Assignment #2 
//
//			      Alex Myers
//				  11-8-19
//	This project displays a .obj file passed on the command line and allows
// transformations of scaling, rotating, and translating.
// perspective view defaults to on and can be toggled with the 'v' key
// wireframe/solid view is toggled with the 'z' key
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
#include <time.h>
#include <algorithm>
#define WIDTH 600
#define HEIGHT 600

/* structs */
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
	color_t c;
}face;

/* Globals to keep track of recent points and the setting */
int x_last,y_last;
char currSetting;

//toggles for viewing setting
bool perspective;
bool color;
vert center;
double r, g, b;


std::vector<vert> myVerts;
std::vector<vert> perspectiveVerts;
std::vector<normal> myNormals;
std::vector<texture> myTextures;
std::vector<face> myFaces;


/* function prototypes */

void parse(std::ifstream& infile);
void normalizeVerts();
void drawFaces();
void drawPerspectiveFaces();
void drawColorFaces();
void drawColorPerspectiveFaces();
void drawColorQuad(int x1, int y1, int x2, int y2, 
						int x3, int y3, int x4, int y4);
void drawColorTriangle(int x1, int y1, int x2, int y2, int x3, int y3);

void updatePerspectiveVerts();
void drawALine(int x1, int y1, int x2, int y2);

void transform(char setting, char direction);
void translate(char direction);
void rotate(char direction);
void scale (char direction);

void vertMult(double matrix[4][4]);

void matrixMult(double mat1[4][4], double mat2[4][4], double result[4][4]);

void write_pixel(int x, int y, double intensity);

bool faceSort(face& left, face& right);





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

	
	//draw functions chosen based on current settings
	if(!perspective){
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



/* Given by Dr. Zordan; altered to record points
	While unused for this project, keeping format for
	future use with clicking within window*/
/***************************************************************************/
void mouse(int button, int state, int x, int y)
{
	
}
 
/* Alters the current tranformation based on key input */
/***************************************************************************/
void keyboard ( unsigned char key, int x, int y )  
{

	switch ( key ) {
		case 'v':
			perspective = !perspective;
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
			color = !color;
			if(!color){
				r = b = g = 1.0;
			}
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
/* Reads in and stores data from .obj file */
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
			char vertCombo[(18 * 3) + 2];
			//while still verts for the face (should be 3 or 4)
			//verts on each face line can be in the form #/#/#,
			// #//#, #/#, or # and this code handles any of those
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
			//setting random color for each face
			newFace.c.r = double((rand() % 1000))/(double)1000;
			newFace.c.g = double((rand() % 1000))/(double)1000;
			newFace.c.b = double((rand() % 1000))/(double)1000;

			myFaces.push_back(newFace);
		}
		
	}
}
/* after parsing the .obj file, need to center to object and scale it
 * to a reasonable size for our window */
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
	//window is 500 x 500. Want to set it up so everything fits inside a 
	//300 x 300 x 300 cube.
	//which means scaling points so that the absolute max value is 150 
	//(to fit in window nicely)
	//recall: scaling is just multiplying points. If all three 
	//coordinates are scaled
	//for every point, then the entire object will scale equally 
	//in all directions

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
	}
	//initializing globals
	center.x = 0.0;
	center.y = 0.0;
	center.z = 0.0;

	updatePerspectiveVerts();

}
/* called after each transform; these verts use the given vert data
 * to store locations with faux perspective */
/***************************************************************************/
void updatePerspectiveVerts(){
	//after initial normalization, the z values range from -150 to 150
	//let's assume a vanishing point of -800 in the z direction.
	//everything along the plane z = 0 will keep it's x and y values.
	//-800 - z/-800 should then be the scale factor for x and y, 
	// because where z = 0, that scale factor is 1.  Where z = 100, 
	//that scale factor is > 1, where z = -800,
	//that factor is 0, and so on. 
	for(int i = 0; i < (int)myVerts.size(); i++){
		double scaleFactor = (-800.0 - myVerts.at(i).z)/-800.0;
		perspectiveVerts.at(i).x = myVerts.at(i).x * scaleFactor;
		perspectiveVerts.at(i).y = myVerts.at(i).y * scaleFactor;
	}
}
/* used for drawing in color. Furthest away faces are checked first */
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
/* handles all transformations, called from keyboard function */
/***************************************************************************/
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

/* rotates object about its center */
/***************************************************************************/
void rotate(char direction){
	double rotationMatrix[4][4];
	double originMatrix[4][4];
	double resetMatrix[4][4];
	double resultMatrix[4][4];
	double finalMatrix[4][4];
	double degree;
	char axis;
	//10 degrees chosen as it's small enough to be precise, but large enough
	//to rotate at a reasonable speed
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

	//building rotation matrix
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

	//Setting up matrices to multiply by in order to rotate object about
	//its center rather than the origin
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
/* moving the object horizontally or vertically */
/***************************************************************************/
void translate(char direction){
	double translationMatrix[4][4];
	int dx = 0, dy = 0, dz = 0;
	//10 pixels chosen as it's small enough to be precise while
	//large enough to translate efficiently
	switch(direction){
		case 'w':
			dy = 10.0;
			//must also update the globals for the object center
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
	//first 3 columns and rows are standard
	for(int i = 0; i < 4; i++){
		for (int j = 0; j < 3; j++){
			translationMatrix[i][j] = (i == j)? 1.0 : 0.0;
		}
	}
	translationMatrix[0][3] = dx;
	translationMatrix[1][3] = dy;
	//unnecessary now, but could one day adapt to move on z axis
	translationMatrix[2][3] = dz;
	translationMatrix[3][3] = 1;

	vertMult(translationMatrix);


}
/* used to scale object and make it larger or smaller */
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
	//building scaling matrix
	for(int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			scaleMatrix[i][j] = (i == j)? scaleFactor : 0.0;
		}
	}
	scaleMatrix[3][3] = 1.0;
	//building movement matrices to scale about the objects center rather
	//than the origin
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
/* performs matrix multiplication for the sake of combining transforms */
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
/* alters all verts based on a transformation matrix */
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
		std::cout << "Error, expected use is <executable> <filename.obj>" << 
		std::endl;
		return 1;
	}
	std::ifstream input(argv[1]);
	if(!input.is_open()){
		std::cout << "Error opening " << argv[1] << std::endl 
		<< "Terminating program" << std::endl;
		return 2;
	}
	//initializing globals
	perspective = true;
	color = true;
	currSetting = 't';
	r = 1.0;
	g = 1.0;
	b = 1.0;

	parse(input);

	normalizeVerts();

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

/* passes the verts for each face to the DDA line drawing function to 
 * produce a non-perspective wireframe */
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
/* Passes the stored perspective verts to the line drawing function to 
 * create a perspective wireframe */
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
/* draws the faces in color using an adapted zbuffer algorithm based on the 
 * fact that no faces can intersect with each other */
/***************************************************************************/
void drawColorFaces(){
	std::sort(myFaces.begin(), myFaces.end(), faceSort);

	for(auto& f: myFaces){
		//setting color to current face
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
/* drawing the colored faces but using the perspective verts */
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
/* Draws a quad in color by drawing two appropriate triangles in color */
/***************************************************************************/
void drawColorQuad(int x1, int y1, int x2, int y2, int x3, int y3, 
												int x4, int y4){
	drawColorTriangle(x1, y1, x2, y2, x3, y3);
	drawColorTriangle(x3, y3, x4, y4, x1, y1);
}
/* fills in a triangle based on an old algorithm based on Barycentric 
 * coordinates */
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
   //only need to check the rectangle within which the triangle is
   for(int x = minX; x <= maxX; x++){
	   for(int y = minY; y <= maxY; y++){
		   	double denomenator = ((y2 - y3) * (x1 - x3) + 
												(x3 - x2) * (y1 - y3));
   			double a = ((y2 - y3) * (x - x3) + (x3 - x2) * 
												(y - y3)) / denomenator;
   			double b = ((y3 - y1) * (x - x3) + (x1 - x3) * 
												(y - y3)) / denomenator;
   			double c = 1 - a - b;
			//if the point is inside the triangle
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

