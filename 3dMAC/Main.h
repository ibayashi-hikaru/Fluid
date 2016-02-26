#ifndef ST_MAIN_H_INCLUDED
#define ST_MAIN_H_INCLUDED
#include <iostream>
#include <sstream>
#include <iomanip>
#include <GLUT/glut.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <eigen3/Eigen/Core>
#include <chrono>
#include <vector>
#include <list> 
#include <cmath>
#include "Field.h"
using namespace std;
using namespace Eigen;

enum DrawMode {
    VELOCITY = 0,
    POINTS = 1,
    MARBLE = 2,
};
auto lastTime = std::chrono::system_clock::now();
double deltaTime; 
Vector2d lastPosition = Vector2d::Zero();
Vector2d currentPosition = Vector2d::Zero();
Field field(16);
vector< vector< vector<Vector3d>>> points;
int marbleCount = 50;
list<Vector2d> marbleEdge;
DrawMode DRAW_MODE;
double theta = 60.0;
Vector3d forceSourcePosition{(field.GridNum() * field.Dx())/2.0,
                              field.GridNum() * field.Dx(),
                              (field.GridNum() * field.Dx())/2.0};
bool saveflg = false;

void drawContainer();
void drawVelocity();
void drawPoints();
void drawMarble();
void drawForceSource();
void myDisplay(void);
void myKeyboard(unsigned char key, int x, int y);
void initPoints();
void initMarble();
void myInit();
void updateDeltaTime();
void updateField(double timeStep);
void updateForce(double timeStep);
void updatePoints(double timeStep);
void updateMarble(double timeStep);
void myIdle(void);
void myMouse(int button, int state, int x, int y);
void myMotion(int x, int y);
void saveImage( const unsigned int imageWidth, const unsigned int imageHeight, const string outImageName);
#endif // ST_MAIN_H_INCLUDED
