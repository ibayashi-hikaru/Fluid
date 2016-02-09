#ifndef ST_MAIN_H_INCLUDED
#define ST_MAIN_H_INCLUDED
#include <iostream>
#include <GLUT/glut.h>
#include <Eigen/Core>
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
Field field(32);
vector< vector<Vector2d>> points;
int marbleCount = 50;
list<Vector2d> marbleEdge;
DrawMode DRAW_MODE;

void drawVelocity();
void drawPoints();
void drawMarble();
void myDisplay(void);
void myKeyboard(unsigned char key, int x, int y);
void initPoints();
void initMarble();
void myInit();
void updateDeltaTime();
void updateForce();
void updateField(double timeStep);
void updatePoints(double timeStep);
void resample();
void updateMarble(double timeStep);
void myIdle(void);
void myMouse(int button, int state, int x, int y);
void myMotion(int x, int y);
#endif // ST_MAIN_H_INCLUDED
