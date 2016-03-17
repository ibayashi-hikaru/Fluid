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
#include "Interface.h"
extern Interface interface;
void myDisplay(void);
void myKeyboard(unsigned char key, int x, int y);
void myIdle(void);
void myMouse(int button, int state, int x, int y);
void myMotion(int x, int y);
void drawContainer();
void drawVelocity();
void drawPoints();
void drawForceSource();
void initPoints();
void myInit();
void updateDeltaTime();
void updateForce(double timeStep);
void updateField(double timeStep);
void updatePoints(double timeStep);
void saveImage(const int imageWidth, const int imageHeight, const std::string outImageName);
void polarview();
#endif // ST_MAIN_H_INCLUDED
