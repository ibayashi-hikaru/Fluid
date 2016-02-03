#include <iostream>
#include <GLUT/glut.h>
#include <Eigen/Core>
#include <chrono>
#include <vector>
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
const double force_k = 10.0;
int windowSize = 512;
int gridNum = 16;
Field field(gridNum, gridNum);
vector< vector<Vector2d>> points;
DrawMode DRAW_MODE;

void drawVelocity() {
	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-1.0, -1.0, 0);
    glScaled(2.0/gridNum, 2.0/gridNum, 1.0);
	glBegin(GL_LINES);
        double windowDx = 2.0 / gridNum;
        for(int i = 0; i < gridNum; i++) {
            for(int j = 0; j < gridNum; j++) {
                Vector2d fieldPosition((i + 0.5) * 1.0, (j + 0.5) * 1.0);
                Vector2d fieldVelocityPosition = fieldPosition + field.GetVelocity(fieldPosition);
			    glVertex2d(fieldPosition.x(), fieldPosition.y());
			    glVertex2d(fieldVelocityPosition.x(), fieldVelocityPosition.y());
            }
        }
	glEnd ();
    glPopMatrix();
    glFlush();
}

void drawPoints() {
	glColor3f(1.0f, 0.0f, 0.0f);
    glPointSize(2.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-1.0, -1.0, 0);
    glScaled(2.0/gridNum, 2.0/gridNum, 1.0);
	glBegin(GL_POINTS);
        for(int i = 0; i < gridNum; i++) {
            for(int j = 0; j < gridNum; j++) {
			    glVertex2d(points.at(i).at(j).x(), points.at(i).at(j).y());
            }
        }
	glEnd ();
    glPopMatrix();
    glFlush();
}

void myDisplay(void) {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
    if(DRAW_MODE == VELOCITY) drawVelocity();
    if(DRAW_MODE == POINTS) drawPoints();
    if(DRAW_MODE == MARBLE) drawMarble();
}

void myKeyboard(unsigned char key, int x, int y) {
    if(key == 27) exit(0);  
    if(key == 'v') DRAW_MODE = VELOCITY;
    if(key == 'p') DRAW_MODE = POINTS;
    if(key == 'm') DRAW_MODE = MARBLE;
}

void initPoints() {
    points = vector< vector<Vector2d>>(gridNum, vector<Vector2d>(gridNum));
    for(int i = 0; i < gridNum; i++) {
        for(int j = 0; j < gridNum; j++) {
            points.at(i).at(j) = Vector2d((i + 0.5) * 1.0, (j + 0.5) * 1.0); 
        } 
    }
}

void myInit() {
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(windowSize, windowSize);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MAC");
    field.Init();
    initPoints();
    DRAW_MODE = VELOCITY;
}

void updateDeltaTime() {
    const auto currentTime = std::chrono::system_clock::now(); 
    deltaTime = (std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime).count())/1000.0;
    lastTime = currentTime;
}

void updateForce() {
    Vector2d lastFieldPosition = field.TransformDisplayToField(lastPosition, windowSize, windowSize);
    Vector2d currentFieldPosition = field.TransformDisplayToField(currentPosition, windowSize, windowSize);
    Vector2d force = force_k * (currentFieldPosition - lastFieldPosition)/deltaTime; //速度に比例した力
    Vector2d position = currentFieldPosition;
    field.SetForce(force, position);
    lastPosition.x() = currentPosition.x();
    lastPosition.y() = currentPosition.y();
}

void updateField() {
    field.Advect(deltaTime);
    field.AddForce(deltaTime);
    field.Project(deltaTime);
}

void updatePoints() {
    for(int i = 0; i < gridNum; i++) {
        for(int j = 0; j < gridNum; j++) {
            points.at(i).at(j) += deltaTime * field.GetVelocity(points.at(i).at(j)); 
        } 
    }
}

void myIdle(void) {
    updateDeltaTime();
    if(currentPosition != Vector2d::Zero() && currentPosition != lastPosition) {
        updateForce();
    }
    updateField();
    updatePoints();
    glutPostRedisplay();
}

void myMouse(int button, int state, int x, int y) {
    if(state == GLUT_DOWN) {
        switch(button) {
        case GLUT_LEFT_BUTTON :
            lastPosition.x() = x;
            lastPosition.y() = y;
            break;
        case GLUT_RIGHT_BUTTON :
            break;
        } 
    } else {
        currentPosition = Vector2d::Zero();
    }
}

void myMotion(int x, int y) {
    currentPosition.x() = x;
    currentPosition.y() = y;
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    myInit();
    glutDisplayFunc(myDisplay);
    glutIdleFunc(myIdle);
    glutMouseFunc(myMouse);
    glutMotionFunc(myMotion);
    glutKeyboardFunc(myKeyboard);
    glutMainLoop();
    return 0;
}
