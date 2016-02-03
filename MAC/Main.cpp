#include <iostream>
#include <GLUT/glut.h>
#include <Eigen/Core>
#include <chrono>
#include "Field.h"
using namespace std;
using namespace Eigen;

auto lastTime = std::chrono::system_clock::now();
double deltaTime; 
Vector2d lastPosition = Vector2d::Zero();
Vector2d currentPosition = Vector2d::Zero();
const double force_k = 1.0;
int windowSize = 512;
int gridNum = 16;
Field field(gridNum, gridNum);

void myDisplay(void) {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);
    glPushMatrix();
    glTranslatef(-1, -1, 0);
	glBegin(GL_LINES);
        double windowDx = 2.0 / gridNum;
        for(int i = 0; i < gridNum; i++) {
            for(int j = 0; j < gridNum; j++) {
                Vector2d displayPosition((i + 0.5) * windowDx, (j + 0.5) * windowDx);
                Vector2d fieldPosition = field.TransformDisplayToField(displayPosition, 2.0, 2.0);
                Vector2d fieldVelocityPosition = fieldPosition + field.GetVelocity(fieldPosition);
                Vector2d displayVelocityPosition = field.TransformFieldToDisplay(fieldVelocityPosition, 2.0, 2.0);

			    glVertex2d(displayPosition.x(), displayPosition.y());
			    glVertex2d(displayVelocityPosition.x(), displayVelocityPosition.y());
            }
        }
	glEnd ();
    glPopMatrix();
    glFlush();
}

void myKeyboard(unsigned char key, int x, int y) {
    if(key == 27) {
        exit(0);  
    }
}

void myInit() {
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(windowSize, windowSize);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MAC");
}

void myIdle(void) {
    const auto currentTime = std::chrono::system_clock::now(); 
    deltaTime = (std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime).count())/1000.0;
    lastTime = currentTime;
    if(currentPosition != Vector2d::Zero() && currentPosition != lastPosition) {
        Vector2d lastFieldPosition = field.TransformDisplayToField(lastPosition, windowSize, windowSize);
        Vector2d currentFieldPosition = field.TransformDisplayToField(currentPosition, windowSize, windowSize);
        Vector2d force = force_k * (currentFieldPosition - lastFieldPosition)/deltaTime; //速度に比例した力
        Vector2d position = currentFieldPosition;
        field.SetForce(force, position);
        lastPosition.x() = currentPosition.x();
        lastPosition.y() = currentPosition.y();
    }
    field.Advect(deltaTime);
    field.AddForce(deltaTime);
    field.Project(deltaTime);
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
    field.Init();
    glutDisplayFunc(myDisplay);
    glutIdleFunc(myIdle);
    glutMouseFunc(myMouse);
    glutMotionFunc(myMotion);
    glutKeyboardFunc(myKeyboard);
    glutMainLoop();
    return 0;
}
