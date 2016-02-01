#include <iostream>
#include <GLUT/glut.h>
#include <Eigen/Core>
#include <chrono>
#include "Field.h"
using namespace std;
using namespace Eigen;

auto lastTime = std::chrono::system_clock::now();
double deltaTime; 
Vector2d lastPosition;
Vector2d currentPosition;
Field field(64, 64);

void myDisplay(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_POLYGON);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.5, 0.0, 0.0);
        glVertex3f(0.5, 0.5, 0.0);
        glVertex3f(0.0, 0.5, 0.0);
    glEnd();
    glFlush();
}

void myKeyboard(unsigned char key, int x, int y) {
    if(key == 27) {
        exit(0);  
    }
}

void myInit() {
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(512, 512);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MAC");
}

void myIdle(void) {
    const auto currentTime = std::chrono::system_clock::now(); 
    deltaTime = (std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime).count())/1000.0;
    lastTime = currentTime;
    if(currentPosition != Vector2d::Zero()) {
        Vector2d force = (currentPosition - lastPosition)/deltaTime; //速度に比例した力
        Vector2d position = field.TransformDisplayToField((currentPosition + lastPosition)/2.0, 512, 512);
        field.SetForce(force, position);    
    }
    glutPostRedisplay();
}

void myMouse(int button, int state, int x, int y) {
    if(state == GLUT_DOWN) {
        switch(button) {
        case GLUT_LEFT_BUTTON :
            std::cout << "left button" << std::endl;
            std::cout << "x :" << x << ", y :" << y << std::endl;
            lastPosition.x() = currentPosition.x();
            lastPosition.y() = currentPosition.y();
            currentPosition.x() = x;
            currentPosition.y() = y;
            break;
        case GLUT_RIGHT_BUTTON :
            std::cout << "right button" << std::endl;
            std::cout << "x :" << x << ", y :" << y << std::endl;
            break;
        } 
    }
}

void myMotion(int x, int y) {
    lastPosition.x() = currentPosition.x();
    lastPosition.y() = currentPosition.y();
    currentPosition.x() = x;
    currentPosition.y() = y;
    std::cout << "x: " << x << "y: " << std::endl;
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