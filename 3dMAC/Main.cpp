#include "Main.h"
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    myInit();
    glutDisplayFunc(myDisplay);
    glutIdleFunc(myIdle);
    glutMouseFunc(myMouse);
    glutMotionFunc(myMotion);
    glutKeyboardFunc(myKeyboard);
    glLoadIdentity();
    gluPerspective(30.0, (double)glutGet(GLUT_WINDOW_WIDTH) / (double)glutGet(GLUT_WINDOW_HEIGHT), 1.0, 100.0);
    glutMainLoop();
    return 0;
}

void drawContainer() {
	glColor3f(0.0f, 0.0f, 0.0f);
	glLineWidth(1.0f);
    glPushMatrix();
    glTranslatef(0.0, 0.0, -10.0);
    glRotatef(theta, 1.0, 0.0, 0.0 );
    glRotatef(theta, 0.0, 1.0, 0.0 );
    glRotatef(theta, 0.0, 0.0, 1.0 );
    glScaled(2.0/(field.GridNum() * field.Dx()), 2.0/(field.GridNum() * field.Dx()), 2.0/(field.GridNum() * field.Dx()));
    glutWireCube(field.GridNum() * field.Dx());
    glPopMatrix();
}

void drawVelocity() {
	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);
    glPushMatrix();
    glTranslatef(0.0, 0.0, -10.0);
    glScaled(2.0/field.GridNum(), 2.0/field.GridNum(), 2.0/field.GridNum());
	glBegin(GL_LINES);
        for(int i = 0; i < field.GridNum(); i++) {
            for(int j = 0; j < field.GridNum(); j++) {
                for(int k = 0; k < field.GridNum(); k++) {
                    Vector3d fieldPosition((i + 0.5) * field.Dx(), (j + 0.5) * field.Dx(), (k + 0.5) * field.Dx());
                    Vector3d fieldVelocityPosition = fieldPosition + field.GetVelocity(fieldPosition);
			        glVertex3d(fieldPosition.x() - (field.GridNum() * field.Dx())/2.0,
                               fieldPosition.y() - (field.GridNum() * field.Dx())/2.0,
                               fieldPosition.z() - (field.GridNum() * field.Dx())/2.0);
			        glVertex3d(fieldVelocityPosition.x() - (field.GridNum() * field.Dx())/2.0,
                               fieldVelocityPosition.y() - (field.GridNum() * field.Dx())/2.0,
                               fieldVelocityPosition.z() - (field.GridNum() * field.Dx())/2.0);
                }
            }
        }
	glEnd ();
    glPopMatrix();
}

void drawPoints() {
	glColor3f(1.0f, 0.0f, 0.0f);
    glPointSize(2.0f);
    glPushMatrix();
    glTranslatef(0.0, 0.0, -10.0);
    glRotatef(theta, 1.0, 0.0, 0.0 );
    glRotatef(theta, 0.0, 1.0, 0.0 );
    glRotatef(theta, 0.0, 0.0, 1.0 );
    glScaled(2.0/(field.GridNum() * field.Dx()), 2.0/(field.GridNum() * field.Dx()), 2.0/(field.GridNum() * field.Dx()));
	glBegin(GL_POINTS);
         for(int i = 0; i < field.GridNum(); i++) {
             for(int j = 0; j < field.GridNum(); j++) {
                 for(int k = 0; k < field.GridNum(); k++) {
 			        glVertex3d(points.at(i).at(j).at(k).x() - (field.GridNum() * field.Dx())/2.0,
                                points.at(i).at(j).at(k).y() - (field.GridNum() * field.Dx())/2.0,
                                points.at(i).at(j).at(k).z() - (field.GridNum() * field.Dx())/2.0);
                 }
             }
         }
	glEnd ();
    glPopMatrix();
}

void drawMarble() {
}

void myDisplay(void) {
	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    drawContainer();
    if(DRAW_MODE == VELOCITY) drawVelocity();
    if(DRAW_MODE == POINTS) drawPoints();
    if(DRAW_MODE == MARBLE) drawMarble();
    glFlush();
}

void myKeyboard(unsigned char key, int x, int y) {
    if(key == 27) exit(0);  
    if(key == 'v') DRAW_MODE = VELOCITY;
    if(key == 'p') DRAW_MODE = POINTS;
    if(key == 'm') DRAW_MODE = MARBLE;
    if(key == 'r') {
        field.Init();
        initPoints();
        initMarble();
    }
}

void initPoints() {
    points = vector< vector< vector<Vector3d>>>(field.GridNum(), vector< vector<Vector3d>>(field.GridNum(), vector<Vector3d>(field.GridNum())));
    for(int i = 0; i < field.GridNum(); i++) {
        for(int j = 0; j < field.GridNum(); j++) {
            for(int k = 0; k < field.GridNum(); k++) {
                points.at(i).at(j).at(k) = Vector3d((i + 0.5) * field.Dx(), (j + 0.5) * field.Dx(), (k + 0.5) * field.Dx()); 
            }
        } 
    }
}

void initMarble() {
}

void myInit() {
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(512, 512);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MAC");
    field.Init();
    initPoints();
    initMarble();
    DRAW_MODE = POINTS;
}

void updateDeltaTime() {
    const auto currentTime = std::chrono::system_clock::now(); 
    deltaTime = (std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime).count())/1000.0;
    lastTime = currentTime;
}

double elapsedTime = 0.0;
void updateForce(double timeStep) {
    double radius = (field.GridNum() * field.Dx())/2.0;
    Vector3d center{radius, radius, radius};
    forceSourcePosition = center 
                          + sin(elapsedTime * 0.01 * (2 * M_PI)) * Vector3d(radius, 0.0, 0.0) 
                          + cos(elapsedTime * 0.01 * (2 * M_PI)) * Vector3d(0.0, radius, 0.0);
    Vector3d force(1.0 * sin(elapsedTime * 0.01 * (2 * M_PI)), -1.0 * sin(elapsedTime * 0.01 * (2 * M_PI)), 0.0);
    elapsedTime += timeStep;
    field.SetForce(force, forceSourcePosition); 
}

void updateField(double timeStep) {
    field.Advect(timeStep);
    field.AddForce(timeStep);
    field.CG_Project(timeStep);
}

void updatePoints(double timeStep) {
    for(int i = 0; i < field.GridNum(); i++) {
        for(int j = 0; j < field.GridNum(); j++) {
            for(int k = 0; k < field.GridNum(); k++) {
                points.at(i).at(j).at(k) += timeStep * field.GetVelocity(points.at(i).at(j).at(k)); 
             }
        } 
    }
}

void updateMarble(double timeStep) {
}

int imageId = 0;
void myIdle(void) {
    double timeStep = 1.0;
    updateDeltaTime();
    updateForce(timeStep);
    updateField(timeStep);
    updatePoints(timeStep);
    glutPostRedisplay();
    theta += 5.0;
    ostringstream sout;
    sout << setfill('0') << setw(5) << imageId;
    string s = sout.str();
    saveImage(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), "images/" + s);
    imageId++;
    cout << "\rdeltaTime: " << deltaTime;
    fflush(stdout);
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

void saveImage(const unsigned int imageWidth, const unsigned int imageHeight, const string outImageName)
{
    std::string fname = outImageName + ".jpg";
    cv::Mat outImage(imageHeight, imageWidth, CV_8UC3);
    glPixelStorei(GL_PACK_ROW_LENGTH, outImage.step/outImage.elemSize()); 
    glReadPixels(0, 0, imageWidth, imageHeight, GL_BGR, GL_UNSIGNED_BYTE, outImage.data);
    cv::imwrite( fname.c_str(), outImage );
}
