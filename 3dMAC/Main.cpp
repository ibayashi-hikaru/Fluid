#include "Main.h"
Interface interface;
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    for(int i = 0; i < argc; i++) {
        if(strcmp(argv[i], "-offline") == 0) interface.saveFlg = true;
        if(strcmp(argv[i], "-noforce") == 0) interface.noForceFlg = true;
    }
    myInit();
    glutDisplayFunc(myDisplay);
    glutIdleFunc(myIdle);
    glutMouseFunc(myMouse);
    glutMotionFunc(myMotion);
    glutKeyboardFunc(myKeyboard);
    glLoadIdentity();
    gluPerspective(30.0, static_cast<double>(glutGet(GLUT_WINDOW_WIDTH)) / static_cast<double>(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 100.0);
    glutMainLoop();
    return 0;
}

void drawContainer() {
	glColor3f(0.0f, 0.0f, 0.0f);
	glLineWidth(1.0f);
    glutWireCube(interface.field.GridNum() * interface.field.Dx());
}

void drawVelocity() {
	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);
	glBegin(GL_LINES);
        for(size_t i = 0; i < interface.field.GridNum(); i++) {
            for(size_t j = 0; j < interface.field.GridNum(); j++) {
                for(size_t k = 0; k < interface.field.GridNum(); k++) {
                    Eigen::Vector3d fieldPosition((i + 0.5) * interface.field.Dx(), (j + 0.5) * interface.field.Dx(), (k + 0.5) * interface.field.Dx());
                    Eigen::Vector3d fieldVelocityPosition = fieldPosition + interface.field.GetVelocity(fieldPosition);
			        glVertex3d(fieldPosition.x() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                               fieldPosition.y() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                               fieldPosition.z() - (interface.field.GridNum() * interface.field.Dx())/2.0);
			        glVertex3d(fieldVelocityPosition.x() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                               fieldVelocityPosition.y() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                               fieldVelocityPosition.z() - (interface.field.GridNum() * interface.field.Dx())/2.0);
                }
            }
        }
	glEnd ();
}

void drawPoints() {
	glColor3f(1.0f, 0.0f, 0.0f);
    glPointSize(2.0f);
	glBegin(GL_POINTS);
        for(auto it = interface.field.sortedMarkersX.begin(); it < interface.field.sortedMarkersX.end(); it++) {
 	        glVertex3d(it->x() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                       it->y() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                       it->z() - (interface.field.GridNum() * interface.field.Dx())/2.0);

        }
	glEnd ();
}

void drawForceSource(){
    struct MaterialStruct {
    	GLfloat ambient[4];
    	GLfloat diffuse[4];
    	GLfloat specular[4];
    	GLfloat shininess;
    };
    //ruby(ルビー)
    MaterialStruct ms_ruby  = {
    	{0.1745f,   0.01175f,  0.01175f,   1.0f},
    	{0.61424f,  0.04136f,  0.04136f,   1.0f},
    	{0.727811f, 0.626959f, 0.626959f,  1.0f},
    	76.8f};
    glEnable(GL_LIGHTING);
    glPushMatrix();
    glTranslated(interface.forceSourcePosition.x() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                 interface.forceSourcePosition.y() - (interface.field.GridNum() * interface.field.Dx())/2.0,
                 interface.forceSourcePosition.z() - (interface.field.GridNum() * interface.field.Dx())/2.0);
    glMaterialfv(GL_FRONT, GL_AMBIENT, ms_ruby.ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, ms_ruby.diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, ms_ruby.specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, &ms_ruby.shininess);
    if(!interface.noForceFlg) glutSolidSphere(interface.field.Dx(), 10, 5);
    glPopMatrix();
    glDisable(GL_LIGHTING);
}

void myDisplay(void) {
	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    glTranslatef(0.0, 0.0, -10.0);
    glScaled(2.0/interface.field.GridNum(), 2.0/interface.field.GridNum(), 2.0/interface.field.GridNum());
    polarview();
    drawContainer();
    drawForceSource();
    if(interface.DRAW_MODE == VELOCITY) drawVelocity();
    if(interface.DRAW_MODE == POINTS) drawPoints();
    glPopMatrix();
    glFlush();
}

void myKeyboard(unsigned char key, int x, int y) {
    if(key == 27) exit(0);  
    if(key == 'v') interface.DRAW_MODE = VELOCITY;
    if(key == 'p') interface.DRAW_MODE = POINTS;
    if(key == 'm') interface.DRAW_MODE = MARBLE;
    if(key == 'r') {
        interface.field.Init();
        initPoints();
        initMarble();
    }
    std::cout << x << "/" << y << std::endl;
}

void initPoints() {
    interface.points = std::vector< std::vector< std::vector<Eigen::Vector3d>>>(static_cast<size_t>(interface.field.GridNum()), std::vector< std::vector<Eigen::Vector3d>>(static_cast<size_t>(interface.field.GridNum()), std::vector<Eigen::Vector3d>(static_cast<size_t>(interface.field.GridNum()))));
    for(size_t i = 0; i < static_cast<size_t>(interface.field.GridNum()); i++) {
        for(size_t j = 0; j < static_cast<size_t>(interface.field.GridNum()); j++) {
            for(size_t k = 0; k < static_cast<size_t>(interface.field.GridNum()); k++) {
                interface.points.at(i).at(j).at(k) = Eigen::Vector3d((i + 0.5) * interface.field.Dx(), (j + 0.5) * interface.field.Dx(), (k + 0.5) * interface.field.Dx()); 
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
    interface.field.Init();
    initPoints();
    initMarble();
    interface.DRAW_MODE = POINTS;

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
     
    glEnable(GL_LIGHT0);
}

void updateDeltaTime() {
    const auto currentTime = std::chrono::system_clock::now(); 
    interface.deltaTime = (std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - interface.lastTime).count())/1000.0;
    interface.lastTime = currentTime;
}

void updateForce(double timeStep) {
    double radius = (interface.field.GridNum() * interface.field.Dx())/2.0;
    Eigen::Vector3d center{radius, radius, radius};
    interface.forceSourcePosition = center 
                          + sin(interface.elapsedTime * 0.01 * (2 * M_PI)) * Eigen::Vector3d(radius, 0.0, 0.0) 
                          + cos(interface.elapsedTime * 0.01 * (2 * M_PI)) * Eigen::Vector3d(0.0, radius, 0.0);
    Eigen::Vector3d force(1.0 * sin(interface.elapsedTime * 0.01 * (2 * M_PI)), -1.0 * sin(interface.elapsedTime * 0.01 * (2 * M_PI)), 0.0);
    interface.elapsedTime += timeStep;
    if(!interface.noForceFlg) interface.field.SetForce(force, interface.forceSourcePosition); 
}

void updateField(double timeStep) {
    interface.field.Advect(timeStep);
    interface.field.AddForce(timeStep);
    interface.field.CoutDiv();
    interface.field.CG_ProjectWithMarker(timeStep);
    //interface.field.CG_Project(timeStep);
    interface.field.CoutDiv();
    interface.field.UpdateMarkers(timeStep);
}

void updatePoints(double timeStep) {
    for(size_t i = 0; i < static_cast<size_t>(interface.field.GridNum()); i++) {
        for(size_t j = 0; j < static_cast<size_t>(interface.field.GridNum()); j++) {
            for(size_t k = 0; k < static_cast<size_t>(interface.field.GridNum()); k++) {
                interface.points.at(i).at(j).at(k) += timeStep * interface.field.GetVelocity(interface.points.at(i).at(j).at(k)); 
            }
        } 
    }
}

void myIdle(void) {
    const double timeStep = 1.0;
    updateDeltaTime();
    updateForce(timeStep);
    updateField(timeStep);
    updatePoints(timeStep);
    glutPostRedisplay();
    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(5) << interface.imageId;
    std::string s = sout.str();
    if(interface.saveFlg) saveImage(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), "images/" + s);
    interface.imageId++;
    std::cout << "\rdeltaTime: " << interface.deltaTime;
    fflush(stdout);
}

void myMouse(int button, int state, int x, int y) {
    if(state == GLUT_DOWN) {
        switch(button) {
        case GLUT_LEFT_BUTTON :
            interface.lastPosition.x() = x;
            interface.lastPosition.y() = y;
            break;
        case GLUT_RIGHT_BUTTON :
            break;
        } 
    } else {
        interface.currentPosition = Eigen::Vector2d::Zero();
    }
}

void myMotion(int x, int y) {
    interface.currentPosition.x() = x;
    interface.currentPosition.y() = y;
	int xDiff, yDiff;
	
	xDiff = interface.currentPosition.x() - interface.lastPosition.x();
	yDiff = interface.currentPosition.y() - interface.lastPosition.y();

	switch (interface.mButton) {
	case GLUT_LEFT_BUTTON:
		interface.azimuth += xDiff/2.0;
		interface.elevation -= yDiff/2.0;
		break;
	case GLUT_RIGHT_BUTTON:
		interface.distance += yDiff/40.0;
		break;
	}
    interface.lastPosition.x() = interface.currentPosition.x();
    interface.lastPosition.y() = interface.currentPosition.y();
}

void saveImage(const int imageWidth, const int imageHeight, const std::string outImageName)
{
    std::string fname = outImageName + ".jpg";
    cv::Mat outImage(imageHeight, imageWidth, CV_8UC3);
    glPixelStorei(GL_PACK_ROW_LENGTH, static_cast<int>(outImage.step/outImage.elemSize())); 
    glReadPixels(0, 0, imageWidth, imageHeight, GL_BGR, GL_UNSIGNED_BYTE, outImage.data);
    cv::flip(outImage, outImage, 0); 
    cv::imwrite( fname.c_str(), outImage );
}

void polarview() {
    glTranslated(0.0, 0.0, -interface.distance);
    glRotated(-interface.twist, 0.0, 0.0, 1.0);
    glRotated(-interface.elevation, 1.0, 0.0, 0.0);
    glRotated(-interface.azimuth, 0.0, 1.0, 0.0);
}
