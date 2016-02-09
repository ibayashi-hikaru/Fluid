#include "Main.h"
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

void drawVelocity() {
	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-1.0, -1.0, 0);
    glScaled(2.0/field.GridNum(), 2.0/field.GridNum(), 1.0);
	glBegin(GL_LINES);
        double windowDx = 2.0 / field.GridNum();
        for(int i = 0; i < field.GridNum(); i++) {
            for(int j = 0; j < field.GridNum(); j++) {
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
    glScaled(2.0/field.GridNum(), 2.0/field.GridNum(), 1.0);
	glBegin(GL_POINTS);
        for(int i = 0; i < field.GridNum(); i++) {
            for(int j = 0; j < field.GridNum(); j++) {
			    glVertex2d(points.at(i).at(j).x(), points.at(i).at(j).y());
            }
        }
	glEnd ();
    glPopMatrix();
    glFlush();
}

void drawMarble() {
    auto forward = marbleEdge.begin();
    auto backward = forward;
	glColor3f(0.0f, 0.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-1.0, -1.0, 0);
    glScaled(2.0/field.GridNum(), 2.0/field.GridNum(), 1.0);
    glEnable(GL_COLOR_LOGIC_OP);
    glLogicOp(GL_INVERT);
    // 三角形ポリゴンを描く
    glBegin(GL_TRIANGLES);
    for(forward++; forward!=marbleEdge.end(); forward++) {
        glVertex2d(0.0,0.0);
        glVertex2d(backward->x(), backward->y());
        glVertex2d(forward->x(), forward->y());
        backward = forward;
    }
        glVertex2d(0.0,0.0);
        glVertex2d(backward->x(), backward->y());
        glVertex2d(marbleEdge.begin()->x(), marbleEdge.begin()->y());
    glEnd();
    glDisable(GL_COLOR_LOGIC_OP);
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
    if(key == 'r') {
        field.Init();
        initPoints();
        initMarble();
    }
}

void initPoints() {
    points = vector< vector<Vector2d>>(field.GridNum(), vector<Vector2d>(field.GridNum()));
    for(int i = 0; i < field.GridNum(); i++) {
        for(int j = 0; j < field.GridNum(); j++) {
            points.at(i).at(j) = Vector2d((i + 0.5) * 1.0, (j + 0.5) * 1.0); 
        } 
    }
}

void initMarble() {
    marbleEdge.clear();
    for(int i = 0; i < marbleCount; i++) {
        marbleEdge.push_back(
                Vector2d(
                    field.GridNum() / 2.0 + (field.GridNum() * 0.25) * sin((2 * M_PI) * (i/(double)(marbleCount))),
                    field.GridNum() / 2.0 + (field.GridNum() * 0.25) * cos((2 * M_PI) * (i/(double)(marbleCount)))
                    )
                ); 
    }
}

void myInit() {
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(512, 512);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MAC");
    field.Init();
    initPoints();
    initMarble();
    DRAW_MODE = VELOCITY;
}

void updateDeltaTime() {
    const auto currentTime = std::chrono::system_clock::now(); 
    deltaTime = (std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime).count())/1000.0;
    lastTime = currentTime;
}

void updateForce() {
    Vector2d lastFieldPosition = field.TransformDisplayToField(lastPosition, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    Vector2d currentFieldPosition = field.TransformDisplayToField(currentPosition, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    Vector2d force = 10.0 * (currentFieldPosition - lastFieldPosition)/deltaTime; //速度に比例した力 10.0は比例定数
    Vector2d position = currentFieldPosition;
    field.SetForce(force, position);
    lastPosition.x() = currentPosition.x();
    lastPosition.y() = currentPosition.y();
}

void updateField(double timeStep) {
    field.Advect(timeStep);
    field.AddForce(timeStep);
    field.GS_Project(timeStep); 
    // field.CG_Project(timeStep);
}

void updatePoints(double timeStep) {
    for(int i = 0; i < field.GridNum(); i++) {
        for(int j = 0; j < field.GridNum(); j++) {
            points.at(i).at(j) += timeStep * field.GetVelocity(points.at(i).at(j)); 
        } 
    }
}
void resample() {
    double min_ds = 0.1;
    double max_ds = 1.0;
    auto forward = marbleEdge.begin();
    auto backward = forward;
    for( forward++; forward!=marbleEdge.end(); ) {
        Vector2d &p0 = *backward;
        Vector2d &p1 = *forward;
        double d = (p0 - p1).norm();
        if( d < min_ds ) { 
        // 頂点間の距離が狭すぎる場合
        // 2つの頂点を結合する
            p0 = 0.5 * (p0 + p1);
            forward = marbleEdge.erase(forward);
        } else if( d > max_ds ) { // 頂点間の幅が広すぎる場合 // 間に新しい頂点を挿入する
            Vector2d p = 0.5 * (p0 + p1);
            forward = marbleEdge.insert(forward,p);
        } else {
            backward = forward; 
            forward ++;
        } 
    }
}

void updateMarble(double timeStep) {
    for(auto itr = marbleEdge.begin(); itr != marbleEdge.end(); ++itr) {
        *itr += timeStep * field.GetVelocity(*itr);
    }
    resample();
}

void myIdle(void) {
    double timeStep = 1.0;
    updateDeltaTime();
    if(currentPosition != Vector2d::Zero() && currentPosition != lastPosition) {
        updateForce();
    }
    updateField(timeStep);
    updatePoints(timeStep);
    updateMarble(timeStep);
    glutPostRedisplay();
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

