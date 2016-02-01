#include "Field.h"
using namespace std;
void
Field::AddForce(double dt) {
}

void
Field::Advect(double dt) {
    // horizontal direction velocity
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            double currentX = i * dx; 
            double currentY = (j + 0.5) * dx;
            double currentVelocityX = getVelocityX(currentX, currentY);
            double currentVelocityY = getVelocityY(currentX, currentY);
            double lastX = currentX - dt * currentVelocityX;
            double lastY = currentY - dt * currentVelocityY;
            ux1.at(i).at(j) = getVelocityX(lastX, lastY);
        }
    }
    // vertical direction velocity
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            double currentX = (i + 0.5) * dx; 
            double currentY = j * dx;
            double currentVelocityX = getVelocityX(currentX, currentY);
            double currentVelocityY = getVelocityY(currentX, currentY);
            double lastX = currentX - dt * currentVelocityX;
            double lastY = currentY - dt * currentVelocityY;
            uy1.at(i).at(j) = getVelocityY(lastX, lastY);
        }
    }
}

void
Field::Project(double dt) {
    double scale = dt / (rho * dx * dx);
    double eps = 1.0e-4;
    double err;
    do {
        err = 0.0;

        for (int j = 0; j < Nx; j++) {
            for (int i = 0; i < Ny; i++) {
                vector<double> D = {1.0, 1.0, -1.0, -1.0}; 
                vector<double> F = {static_cast<double>(i < Nx - 1),
                                    static_cast<double>(j < Ny - 1),
                                    static_cast<double>(i > 0),
                                    static_cast<double>(j > 0)};
                vector<double> P = {F.at(0) ? p.at(i + 1).at(j) : 0.0,
                                    F.at(1) ? p.at(i).at(j + 1) : 0.0,
                                    F.at(2) ? p.at(i - 1).at(j) : 0.0,
                                    F.at(3) ? p.at(i).at(j - 1) : 0.0 };
                vector<double> U = {ux1.at(i + 1).at(j), uy1.at(i).at(j + 1), ux1.at(i).at(j), uy1.at(i).at(j)};
                
                double det = 0.0;
                double sum_L = 0.0;
                double sum_R = 0.0;
                for(int n = 0; n < 4; n++) {
                    det += F.at(n) * scale;
                    sum_L += F.at(n) * P.at(n) * scale;
                    sum_R += F.at(n) * D.at(n) * U.at(n)/dx;
                }
                err = fmax(err, fabs(det*p.at(i).at(j) - sum_L + sum_R)); 
                p.at(i).at(j) = (sum_L - sum_R)/det;
            } 
        }
    } while(eps < err);

    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
             ux1.at(i).at(j) = ux1.at(i).at(j) - (dt/rho) * (p.at(i).at(j) - p.at(i-1).at(j));
        } 
    }    
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
             uy1.at(i).at(j) = uy1.at(i).at(j) - (dt/rho) * (p.at(i).at(j) - p.at(i).at(j-1));
        } 
    }    
}

void
Field::SetForce(Vector2d force, Vector2d position) {
}

Vector2d
Field::TransformDisplayToField(Vector2d displayPosition, int width, int height) const {
    displayPosition.x() /= (Nx * dx)/width;
    displayPosition.y() /= (Ny * dx)/height;
    return displayPosition;
}
// grid外は境界と同じ値
double
Field::getVelocityX(double x, double y) const {
    y -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    unsigned long i = x;
    unsigned long j = y;
    vector<double> f = {ux1.at(i).at(j), ux1.at(i).at(j + 1), ux1.at(i + 1).at(j), ux1.at(i + 1).at(j + 1)};
    x = x - i;
    y = y - j;
    vector<double> c = {(1.0 - x) * (1.0 - y), (1.0 - x) * y, x * (1.0 - y), x * y};
    return c.at(0) * f.at(0) + c.at(1) * f.at(1) + c.at(2) * f.at(2) + c.at(3) * f.at(3); 
}

double
Field::getVelocityY(double x, double y) const {
    x -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1e-6, y/dx));
    unsigned long i = x;
    unsigned long j = y;
    vector<double> f = {uy1.at(i).at(j), uy1.at(i).at(j + 1), uy1.at(i + 1).at(j), uy1.at(i + 1).at(j + 1)};
    x = x - i;
    y = y - j;
    vector<double> c = {(1.0 - x) * (1.0 - y), (1.0 - x) * y, x * (1.0 - y), x * y};
    return c.at(0) * f.at(0) + c.at(1) * f.at(1) + c.at(2) * f.at(2) + c.at(3) * f.at(3); 
}

void
Field::makeBoundary() {
    for(int i = 0; i < Ny; i++) {
        ux0.at(0).at(i) = 0.0; 
        ux0.at(Nx).at(i) = 0.0; 
        ux1.at(0).at(i) = 0.0; 
        ux1.at(Nx).at(i) = 0.0; 
    }
    for(int i = 0; i < Nx; i++) {
        uy0.at(i).at(0) = 0.0; 
        uy0.at(i).at(Ny) = 0.0; 
        uy1.at(i).at(0) = 0.0; 
        uy1.at(i).at(Ny) = 0.0; 
    }
}