#include "Field.h"
using namespace std;

void
Field::Init() {
    makeBoundary();
    clearForce();
    initVelocity();
    initPressure();
    initMarkers();
}

void
Field::AddForce(double dt) {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                ux[i][j][k] += dt * forcex[i][j][k];
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                uy[i][j][k] += dt * forcey[i][j][k];
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                uz[i][j][k] += dt * forcez[i][j][k];
            }
        }
    }
    addGravityForce(dt);
    clearForce(); 
}

void
Field::Advect(double dt) {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                Vector3d currentPosition(i * dx, (j + 0.5) * dx, (k + 0.5) * dx);
                ux[i][j][k] = getVelocityX(getLastPosition(currentPosition, dt));
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                Vector3d currentPosition((i + 0.5) * dx, j * dx, (k + 0.5) * dx);
                uy[i][j][k] = getVelocityY(getLastPosition(currentPosition, dt));
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                Vector3d currentPosition((i + 0.5) * dx, (j + 0.5)* dx, k * dx);
                uz[i][j][k] = getVelocityZ(getLastPosition(currentPosition, dt));
            }
        }
    }
}

void
Field::CG_Project(double dt) {
     VectorXd x(Nx*Ny*Nz), b(Nx*Ny*Nz);
     SparseMatrix<double> A(Nx*Ny*Nz, Nx*Ny*Nz);
     A.reserve(allocator);
     double invScale = (rho * dx * dx) / dt;
     for(int k = 0; k < Nz; k++) {
         for(int j = 0; j < Ny; j++) {
             for(int i = 0; i < Nx; i++) {
                 vector<double> D = {1.0, 1.0, 1.0, -1.0, -1.0, -1.0}; 
                 vector<double> F = {static_cast<double>(i < Nx - 1),
                                     static_cast<double>(j < Ny - 1),
                                     static_cast<double>(k < Nz - 1),
                                     static_cast<double>(i > 0),
                                     static_cast<double>(j > 0),
                                     static_cast<double>(k > 0)};
                 vector<double> U = {ux[i + 1][j][k], uy[i][j + 1][k], uz[i][j][k + 1], ux[i][j][k], uy[i][j][k], uz[i][j][k]};
                 double sum_R = 0.0;
                 for(int n = 0; n < 6; n++) {
                     sum_R += invScale * F[n] * D[n] * U[n]/dx;
                 }
                 b[k*(Nx*Ny) + j*Ny + i] = sum_R;
                 if(F[0]) A.insert(k*(Nx*Ny) + j*Ny + i, (k + 0)*(Nx*Ny) + (j + 0)*Ny + (i + 1)) += 1.0;
                 if(F[1]) A.insert(k*(Nx*Ny) + j*Ny + i, (k + 0)*(Nx*Ny) + (j + 1)*Ny + (i + 0)) += 1.0;
                 if(F[2]) A.insert(k*(Nx*Ny) + j*Ny + i, (k + 1)*(Nx*Ny) + (j + 0)*Ny + (i + 0)) += 1.0;
                 if(F[3]) A.insert(k*(Nx*Ny) + j*Ny + i, (k + 0)*(Nx*Ny) + (j + 0)*Ny + (i - 1)) += 1.0;
                 if(F[4]) A.insert(k*(Nx*Ny) + j*Ny + i, (k + 0)*(Nx*Ny) + (j - 1)*Ny + (i + 0)) += 1.0;
                 if(F[5]) A.insert(k*(Nx*Ny) + j*Ny + i, (k - 1)*(Nx*Ny) + (j + 0)*Ny + (i + 0)) += 1.0;
                 A.insert(k*(Nx*Ny) + j*Ny + i, k*(Nx*Ny) + j*Ny + i) -= 6.0;
             }
         } 
     }
     ConjugateGradient<SparseMatrix<double> > cg;
     // cg.setTolerance(1.0e-4);
     cg.setMaxIterations(20);
     cg.compute(A);
     x = cg.solve(b);
     for(int k = 0; k < Nz; k++) {
         for(int j = 0; j < Ny; j++) {
             for(int i = 0; i < Nx; i++) {
                  p[i][j][k] = x[index(i, j, k)];
             }
         } 
     }

     for(int i = 1; i < Nx; i++) {
         for(int j = 0; j < Ny; j++) {
             for(int k = 0; k < Nz; k++) {
                  ux[i][j][k] = ux[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i-1][j][k])/dx);
             }
         } 
     }    
     for(int i = 0; i < Nx; i++) {
         for(int j = 1; j < Ny; j++) {
             for(int k = 0; k < Nz; k++) {
                 uy[i][j][k] = uy[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i][j-1][k])/dx);
             }
         } 
     }    
     for(int i = 0; i < Nx; i++) {
         for(int j = 0; j < Ny; j++) {
             for(int k = 1; k < Nz; k++) {
                 uz[i][j][k] = uz[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i][j][k-1])/dx);
             }
         } 
     }    
}

void
Field::UpdateMarkers(double dt) {
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                sortedMarkersX[index(i, j, k)] += dt * GetVelocity(sortedMarkersX[index(i, j, k)]);
            }
        }
    }
    sortMarkers();
}

Vector3d
Field::getLastPosition(const Vector3d& currentPosition, double dt) {
   Vector3d k1 = GetVelocity(currentPosition); 
   Vector3d k2 = GetVelocity(currentPosition - (dt/2.0) * k1);
   Vector3d k3 = GetVelocity(currentPosition - (dt/2.0) * k2);
   Vector3d k4 = GetVelocity(currentPosition - dt * k3);
   return currentPosition - (dt/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

bool
Field::isInside(const Vector3d& position) const {
    return position.x() > 0.0 && position.x() < Nx * dx 
        && position.y() > 0.0 && position.y() < Ny * dx
        && position.z() > 0.0 && position.z() < Nz * dx;
}

void
Field::SetForce(const Vector3d& force, const Vector3d& position) {
    if(isInside(position)) {
        setForceX(force.x(), position);
        setForceY(force.y(), position);
        setForceZ(force.z(), position);
    }
}

void
Field::setForceX(double fx, const Vector3d& position) {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    y -= dx/2.0;
    z -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1 - 1e-6, z/dx));
    unsigned long i = x;
    unsigned long j = y;
    unsigned long k = z;
    x = x - i;
    y = y - j;
    z = z - k;
    vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
                        (1.0 - x) * (1.0 - y) * z,
                        (1.0 - x) * y * (1.0 - z),
                        x * (1.0 - y) * (1.0 - z),
                        x * y * (1.0 - z),
                        x * (1.0 - y) * z,
                        (1.0 - x) * y * z,
                        x * y * z};
    forcex[i + 0][j + 0][k + 0] = c[0] * fx;
    forcex[i + 0][j + 0][k + 1] = c[1] * fx;
    forcex[i + 0][j + 1][k + 0] = c[2] * fx;
    forcex[i + 1][j + 0][k + 0] = c[3] * fx;
    forcex[i + 1][j + 1][k + 0] = c[4] * fx;
    forcex[i + 1][j + 0][k + 1] = c[5] * fx;
    forcex[i + 0][j + 1][k + 1] = c[6] * fx;
    forcex[i + 1][j + 1][k + 1] = c[7] * fx;
}

void
Field::setForceY(double fy, const Vector3d& position) {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    z -= dx/2.0;
    x -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1 - 1e-6, z/dx));
    unsigned long i = x;
    unsigned long j = y;
    unsigned long k = z;
    x = x - i;
    y = y - j;
    z = z - k;
    vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
                        (1.0 - x) * (1.0 - y) * z,
                        (1.0 - x) * y * (1.0 - z),
                        x * (1.0 - y) * (1.0 - z),
                        x * y * (1.0 - z),
                        x * (1.0 - y) * z,
                        (1.0 - x) * y * z,
                        x * y * z};
    forcey[i + 0][j + 0][k + 0] = c[0] * fy;
    forcey[i + 0][j + 0][k + 1] = c[1] * fy;
    forcey[i + 0][j + 1][k + 0] = c[2] * fy;
    forcey[i + 1][j + 0][k + 0] = c[3] * fy;
    forcey[i + 1][j + 1][k + 0] = c[4] * fy;
    forcey[i + 1][j + 0][k + 1] = c[5] * fy;
    forcey[i + 0][j + 1][k + 1] = c[6] * fy;
    forcey[i + 1][j + 1][k + 1] = c[7] * fy;
}

void
Field::setForceZ(double fz, const Vector3d& position) {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    x -= dx/2.0;
    z -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Ny - 1e-6, z/dx));
    unsigned long i = x;
    unsigned long j = y;
    unsigned long k = z;
    x = x - i;
    y = y - j;
    z = z - k;
    vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
                        (1.0 - x) * (1.0 - y) * z,
                        (1.0 - x) * y * (1.0 - z),
                        x * (1.0 - y) * (1.0 - z),
                        x * y * (1.0 - z),
                        x * (1.0 - y) * z,
                        (1.0 - x) * y * z,
                        x * y * z};
    forcez[i + 0][j + 0][k + 0] = c[0] * fz;
    forcez[i + 0][j + 0][k + 1] = c[1] * fz;
    forcez[i + 0][j + 1][k + 0] = c[2] * fz;
    forcez[i + 1][j + 0][k + 0] = c[3] * fz;
    forcez[i + 1][j + 1][k + 0] = c[4] * fz;
    forcez[i + 1][j + 0][k + 1] = c[5] * fz;
    forcez[i + 0][j + 1][k + 1] = c[6] * fz;
    forcez[i + 1][j + 1][k + 1] = c[7] * fz;
}

Vector3d
Field::GetVelocity(const Vector3d& position) const {
    return Vector3d(getVelocityX(position), getVelocityY(position), getVelocityZ(position));
}

// grid外は境界と同じ値
double
Field::getVelocityX(const Vector3d& position) const {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    y -= dx/2.0;
    z -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1 - 1e-6, z/dx));
    unsigned long i = x;
    unsigned long j = y;
    unsigned long k = z;
    vector<double> f = {ux[i][j][k],
                        ux[i][j][k + 1], ux[i][j + 1][k], ux[i + 1][j][k],
                        ux[i + 1][j + 1][k], ux[i + 1][j][k + 1], ux[i][j + 1][k + 1],
                        ux[i + 1][j + 1][k + 1]};
    x = x - i;
    y = y - j;
    z = z - k;
    vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
                        (1.0 - x) * (1.0 - y) * z,
                        (1.0 - x) * y * (1.0 - z),
                        x * (1.0 - y) * (1.0 - z),
                        x * y * (1.0 - z),
                        x * (1.0 - y) * z,
                        (1.0 - x) * y * z,
                        x * y * z};
    return c[0] * f[0]
           + c[1] * f[1] + c[2] * f[2] + c[3] * f[3]
           + c[4] * f[4] + c[5] * f[5] + c[6] * f[6]
           + c[7] * f[7]; 
}

double
Field::getVelocityY(const Vector3d& position) const {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    z -= dx/2.0;
    x -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1 - 1e-6, z/dx));
    unsigned long i = x;
    unsigned long j = y;
    unsigned long k = z;
    vector<double> f = {uy[i][j][k],
                        uy[i][j][k + 1], uy[i][j + 1][k], uy[i + 1][j][k],
                        uy[i + 1][j + 1][k], uy[i + 1][j][k + 1], uy[i][j + 1][k + 1],
                        uy[i + 1][j + 1][k + 1]};
    x = x - i;
    y = y - j;
    z = z - k;
    vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
                        (1.0 - x) * (1.0 - y) * z,
                        (1.0 - x) * y * (1.0 - z),
                        x * (1.0 - y) * (1.0 - z),
                        x * y * (1.0 - z),
                        x * (1.0 - y) * z,
                        (1.0 - x) * y * z,
                        x * y * z};
    return c[0] * f[0]
           + c[1] * f[1] + c[2] * f[2] + c[3] * f[3]
           + c[4] * f[4] + c[5] * f[5] + c[6] * f[6]
           + c[7] * f[7]; 
}

double
Field::getVelocityZ(const Vector3d& position) const {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    x -= dx/2.0;
    y -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1e-6, z/dx));
    unsigned long i = x;
    unsigned long j = y;
    unsigned long k = z;
    vector<double> f = {uz[i][j][k],
                        uz[i][j][k + 1], uz[i][j + 1][k], uz[i + 1][j][k],
                        uz[i + 1][j + 1][k], uz[i + 1][j][k + 1], uz[i][j + 1][k + 1],
                        uz[i + 1][j + 1][k + 1]};
    x = x - i;
    y = y - j;
    z = z - k;
    vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
                        (1.0 - x) * (1.0 - y) * z,
                        (1.0 - x) * y * (1.0 - z),
                        x * (1.0 - y) * (1.0 - z),
                        x * y * (1.0 - z),
                        x * (1.0 - y) * z,
                        (1.0 - x) * y * z,
                        x * y * z};
    return c[0] * f[0]
           + c[1] * f[1] + c[2] * f[2] + c[3] * f[3]
           + c[4] * f[4] + c[5] * f[5] + c[6] * f[6]
           + c[7] * f[7]; 
}

void
Field::makeBoundary() {
    for(int j = 0; j < Ny; j++) {
        for(int k = 0; k < Nz; k++) {
            ux[0][j][k] = 0.0;
            ux[Nx][j][k] = 0.0;
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int k = 0; k < Nz; k++) {
            uy[i][0][k] = 0.0;
            uy[i][Ny][k] = 0.0;
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            uz[i][j][0] = 0.0;
            uz[i][j][Nz] = 0.0;
        }
    }
}

void
Field::clearForce() {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                forcex[i][j][k] = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                forcey[i][j][k] = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                forcez[i][j][k] = 0.0;
            }
        }
    }
}

void
Field::initVelocity() {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                ux[i][j][k] = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                uy[i][j][k] = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                uz[i][j][k] = 0.0;
            }
        }
    }
}

void
Field::initPressure() {
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                p[i][j][k] = 1.0;
            }
        }
    }
}

void
Field::initMarkers() {
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                sortedMarkersX[index(i, j, k)] = Vector3d((i + 0.5) * dx, (j + 0.5) * dx, (k + 0.5) * dx);
            }
        }
    }
    sortMarkers();
}

void
Field::addGravityForce(double dt) {
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                uz[i][j][k] += dt * g;
            }
        }
    }
}

void
Field::sortMarkers() {
    sort(sortedMarkersX.begin(),
         sortedMarkersX.end(),
         [](const Vector3d& a, const Vector3d& b){return a.x() < b.x();}
        );
}

bool
Field::existsMarker(int cellIndex_x, int cellIndex_y, int cellIndex_z) {
    // Returns an iterator pointing to the first element in the range [first,last) which does not compare less than val.
    auto lower_it = lower_bound(sortedMarkersX.begin(),
                                sortedMarkersX.end(),
                                Vector3d(cellIndex_x * dx, 0, 0),
                                [](const Vector3d& a, const Vector3d& b){return a.x() < b.x();}
                                );
    // Returns an iterator pointing to the first element in the range [first,last) which compares greater than val.
    auto upper_it = upper_bound(sortedMarkersX.begin(),
                                sortedMarkersX.end(),
                                Vector3d((cellIndex_x + 1) * dx, 0, 0),
                                [](const Vector3d& a, const Vector3d& b){return a.x() < b.x();}
                               );
    if(upper_it == ++lower_it) return false;
    for(auto it = lower_it; it < upper_it; it++) {
        if(   cellIndex_y <= it->y() && it->y() <= cellIndex_y + Ny 
           && cellIndex_z <= it->z() && it->z() <= cellIndex_z + Nz) return true;
    }
    return false;
}
