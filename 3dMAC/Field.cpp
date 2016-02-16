#include "Field.h"
using namespace std;

void
Field::Init() {
    makeBoundary();
    clearForce();
    initVelocity();
    initPressure();
}

void
Field::AddForce(double dt) {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                ux.at(i).at(j).at(k) += dt * forcex.at(i).at(j).at(k);
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                uy.at(i).at(j).at(k) += dt * forcey.at(i).at(j).at(k);
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                uz.at(i).at(j).at(k) += dt * forcez.at(i).at(j).at(k);
            }
        }
    }
    clearForce(); 
}

void
Field::Advect(double dt) {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                Vector3d currentPosition(i * dx, (j + 0.5) * dx, (k + 0.5) * dx);
                Vector3d lastPosition = getLastPosition(currentPosition, dt);
                ux.at(i).at(j).at(k) = getVelocityX(lastPosition);
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                Vector3d currentPosition((i + 0.5) * dx, j * dx, (k + 0.5) * dx);
                Vector3d lastPosition = getLastPosition(currentPosition, dt);
                uy.at(i).at(j).at(k) = getVelocityY(lastPosition);
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                Vector3d currentPosition((i + 0.5) * dx, (j + 0.5)* dx, k * dx);
                Vector3d lastPosition = getLastPosition(currentPosition, dt);
                uz.at(i).at(j).at(k) = getVelocityY(lastPosition);
            }
        }
    }
}

void
Field::GS_Project(double dt) {
    // double scale = dt / (rho * dx * dx);
    // double eps = 1.0e-1;
    // double err;
    // double sorParam = 1.8;
    // do {
    //     err = 0.0;

    //     for (int j = 0; j < Ny; j++) {
    //         for (int i = 0; i < Nx; i++) {
    //             vector<double> D = {1.0, 1.0, -1.0, -1.0}; 
    //             vector<double> F = {static_cast<double>(i < Nx - 1),
    //                                 static_cast<double>(j < Ny - 1),
    //                                 static_cast<double>(i > 0),
    //                                 static_cast<double>(j > 0)};
    //             vector<double> P = {F[0] ? p[i + 1][j + 0] : 0.0,
    //                                 F[1] ? p[i + 0][j + 1] : 0.0,
    //                                 F[2] ? p[i - 1][j + 0] : 0.0,
    //                                 F[3] ? p[i + 0][j - 1] : 0.0 };
    //             vector<double> U = {ux[i + 1][j], uy[i][j + 1], ux[i][j], uy[i][j]};
    //             
    //             double det = 0.0;
    //             double sum_L = 0.0;
    //             double sum_R = 0.0;
    //             for(int n = 0; n < 4; n++) {
    //                 det += F[n] * scale;
    //                 sum_L += F[n] * P[n] * scale;
    //                 sum_R += F[n] * D[n] * U[n]/dx;
    //             }
    //             err = fmax(err, fabs(det*p[i][j] - sum_L + sum_R)); 
    //             p[i][j] = (1.0 - sorParam) * p[i][j] + sorParam * (sum_L - sum_R)/det;
    //         } 
    //     }
    // } while(eps < err);

    // for(int i = 1; i < Nx; i++) {
    //     for(int j = 0; j < Ny; j++) {
    //          ux.at(i).at(j) = ux.at(i).at(j) - (dt/rho) * ((p.at(i).at(j) - p.at(i-1).at(j))/dx);
    //     } 
    // }    
    // for(int i = 0; i < Nx; i++) {
    //     for(int j = 1; j < Ny; j++) {
    //          uy.at(i).at(j) = uy.at(i).at(j) - (dt/rho) * ((p.at(i).at(j) - p.at(i).at(j-1))/dx);
    //     } 
    // }    
}

void
Field::CG_Project(double dt) {
    // VectorXd x(Nx*Ny), b(Nx*Ny);
    // SparseMatrix<double> A(Nx*Ny, Nx*Ny);
    // A.reserve(static_cast<int> (5*Nx*Ny));
    // double invScale = (rho * dx * dx) / dt;
    // for(int i = 0; i < Nx; i++) {
    //     for(int j = 0; j < Ny; j++) {
    //         vector<double> D = {1.0, 1.0, -1.0, -1.0}; 
    //         vector<double> F = {static_cast<double>(i < Nx - 1),
    //                             static_cast<double>(j < Ny - 1),
    //                             static_cast<double>(i > 0),
    //                             static_cast<double>(j > 0)};
    //         vector<double> U = {ux[i + 1][j], uy[i][j + 1], ux[i][j], uy[i][j]};
    //         double sum_R = 0.0;
    //         for(int n = 0; n < 4; n++) {
    //             sum_R += invScale * F[n] * D[n] * U[n]/dx;
    //         }
    //         b[i*Ny + j] = sum_R;
    //         if(F.at(0)) A.insert(i*Ny + j, (i + 1)*Ny + (j + 0)) += 1.0;
    //         if(F.at(1)) A.insert(i*Ny + j, (i + 0)*Ny + (j + 1)) += 1.0;
    //         if(F.at(2)) A.insert(i*Ny + j, (i - 1)*Ny + (j + 0)) += 1.0;
    //         if(F.at(3)) A.insert(i*Ny + j, (i + 0)*Ny + (j - 1)) += 1.0;
    //         A.insert(i*Ny + j, (i + 0)*Ny + (j + 0)) -= 4.0;
    //     } 
    // }
    // ConjugateGradient<SparseMatrix<double> > cg;
    // cg.setTolerance(1.0e-1);
    // cg.compute(A);
    // x = cg.solve(b);
    // for(int i = 0; i < Nx; i++) {
    //     for(int j = 0; j < Ny; j++) {
    //         p[i][j] = x[i*Ny + j];
    //     } 
    // }

    // for(int i = 1; i < Nx; i++) {
    //     for(int j = 0; j < Ny; j++) {
    //          ux.at(i).at(j) = ux.at(i).at(j) - (dt/rho) * ((p.at(i).at(j) - p.at(i-1).at(j))/dx);
    //     } 
    // }    
    // for(int i = 0; i < Nx; i++) {
    //     for(int j = 1; j < Ny; j++) {
    //          uy.at(i).at(j) = uy.at(i).at(j) - (dt/rho) * ((p.at(i).at(j) - p.at(i).at(j-1))/dx);
    //     } 
    // }    
}

Vector3d
Field::getLastPosition(Vector3d currentPosition, double dt) {
   Vector3d k1 = GetVelocity(currentPosition); 
   Vector3d k2 = GetVelocity(currentPosition - (dt/2.0) * k1);
   Vector3d k3 = GetVelocity(currentPosition - (dt/2.0) * k2);
   Vector3d k4 = GetVelocity(currentPosition - dt * k3);
   return currentPosition - (dt/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

bool
Field::isInside(Vector3d position) const {
    return position.x() >= 0.0 && position.x() <= Nx * dx 
        && position.y() >= 0.0 && position.y() <= Ny * dx
        && position.z() >= 0.0 && position.z() <= Nz * dx;
}

void
Field::SetForce(Vector3d force, Vector3d position) {
    if(isInside(position)) {
        setForceX(force.x(), position);
        setForceY(force.y(), position);
        setForceZ(force.z(), position);
    }
}

void
Field::setForceX(double fx, Vector3d position) {
    position.y() -= dx/2.0;
    position.z() -= dx/2.0;
    position.x() = fmax(0.0, fmin(Nx - 1e-6, position.x()/dx));
    position.y() = fmax(0.0, fmin(Ny - 1 - 1e-6, position.y()/dx));
    position.z() = fmax(0.0, fmin(Nz - 1 - 1e-6, position.z()/dx));
    unsigned long i = position.x();
    unsigned long j = position.y();
    unsigned long k = position.z();
    position.x() = position.x() - i;
    position.y() = position.y() - j;
    position.z() = position.z() - k;
    vector<double> c = {(1.0 - position.x()) * (1.0 - position.y()) * (1.0 - position.z()),
                        (1.0 - position.x()) * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * (1.0 - position.z()),
                        position.x() * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * position.z(),
                        position.x() * position.y() * position.z()};
    forcex.at(i + 0).at(j + 0).at(k + 0) = c.at(0) * fx;
    forcex.at(i + 0).at(j + 0).at(k + 1) = c.at(1) * fx;
    forcex.at(i + 0).at(j + 1).at(k + 0) = c.at(2) * fx;
    forcex.at(i + 1).at(j + 0).at(k + 0) = c.at(3) * fx;
    forcex.at(i + 1).at(j + 1).at(k + 0) = c.at(4) * fx;
    forcex.at(i + 1).at(j + 0).at(k + 1) = c.at(5) * fx;
    forcex.at(i + 0).at(j + 1).at(k + 1) = c.at(6) * fx;
    forcex.at(i + 1).at(j + 1).at(k + 1) = c.at(7) * fx;
}

void
Field::setForceY(double fy, Vector3d position) {
    position.z() -= dx/2.0;
    position.x() -= dx/2.0;
    position.x() = fmax(0.0, fmin(Nx - 1 - 1e-6, position.x()/dx));
    position.y() = fmax(0.0, fmin(Ny - 1e-6, position.y()/dx));
    position.z() = fmax(0.0, fmin(Nz - 1 - 1e-6, position.z()/dx));
    unsigned long i = position.x();
    unsigned long j = position.y();
    unsigned long k = position.z();
    position.x() = position.x() - i;
    position.y() = position.y() - j;
    position.z() = position.z() - k;
    vector<double> c = {(1.0 - position.x()) * (1.0 - position.y()) * (1.0 - position.z()),
                        (1.0 - position.x()) * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * (1.0 - position.z()),
                        position.x() * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * position.z(),
                        position.x() * position.y() * position.z()};
    forcex.at(i + 0).at(j + 0).at(k + 0) = c.at(0) * fy;
    forcex.at(i + 0).at(j + 0).at(k + 1) = c.at(1) * fy;
    forcex.at(i + 0).at(j + 1).at(k + 0) = c.at(2) * fy;
    forcex.at(i + 1).at(j + 0).at(k + 0) = c.at(3) * fy;
    forcex.at(i + 1).at(j + 1).at(k + 0) = c.at(4) * fy;
    forcex.at(i + 1).at(j + 0).at(k + 1) = c.at(5) * fy;
    forcex.at(i + 0).at(j + 1).at(k + 1) = c.at(6) * fy;
    forcex.at(i + 1).at(j + 1).at(k + 1) = c.at(7) * fy;
}

void
Field::setForceZ(double fz, Vector3d position) {
    position.x() -= dx/2.0;
    position.z() -= dx/2.0;
    position.x() = fmax(0.0, fmin(Nx - 1 - 1e-6, position.x()/dx));
    position.y() = fmax(0.0, fmin(Ny - 1 - 1e-6, position.y()/dx));
    position.z() = fmax(0.0, fmin(Ny - 1e-6, position.z()/dx));
    unsigned long i = position.x();
    unsigned long j = position.y();
    unsigned long k = position.z();
    position.x() = position.x() - i;
    position.y() = position.y() - j;
    position.z() = position.z() - k;
    vector<double> c = {(1.0 - position.x()) * (1.0 - position.y()) * (1.0 - position.z()),
                        (1.0 - position.x()) * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * (1.0 - position.z()),
                        position.x() * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * position.z(),
                        position.x() * position.y() * position.z()};
    forcex.at(i + 0).at(j + 0).at(k + 0) = c.at(0) * fz;
    forcex.at(i + 0).at(j + 0).at(k + 1) = c.at(1) * fz;
    forcex.at(i + 0).at(j + 1).at(k + 0) = c.at(2) * fz;
    forcex.at(i + 1).at(j + 0).at(k + 0) = c.at(3) * fz;
    forcex.at(i + 1).at(j + 1).at(k + 0) = c.at(4) * fz;
    forcex.at(i + 1).at(j + 0).at(k + 1) = c.at(5) * fz;
    forcex.at(i + 0).at(j + 1).at(k + 1) = c.at(6) * fz;
    forcex.at(i + 1).at(j + 1).at(k + 1) = c.at(7) * fz;
}

Vector3d
Field::GetVelocity(Vector3d position) const {
    return Vector3d(getVelocityX(position), getVelocityY(position), getVelocityZ(position));
}

// grid外は境界と同じ値
double
Field::getVelocityX(Vector3d position) const {
    position.y() -= dx/2.0;
    position.z() -= dx/2.0;
    position.x() = fmax(0.0, fmin(Nx - 1e-6, position.x()/dx));
    position.y() = fmax(0.0, fmin(Ny - 1 - 1e-6, position.y()/dx));
    position.z() = fmax(0.0, fmin(Nz - 1 - 1e-6, position.z()/dx));
    unsigned long i = position.x();
    unsigned long j = position.y();
    unsigned long k = position.z();
    vector<double> f = {ux.at(i).at(j).at(k),
                        ux.at(i).at(j).at(k + 1), ux.at(i).at(j + 1).at(k), ux.at(i + 1).at(j).at(k),
                        ux.at(i + 1).at(j + 1).at(k), ux.at(i + 1).at(j).at(k + 1), ux.at(i).at(j + 1).at(k + 1),
                        ux.at(i + 1).at(j + 1).at(k + 1)};
    position.x() = position.x() - i;
    position.y() = position.y() - j;
    position.z() = position.z() - k;
    vector<double> c = {(1.0 - position.x()) * (1.0 - position.y()) * (1.0 - position.z()),
                        (1.0 - position.x()) * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * (1.0 - position.z()),
                        position.x() * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * position.z(),
                        position.x() * position.y() * position.z()};
    return c.at(0) * f.at(0)
           + c.at(1) * f.at(1) + c.at(2) * f.at(2) + c.at(3) * f.at(3)
           + c.at(4) * f.at(4) + c.at(5) * f.at(5) + c.at(6) * f.at(6)
           + c.at(7) * f.at(7); 
}

double
Field::getVelocityY(Vector3d position) const {
    position.z() -= dx/2.0;
    position.x() -= dx/2.0;
    position.x() = fmax(0.0, fmin(Nx - 1 - 1e-6, position.x()/dx));
    position.y() = fmax(0.0, fmin(Ny - 1e-6, position.y()/dx));
    position.z() = fmax(0.0, fmin(Nz - 1 - 1e-6, position.z()/dx));
    unsigned long i = position.x();
    unsigned long j = position.y();
    unsigned long k = position.z();
    vector<double> f = {uy.at(i).at(j).at(k),
                        uy.at(i).at(j).at(k + 1), uy.at(i).at(j + 1).at(k), uy.at(i + 1).at(j).at(k),
                        uy.at(i + 1).at(j + 1).at(k), uy.at(i + 1).at(j).at(k + 1), uy.at(i).at(j + 1).at(k + 1),
                        uy.at(i + 1).at(j + 1).at(k + 1)};
    position.x() = position.x() - i;
    position.y() = position.y() - j;
    position.z() = position.z() - k;
    vector<double> c = {(1.0 - position.x()) * (1.0 - position.y()) * (1.0 - position.z()),
                        (1.0 - position.x()) * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * (1.0 - position.z()),
                        position.x() * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * position.z(),
                        position.x() * position.y() * position.z()};
    return c.at(0) * f.at(0)
           + c.at(1) * f.at(1) + c.at(2) * f.at(2) + c.at(3) * f.at(3)
           + c.at(4) * f.at(4) + c.at(5) * f.at(5) + c.at(6) * f.at(6)
           + c.at(7) * f.at(7); 
}

double
Field::getVelocityZ(Vector3d position) const {
    position.x() -= dx/2.0;
    position.y() -= dx/2.0;
    position.x() = fmax(0.0, fmin(Nx - 1 - 1e-6, position.x()/dx));
    position.y() = fmax(0.0, fmin(Ny - 1 - 1e-6, position.y()/dx));
    position.z() = fmax(0.0, fmin(Nz - 1e-6, position.z()/dx));
    unsigned long i = position.x();
    unsigned long j = position.y();
    unsigned long k = position.z();
    vector<double> f = {uz.at(i).at(j).at(k),
                        uz.at(i).at(j).at(k + 1), uz.at(i).at(j + 1).at(k), uz.at(i + 1).at(j).at(k),
                        uz.at(i + 1).at(j + 1).at(k), uz.at(i + 1).at(j).at(k + 1), uz.at(i).at(j + 1).at(k + 1),
                        uz.at(i + 1).at(j + 1).at(k + 1)};
    position.x() = position.x() - i;
    position.y() = position.y() - j;
    position.z() = position.z() - k;
    vector<double> c = {(1.0 - position.x()) * (1.0 - position.y()) * (1.0 - position.z()),
                        (1.0 - position.x()) * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * (1.0 - position.z()),
                        position.x() * position.y() * (1.0 - position.z()),
                        position.x() * (1.0 - position.y()) * position.z(),
                        (1.0 - position.x()) * position.y() * position.z(),
                        position.x() * position.y() * position.z()};
    return c.at(0) * f.at(0)
           + c.at(1) * f.at(1) + c.at(2) * f.at(2) + c.at(3) * f.at(3)
           + c.at(4) * f.at(4) + c.at(5) * f.at(5) + c.at(6) * f.at(6)
           + c.at(7) * f.at(7); 
}

void
Field::makeBoundary() {
    for(int j = 0; j < Ny; j++) {
        for(int k = 0; k < Nz; k++) {
            ux.at(0).at(j).at(k) = 0.0;
            ux.at(Nx).at(j).at(k) = 0.0;
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int k = 0; k < Nz; k++) {
            uy.at(i).at(0).at(k) = 0.0;
            uy.at(i).at(Ny).at(k) = 0.0;
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            uz.at(i).at(j).at(0) = 0.0;
            uz.at(i).at(j).at(Nz) = 0.0;
        }
    }
}

void
Field::clearForce() {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                forcex.at(i).at(j).at(k) = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                forcey.at(i).at(j).at(k) = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                forcez.at(i).at(j).at(k) = 0.0;
            }
        }
    }
}
void
Field::initVelocity() {
    for(int i = 1; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                ux.at(i).at(j).at(k) = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 1; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                uy.at(i).at(j).at(k) = 0.0;
            }
        }
    }
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 1; k < Nz; k++) {
                uz.at(i).at(j).at(k) = 0.0;
            }
        }
    }
}

void
Field::initPressure() {
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < Nz; k++) {
                uz.at(i).at(j).at(k) = 1.0;
            }
        }
    }
}
