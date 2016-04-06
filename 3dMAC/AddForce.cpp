#include "Field.h"

void
Field::SetForce(const Eigen::Vector3d& force, const Eigen::Vector3d& position) {
    if(isInside(position)) {
        setForceX(force.x(), position);
        setForceY(force.y(), position);
        setForceZ(force.z(), position);
    }
}

void
Field::AddForce(double dt) {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                if(existsMarker(i - 1, j, k) || existsMarker(i, j, k)) {
                    ux[i][j][k] += dt * forcex[i][j][k];
                } else {
                    ux[i][j][k] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                if(existsMarker(i, j - 1, k) || existsMarker(i, j, k)) {
                    uy[i][j][k] += dt * forcey[i][j][k];
                } else {
                    uy[i][j][k] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                if(existsMarker(i, j, k - 1) || existsMarker(i, j, k)) {
                    uz[i][j][k] += dt * forcez[i][j][k];
                } else {
                    uz[i][j][k] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    addGravityForce(dt);
    clearForce(); 
}

void
Field::setForceX(double fx, const Eigen::Vector3d& position) {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    y -= dx/2.0;
    z -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1 - 1e-6, z/dx));
    unsigned long i = static_cast<unsigned long>(x);
    unsigned long j = static_cast<unsigned long>(y);
    unsigned long k = static_cast<unsigned long>(z);
    x = x - i;
    y = y - j;
    z = z - k;
    std::vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
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
Field::setForceY(double fy, const Eigen::Vector3d& position) {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    z -= dx/2.0;
    x -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1 - 1e-6, z/dx));
    unsigned long i = static_cast<unsigned long>(x);
    unsigned long j = static_cast<unsigned long>(y);
    unsigned long k = static_cast<unsigned long>(z);
    x = x - i;
    y = y - j;
    z = z - k;
    std::vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
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
Field::setForceZ(double fz, const Eigen::Vector3d& position) {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    x -= dx/2.0;
    z -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Ny - 1e-6, z/dx));
    unsigned long i = static_cast<unsigned long>(x);
    unsigned long j = static_cast<unsigned long>(y);
    unsigned long k = static_cast<unsigned long>(z);
    x = x - i;
    y = y - j;
    z = z - k;
    std::vector<double> c = {(1.0 - x) * (1.0 - y) * (1.0 - z),
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

void
Field::clearForce() {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                forcex[i][j][k] = 0.0;
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                forcey[i][j][k] = 0.0;
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                forcez[i][j][k] = 0.0;
            }
        }
    }
}

void
Field::addGravityForce(double dt) {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                if(existsMarker(i - 1, j, k) || existsMarker(i, j, k)) {
                    ux[i][j][k] -= dt * g;
                }
            }
        }
    }
}

