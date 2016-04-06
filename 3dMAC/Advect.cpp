#include "Field.h"

void
Field::Advect(double dt) {
    copyVelocity();
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                Eigen::Vector3d currentPosition(i * dx, (j + 0.5) * dx, (k + 0.5) * dx);
                xSwap[i][j][k] = getVelocityX(getLastPosition(currentPosition, dt));
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                Eigen::Vector3d currentPosition((i + 0.5) * dx, j * dx, (k + 0.5) * dx);
                ySwap[i][j][k] = getVelocityY(getLastPosition(currentPosition, dt));
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                Eigen::Vector3d currentPosition((i + 0.5) * dx, (j + 0.5)* dx, k * dx);
                zSwap[i][j][k] = getVelocityZ(getLastPosition(currentPosition, dt));
            }
        }
    }
    updateVelocityBySwap();
}

