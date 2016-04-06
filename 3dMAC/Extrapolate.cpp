#include "Field.h"

void
Field::Extrapolate() {
    copyVelocity();
    bool existNan = true;
    while(existNan) {
        existNan = false;
        for(size_t i = 1; i < Nx; i++) {
            for(size_t j = 0; j < Ny; j++) {
                for(size_t k = 0; k < Nz; k++) {
                    if(std::isnan(ux[i][j][k])) {
                        xSwap[i][j][k] = getAveVelocityX(i, j, k); 
                        existNan = true;
                    }
                }
            }
        }
        for(size_t i = 0; i < Nx; i++) {
            for(size_t j = 1; j < Ny; j++) {
                for(size_t k = 0; k < Nz; k++) {
                    if(std::isnan(uy[i][j][k])) {
                        ySwap[i][j][k] = getAveVelocityY(i, j, k); 
                        existNan = true;
                    }
                }
            }
        }
        for(size_t i = 0; i < Nx; i++) {
            for(size_t j = 0; j < Ny; j++) {
                for(size_t k = 1; k < Nz; k++) {
                    if(std::isnan(uz[i][j][k])) {
                        zSwap[i][j][k] = getAveVelocityZ(i, j, k); 
                        existNan = true;
                    }
                }
            }
        }
        updateVelocityBySwap();
    }
}

double
Field::getAveVelocityX(int i, int j, int k) const {
    std::vector<double> F = {static_cast<double>(i < Nx),
                             static_cast<double>(j < Ny - 1),
                             static_cast<double>(k < Nz - 1),
                             static_cast<double>(i > 0),
                             static_cast<double>(j > 0),
                             static_cast<double>(k > 0)};
    double over = 0.0;
    double under = 0.0;
    if(F[0] && !std::isnan(ux[i + 1][j + 0][k + 0])) {
        over += ux[i + 1][j + 0][k + 0];
        under += 1.0;
    }
    if(F[1] && !std::isnan(ux[i + 0][j + 1][k + 0])) {
        over += ux[i + 0][j + 1][k + 0];
        under += 1.0;
    }
    if(F[2] && !std::isnan(ux[i + 0][j + 0][k + 1])) {
        over += ux[i + 0][j + 0][k + 1];
        under += 1.0;
    }
    if(F[3] && !std::isnan(ux[i - 1][j + 0][k + 0])) {
        over += ux[i - 1][j + 0][k + 0];
        under += 1.0;
    }
    if(F[4] && !std::isnan(ux[i + 0][j - 1][k + 0])) {
        over += ux[i + 0][j - 1][k + 0];
        under += 1.0;
    }
    if(F[5] && !std::isnan(ux[i + 0][j + 0][k - 1])) {
        over += ux[i + 0][j + 0][k - 1];
        under += 1.0;
    }
    if(under != 0.0) {
        return over/under;
    } else {
        return std::numeric_limits<double>::quiet_NaN(); 
    }
}

double
Field::getAveVelocityY(int i, int j, int k) const {
    std::vector<double> F = {static_cast<double>(i < Nx - 1),
                             static_cast<double>(j < Ny),
                             static_cast<double>(k < Nz - 1),
                             static_cast<double>(i > 0),
                             static_cast<double>(j > 0),
                             static_cast<double>(k > 0)};
    double over = 0.0;
    double under = 0.0;
    if(F[0] && !std::isnan(uy[i + 1][j + 0][k + 0])) {
        over += uy[i + 1][j + 0][k + 0];
        under += 1.0;
    }
    if(F[1] && !std::isnan(uy[i + 0][j + 1][k + 0])) {
        over += uy[i + 0][j + 1][k + 0];
        under += 1.0;
    }
    if(F[2] && !std::isnan(uy[i + 0][j + 0][k + 1])) {
        over += uy[i + 0][j + 0][k + 1];
        under += 1.0;
    }
    if(F[3] && !std::isnan(uy[i - 1][j + 0][k + 0])) {
        over += uy[i - 1][j + 0][k + 0];
        under += 1.0;
    }
    if(F[4] && !std::isnan(uy[i + 0][j - 1][k + 0])) {
        over += uy[i + 0][j - 1][k + 0];
        under += 1.0;
    }
    if(F[5] && !std::isnan(uy[i + 0][j + 0][k - 1])) {
        over += uy[i + 0][j + 0][k - 1];
        under += 1.0;
    }
    if(under != 0.0) {
        return over/under;
    } else {
        return std::numeric_limits<double>::quiet_NaN(); 
    }
}

double
Field::getAveVelocityZ(int i, int j, int k) const {
    std::vector<double> F = {static_cast<double>(i < Nx - 1),
                             static_cast<double>(j < Ny - 1),
                             static_cast<double>(k < Nz),
                             static_cast<double>(i > 0),
                             static_cast<double>(j > 0),
                             static_cast<double>(k > 0)};
    double over = 0.0;
    double under = 0.0;
    if(F[0] && !std::isnan(uz[i + 1][j + 0][k + 0])) {
        over += uz[i + 1][j + 0][k + 0];
        under += 1.0;
    }
    if(F[1] && !std::isnan(uz[i + 0][j + 1][k + 0])) {
        over += uz[i + 0][j + 1][k + 0];
        under += 1.0;
    }
    if(F[2] && !std::isnan(uz[i + 0][j + 0][k + 1])) {
        over += uz[i + 0][j + 0][k + 1];
        under += 1.0;
    }
    if(F[3] && !std::isnan(uz[i - 1][j + 0][k + 0])) {
        over += uz[i - 1][j + 0][k + 0];
        under += 1.0;
    }
    if(F[4] && !std::isnan(uz[i + 0][j - 1][k + 0])) {
        over += uz[i + 0][j - 1][k + 0];
        under += 1.0;
    }
    if(F[5] && !std::isnan(uz[i + 0][j + 0][k - 1])) {
        over += uz[i + 0][j + 0][k - 1];
        under += 1.0;
    }
    if(under != 0.0) {
        return over/under;
    } else {
        return std::numeric_limits<double>::quiet_NaN(); 
    }
}

