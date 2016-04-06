#include "Field.h"
Field::Field(unsigned long gridNum) {
    this->Nx = gridNum;
    this->Ny = gridNum;
    this->Nz = gridNum;
    this->dx = 1.0;
    this->div = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum))));
    this->p = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum))));

    this->ux = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum + 1), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum))));
    this->uy = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum + 1), std::vector<double>(static_cast<size_t>(gridNum))));
    this->uz = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum + 1))));
    this->xSwap = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum + 1), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum))));
    this->ySwap = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum + 1), std::vector<double>(static_cast<size_t>(gridNum))));
    this->zSwap = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum + 1))));

    this->forcex = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum + 1), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum))));
    this->forcey = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum + 1), std::vector<double>(static_cast<size_t>(gridNum))));
    this->forcez = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(gridNum), std::vector<std::vector<double>>(static_cast<size_t>(gridNum), std::vector<double>(static_cast<size_t>(gridNum + 1))));
    this->sortedMarkersX.resize(static_cast<size_t>(gridNum * gridNum * gridNum));
    allocator = std::vector<int>(static_cast<size_t>(Nx*Ny*Nz), 7);
    allocator.at(0) = 4; 
    allocator.at(1) = 5; 
    allocator.at(2) = 6; 
    allocator.at(static_cast<size_t>(Nx*Ny*Nz - 3)) = 6; 
    allocator.at(static_cast<size_t>(Nx*Ny*Nz - 2)) = 5; 
    allocator.at(static_cast<size_t>(Nx*Ny*Nz - 1)) = 4; 
    tripletList.reserve(static_cast<size_t>(Nx*Ny*Nz));
    newTripletList.reserve(static_cast<size_t>(Nx*Ny*Nz));
}

Field& Field::operator=(const Field& other) {
    this->Nx = other.Nx;
    this->Ny = other.Ny;
    this->Nz = other.Nz;
    this->dx = 1.0;
    this->div = other.div;
    this->p = other.p;
    this->ux = other.ux;
    this->uy = other.uy;
    this->uz = other.uz;
    this->xSwap = other.xSwap;
    this->ySwap = other.ySwap;
    this->zSwap = other.zSwap;
    this->forcex = other.forcex;
    this->forcey = other.forcey;
    this->forcez = other.forcez;
    std::copy(other.sortedMarkersX.begin(), other.sortedMarkersX.end(), back_inserter(this->sortedMarkersX));
    std::copy(other.allocator.begin(), other.allocator.end(), back_inserter(this->allocator));
    std::copy(other.tripletList.begin(), other.tripletList.end(), back_inserter(this->tripletList));
    std::copy(other.newTripletList.begin(), other.newTripletList.end(), back_inserter(this->newTripletList));
    return *this;
}

void
Field::Init() {
    makeBoundary();
    clearForce();
    initVelocity();
    initPressure();
    initMarkers();
}

Eigen::Vector3d
Field::getLastPosition(const Eigen::Vector3d& currentPosition, double dt) {
   Eigen::Vector3d k1 = GetVelocity(currentPosition); 
   Eigen::Vector3d k2 = GetVelocity(currentPosition - (dt/2.0) * k1);
   Eigen::Vector3d k3 = GetVelocity(currentPosition - (dt/2.0) * k2);
   Eigen::Vector3d k4 = GetVelocity(currentPosition - dt * k3);
   return currentPosition - (dt/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

bool
Field::isInside(const Eigen::Vector3d& position) const {
    return position.x() > 0.0 && position.x() < Nx * dx 
        && position.y() > 0.0 && position.y() < Ny * dx
        && position.z() > 0.0 && position.z() < Nz * dx;
}


void
Field::copyVelocity() {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                xSwap[i][j][k] = ux[i][j][k];
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                ySwap[i][j][k] = uy[i][j][k];
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                zSwap[i][j][k] = uz[i][j][k];
            }
        }
    }
}

void
Field::updateVelocityBySwap() {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                ux[i][j][k] = xSwap[i][j][k];
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                uy[i][j][k] = ySwap[i][j][k];
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                uz[i][j][k] = zSwap[i][j][k];
            }
        }
    }

}


Eigen::Vector3d
Field::GetVelocity(const Eigen::Vector3d& position) const {
    return Eigen::Vector3d(getVelocityX(position), getVelocityY(position), getVelocityZ(position));
}

// grid外は境界と同じ値
double
Field::getVelocityX(const Eigen::Vector3d& position) const {
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
    std::vector<double> f = {ux[i][j][k],
                        ux[i][j][k + 1], ux[i][j + 1][k], ux[i + 1][j][k],
                        ux[i + 1][j + 1][k], ux[i + 1][j][k + 1], ux[i][j + 1][k + 1],
                        ux[i + 1][j + 1][k + 1]};
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
    return c[0] * f[0]
           + c[1] * f[1] + c[2] * f[2] + c[3] * f[3]
           + c[4] * f[4] + c[5] * f[5] + c[6] * f[6]
           + c[7] * f[7]; 
}

double
Field::getVelocityY(const Eigen::Vector3d& position) const {
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
    std::vector<double> f = {uy[i][j][k],
                        uy[i][j][k + 1], uy[i][j + 1][k], uy[i + 1][j][k],
                        uy[i + 1][j + 1][k], uy[i + 1][j][k + 1], uy[i][j + 1][k + 1],
                        uy[i + 1][j + 1][k + 1]};
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
    return c[0] * f[0]
           + c[1] * f[1] + c[2] * f[2] + c[3] * f[3]
           + c[4] * f[4] + c[5] * f[5] + c[6] * f[6]
           + c[7] * f[7]; 
}

double
Field::getVelocityZ(const Eigen::Vector3d& position) const {
    double x = position.x();
    double y = position.y();
    double z = position.z();
    x -= dx/2.0;
    y -= dx/2.0;
    x = fmax(0.0, fmin(Nx - 1 - 1e-6, x/dx));
    y = fmax(0.0, fmin(Ny - 1 - 1e-6, y/dx));
    z = fmax(0.0, fmin(Nz - 1e-6, z/dx));
    unsigned long i = static_cast<unsigned long>(x);
    unsigned long j = static_cast<unsigned long>(y);
    unsigned long k = static_cast<unsigned long>(z);
    std::vector<double> f = {uz[i][j][k],
                        uz[i][j][k + 1], uz[i][j + 1][k], uz[i + 1][j][k],
                        uz[i + 1][j + 1][k], uz[i + 1][j][k + 1], uz[i][j + 1][k + 1],
                        uz[i + 1][j + 1][k + 1]};
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
    return c[0] * f[0]
           + c[1] * f[1] + c[2] * f[2] + c[3] * f[3]
           + c[4] * f[4] + c[5] * f[5] + c[6] * f[6]
           + c[7] * f[7]; 
}

void
Field::makeBoundary() {
    for(size_t j = 0; j < Ny; j++) {
        for(size_t k = 0; k < Nz; k++) {
            ux[0][j][k] = 0.0;
            ux[Nx][j][k] = 0.0;
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t k = 0; k < Nz; k++) {
            uy[i][0][k] = 0.0;
            uy[i][Ny][k] = 0.0;
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            uz[i][j][0] = 0.0;
            uz[i][j][Nz] = 0.0;
        }
    }
}


void
Field::initVelocity() {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                ux[i][j][k] = 0.0;
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                uy[i][j][k] = 0.0;
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                uz[i][j][k] = 0.0;
            }
        }
    }
}

void
Field::initPressure() {
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                p[i][j][k] = 1.0;
            }
        }
    }
}


void
Field::waterDrop(double x, double y, double z, double radius) {
    Eigen::Vector3d center(x, y, z);
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                if((centerPosition(i, j, k) - center).norm() < radius) {
                    sortedMarkersX.push_back(centerPosition(i, j, k));
                }
            }
        }
    }
}

void
Field::storeWater(double rate) {
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                if(centerPosition(i, j, k).x() < Nx*dx*rate) {
                    sortedMarkersX.push_back(centerPosition(i, j, k));
                }
            }
        }
    }
}

