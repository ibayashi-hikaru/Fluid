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
    this->xSwap = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(Nx + 1), std::vector<std::vector<double>>(static_cast<size_t>(Ny), std::vector<double>(static_cast<size_t>(Nz))));
    this->ySwap = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(Nx), std::vector<std::vector<double>>(static_cast<size_t>(Ny + 1), std::vector<double>(static_cast<size_t>(Nz))));
    this->zSwap = std::vector<std::vector<std::vector<double>>>(static_cast<size_t>(Nx), std::vector<std::vector<double>>(static_cast<size_t>(Ny), std::vector<double>(static_cast<size_t>(Nz + 1))));
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
Field::Advect(double dt) {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                Eigen::Vector3d currentPosition(i * dx, (j + 0.5) * dx, (k + 0.5) * dx);
                ux[i][j][k] = getVelocityX(getLastPosition(currentPosition, dt));
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 1; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                Eigen::Vector3d currentPosition((i + 0.5) * dx, j * dx, (k + 0.5) * dx);
                uy[i][j][k] = getVelocityY(getLastPosition(currentPosition, dt));
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 1; k < Nz; k++) {
                Eigen::Vector3d currentPosition((i + 0.5) * dx, (j + 0.5)* dx, k * dx);
                uz[i][j][k] = getVelocityZ(getLastPosition(currentPosition, dt));
            }
        }
    }
}

void
Field::CG_Project(double dt) {
     const unsigned int fullMatrixSize = static_cast<unsigned int>(Nx*Ny*Nz);
     Eigen::VectorXd x(fullMatrixSize), b(fullMatrixSize);
     Eigen::SparseMatrix<double> A(fullMatrixSize, fullMatrixSize);
     tripletList.clear();
     A.reserve(allocator);
     double invScale = (rho * dx * dx) / dt;
     for(size_t k = 0; k < Nz; k++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t i = 0; i < Nx; i++) {
                 std::vector<double> D = {1.0, 1.0, 1.0, -1.0, -1.0, -1.0}; 
                 std::vector<double> F = {static_cast<double>(i < Nx - 1),
                                     static_cast<double>(j < Ny - 1),
                                     static_cast<double>(k < Nz - 1),
                                     static_cast<double>(i > 0),
                                     static_cast<double>(j > 0),
                                     static_cast<double>(k > 0)};
                 std::vector<double> U = {ux[i + 1][j][k], uy[i][j + 1][k], uz[i][j][k + 1], ux[i][j][k], uy[i][j][k], uz[i][j][k]};
                 double sum_R = 0.0;
                 for(size_t n = 0; n < 6; n++) {
                     sum_R += invScale * F[n] * D[n] * U[n]/dx;
                 }
                 b[k*(Nx*Ny) + j*Ny + i] = sum_R;
                 double diagVal = 0.0;
                 if(static_cast<bool>(F[0])) {
                     tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>((k + 0)*(Nx*Ny) + (j + 0)*Ny + (i + 1)), 1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[1])) {
                     tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>((k + 0)*(Nx*Ny) + (j + 1)*Ny + (i + 0)), 1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[2])) {
                     tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>((k + 1)*(Nx*Ny) + (j + 0)*Ny + (i + 0)), 1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[3])) {
                     tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>((k + 0)*(Nx*Ny) + (j + 0)*Ny + (i - 1)), 1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[4])) {
                     tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>((k + 0)*(Nx*Ny) + (j - 1)*Ny + (i + 0)), 1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[5])) {
                     tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>((k - 1)*(Nx*Ny) + (j + 0)*Ny + (i + 0)), 1.0));
                     diagVal -= 1.0;
                 }
                 tripletList.push_back(T(static_cast<const int>(k*(Nx*Ny) + j*Ny + i), static_cast<const int>(k*(Nx*Ny) + j*Ny + i), diagVal));
             }
         } 
     }
     A.setFromTriplets(tripletList.begin(), tripletList.end());
     Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
     // cg.setTolerance(1.0e-4);
     cg.setMaxIterations(20);
     cg.compute(A);
     x = cg.solve(b);
     for(size_t k = 0; k < Nz; k++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t i = 0; i < Nx; i++) {
                  p[i][j][k] = x[index(i, j, k)];
             }
         } 
     }

     for(size_t i = 1; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                  ux[i][j][k] = ux[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i-1][j][k])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 1; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                 uy[i][j][k] = uy[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i][j-1][k])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 1; k < Nz; k++) {
                 uz[i][j][k] = uz[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i][j][k-1])/dx);
             }
         } 
     }    
}

void
Field::CG_ProjectWithMarker(double dt) {
     const unsigned int fullMatrixSize = static_cast<unsigned int>(Nx*Ny*Nz);
     Eigen::VectorXd b(fullMatrixSize);
     Eigen::SparseMatrix<double> A(fullMatrixSize, fullMatrixSize);
     tripletList.clear();
     newTripletList.clear();
     A.reserve(allocator);
     double invScale = (rho * dx * dx) / dt;
     size_t newIndex = 0;
     std::vector<int> m(fullMatrixSize, -1);
     for(size_t k = 0; k < Nz; k++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t i = 0; i < Nx; i++) {
                 std::vector<double> D = {1.0, 1.0, 1.0, -1.0, -1.0, -1.0}; 
                 std::vector<double> F = {static_cast<double>(i < Nx - 1),
                                          static_cast<double>(j < Ny - 1),
                                          static_cast<double>(k < Nz - 1),
                                          static_cast<double>(i > 0),
                                          static_cast<double>(j > 0),
                                          static_cast<double>(k > 0)};
                 std::vector<double> U = {ux[i + 1][j][k],
                                          uy[i][j + 1][k],
                                          uz[i][j][k + 1],
                                          ux[i][j][k],
                                          uy[i][j][k],
                                          uz[i][j][k]};
                 double sum_R = 0.0;
                 for(size_t n = 0; n < 6; n++) {
                     sum_R += invScale * D[n] * U[n]/dx;
                 }
                 b[index(i, j, k)] = sum_R;
                 double diagVal = 0.0;
                 if(static_cast<bool>(F[0])) {
                     tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                             static_cast<const int> (index(i + 1, j + 0, k + 0)),
                                             1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[1])) {
                     tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                             static_cast<const int> (index(i + 0, j + 1, k + 0)),
                                             1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[2])) {
                     tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                             static_cast<const int> (index(i + 0, j + 0, k + 1)),
                                             1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[3])) {
                     tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                             static_cast<const int> (index(i - 1, j + 0, k - 0)),
                                             1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[4])) {
                     tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                             static_cast<const int> (index(i + 0, j - 1, k + 0)),
                                             1.0));
                     diagVal -= 1.0;
                 }
                 if(static_cast<bool>(F[5])) {
                     tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                             static_cast<const int> (index(i + 0, j + 0, k - 1)),
                                             1.0));
                     diagVal -= 1.0;
                 }
                 tripletList.push_back(T(static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                         static_cast<const int> (index(i + 0, j + 0, k + 0)),
                                         diagVal));
                 if(existsMarker(i, j, k)) {
                     m[index(i, j, k)] = static_cast<int>(newIndex++);
                 }
             }
         } 
     }
     A.setFromTriplets(tripletList.begin(), tripletList.end());
     unsigned int shrinkedSize = static_cast<unsigned int>(newIndex);
     Eigen::VectorXd x(shrinkedSize), newb(shrinkedSize);
     Eigen::SparseMatrix<double> newA(shrinkedSize, shrinkedSize);
     for(auto it = tripletList.begin(); it < tripletList.end(); it++) {
        if(m[static_cast<size_t>(it->row())] != -1 && m[static_cast<size_t>(it->col())] != -1) {
            newTripletList.push_back(T(m[static_cast<size_t>(it->row())], m[static_cast<size_t>(it->col())], it->value())); 
        }// definition of row and col need to survey
     }
     newA.setFromTriplets(newTripletList.begin(), newTripletList.end());
     for(size_t i = 0; i < fullMatrixSize; i++) {
         if(m[i] != -1) {
            newb[m[i]] = b[i];
         }
     }
     Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
     // cg.setTolerance(1.0e-4);
     cg.setMaxIterations(20);
     cg.compute(newA);
     x = cg.solve(newb);
     newIndex = 0;
     for(size_t k = 0; k < Nz; k++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t i = 0; i < Nx; i++) {
                 p[i][j][k] = existsMarker(i, j, k) ? x(newIndex++) : 0.0;
             }
         } 
     }

     for(size_t i = 1; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                  ux[i][j][k] = ux[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i-1][j][k])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 1; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                 uy[i][j][k] = uy[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i][j-1][k])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 1; k < Nz; k++) {
                 uz[i][j][k] = uz[i][j][k] - (dt/rho) * ((p[i][j][k] - p[i][j][k-1])/dx);
             }
         } 
     }    
}

void
Field::UpdateMarkers(double dt) {
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                sortedMarkersX[index(i, j, k)] += dt * GetVelocity(sortedMarkersX[index(i, j, k)]);
            }
        }
    }
    sortMarkers();
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
Field::SetForce(const Eigen::Vector3d& force, const Eigen::Vector3d& position) {
    if(isInside(position)) {
        setForceX(force.x(), position);
        setForceY(force.y(), position);
        setForceZ(force.z(), position);
    }
}

void copyVelocity() {
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

void swapVelocity() {

}

void
Field::Extrapolate() {
    saveVelocity();
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
                    if(isnan(uz[i][j][k])) {
                        zSwap[i][j][k] = getAveVelocityZ(i, j, k); 
                        existNan = true;
                    }
                }
            }
        }
        replaceVelocity();
    }
}

double
Field::getAveVelocityX(int i, int j, int k) const {
    std::vector<double> F = {static_cast<double>(i < Nx - 1),
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
                             static_cast<double>(j < Ny - 1),
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
                             static_cast<double>(k < Nz - 1),
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
Field::initMarkers() {
    sortedMarkersX.clear();
    sortedMarkersX.shrink_to_fit();
    waterDrop(Nx*(5.0/6.0), Ny/2.0, Nz/2.0, Nx/6.0);
    storeWater(1.0/3.0);
    sortMarkers();
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

void
Field::addGravityForce(double dt) {
    for(size_t i = 1; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                if(existsMarker(i, j, k) || existsMarker(i + 1, j, k)) {
                    ux[i][j][k] -= dt * g;
                }
            }
        }
    }
}

void
Field::sortMarkers() {
    sort(sortedMarkersX.begin(),
         sortedMarkersX.end(),
         [](const Eigen::Vector3d& a, const Eigen::Vector3d& b){return a.x() < b.x();}
        );
}

bool
Field::existsMarker(size_t cellIndex_x, size_t cellIndex_y, size_t cellIndex_z) {
    if(sortedMarkersX.begin()->x() > (cellIndex_x + 1) * dx) return false;
    if((sortedMarkersX.end() - 1)->x() < cellIndex_x * dx) return false;
    // Returns an iterator pointing to the first element in the range [first,last) which does not compare less than val.
    auto lower_it = lower_bound(sortedMarkersX.begin(),
                                sortedMarkersX.end(),
                                Eigen::Vector3d(cellIndex_x * dx, 0, 0),
                                [](const Eigen::Vector3d& a, const Eigen::Vector3d& b){return a.x() < b.x();}
                                );
    // Returns an iterator pointing to the first element in the range [first,last) which compares greater than val.
    auto upper_it = upper_bound(sortedMarkersX.begin(),
                                sortedMarkersX.end(),
                                Eigen::Vector3d((cellIndex_x + 1) * dx, 0, 0),
                                [](const Eigen::Vector3d& a, const Eigen::Vector3d& b){return a.x() <= b.x();}
                               );
    if(upper_it == lower_it) return false;
    for(auto it = lower_it; it < upper_it; it++) {
        if(   (cellIndex_y) * dx <= it->y() && it->y() <= (cellIndex_y + 1) * dx 
           && (cellIndex_z) * dx <= it->z() && it->z() <= (cellIndex_z + 1) * dx) return true;
    }
    return false;
}

void
Field::CoutDiv() {
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                div[i][j][k] = (ux[i + 1][j][k] - ux[i][j][k])/dx
                             + (uy[i][j + 1][k] - uy[i][j][k])/dx
                             + (uz[i][j][k + 1] - uz[i][j][k])/dx;
            }
        }
    }
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
            }
        }
    }
}
