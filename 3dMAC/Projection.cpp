#include "Field.h"
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
                  p[index(i, j, k)] = x[index(i, j, k)];
             }
         } 
     }

     for(size_t i = 1; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                  ux[i][j][k] = ux[i][j][k] - (dt/rho) * ((p[index(i, j, k)] - p[index(i-1, j, k)])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 1; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                 uy[i][j][k] = uy[i][j][k] - (dt/rho) * ((p[index(i, j, k)] - p[index(i, j-1, k)])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 1; k < Nz; k++) {
                 uz[i][j][k] = uz[i][j][k] - (dt/rho) * ((p[index(i, j, k)] - p[index(i, j, k-1)])/dx);
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
                 p[index(i, j, k)] = existsMarker(i, j, k) ? x(newIndex++) : 0.0;
             }
         } 
     }

     for(size_t i = 1; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                  ux[i][j][k] = ux[i][j][k] - (dt/rho) * ((p[index(i, j, k)] - p[index(i-1, j, k)])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 1; j < Ny; j++) {
             for(size_t k = 0; k < Nz; k++) {
                 uy[i][j][k] = uy[i][j][k] - (dt/rho) * ((p[index(i, j, k)] - p[index(i, j-1, k)])/dx);
             }
         } 
     }    
     for(size_t i = 0; i < Nx; i++) {
         for(size_t j = 0; j < Ny; j++) {
             for(size_t k = 1; k < Nz; k++) {
                 uz[i][j][k] = uz[i][j][k] - (dt/rho) * ((p[index(i, j, k)] - p[index(i, j, k-1)])/dx);
             }
         } 
     }    
}

