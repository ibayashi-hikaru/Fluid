#ifndef ST_FIELD_H_INCLUDED
#define ST_FIELD_H_INCLUDED
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <vector>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
using namespace std;
using namespace Eigen;
using T = Eigen::Triplet<double>;
/*
    Class's description
*/

class Field {
    public:
       Field(unsigned long gridNum) {
            this->Nx = gridNum;
            this->Ny = gridNum;
            this->Nz = gridNum;
            this->dx = 1.0;
            this->div = vector<vector<vector<double>>>(gridNum, vector<vector<double>>(gridNum, vector<double>(gridNum)));
            this->p = vector<vector<vector<double>>>(gridNum, vector<vector<double>>(gridNum, vector<double>(gridNum)));
            this->ux = vector<vector<vector<double>>>(gridNum + 1, vector<vector<double>>(gridNum, vector<double>(gridNum)));
            this->uy = vector<vector<vector<double>>>(gridNum, vector<vector<double>>(gridNum + 1, vector<double>(gridNum)));
            this->uz = vector<vector<vector<double>>>(gridNum, vector<vector<double>>(gridNum, vector<double>(gridNum + 1)));
            this->forcex = vector<vector<vector<double>>>(gridNum + 1, vector<vector<double>>(gridNum, vector<double>(gridNum)));
            this->forcey = vector<vector<vector<double>>>(gridNum, vector<vector<double>>(gridNum + 1, vector<double>(gridNum)));
            this->forcez = vector<vector<vector<double>>>(gridNum, vector<vector<double>>(gridNum, vector<double>(gridNum + 1)));
            this->sortedMarkersX.resize(gridNum * gridNum * gridNum);
            allocator = vector<int>(Nx*Ny*Nz, 7);
            allocator.at(0) = 4; 
            allocator.at(1) = 5; 
            allocator.at(2) = 6; 
            allocator.at(Nx*Ny*Nz - 3) = 6; 
            allocator.at(Nx*Ny*Nz - 2) = 5; 
            allocator.at(Nx*Ny*Nz - 1) = 4; 
            tripletList.reserve(Nx*Ny*Nz);
            newTripletList.reserve(Nx*Ny*Nz);
       }
       int GridNum() const {return Nx;};
       double Dx() const {return dx;};
       void Init();
       void Advect(double dt);
       void AddForce(double dt);      
       void CG_Project(double dt);
       void CG_ProjectWithMarker(double dt);
       void UpdateMarkers(double dt);
       void SetForce(const Vector3d& force, const Vector3d& position);
       Vector3d GetVelocity(const Vector3d& position) const;
       vector< Vector3d> sortedMarkersX;
       void CoutDiv();
    private:
       const double rho = 1.0;
       const double g = 0.00098;
       unsigned long Nx, Ny, Nz;
       double dx;
       vector< vector< vector< double>>> div;
       vector< vector< vector< double>>> p;
       vector< vector< vector< double>>> ux;
       vector< vector< vector< double>>> uy;
       vector< vector< vector< double>>> uz;
       vector< vector< vector< double>>> forcex;
       vector< vector< vector< double>>> forcey;
       vector< vector< vector< double>>> forcez;
       double getVelocityX(const Vector3d& position) const;
       double getVelocityY(const Vector3d& position) const;
       double getVelocityZ(const Vector3d& position) const;
       void makeBoundary();
       bool isInside(const Vector3d& position) const;
       void setForceX(double fx, const Vector3d& position);
       void setForceY(double fy, const Vector3d& position);
       void setForceZ(double fz, const Vector3d& position);
       void clearForce();
       void initVelocity();
       void initPressure();
       void initMarkers();
       Vector3d getLastPosition(const Vector3d& currentPosition, double dt);
       vector<int> allocator;
       void addGravityForce(double dt);
       void sortMarkers();
       bool existsMarker(int cellIndex_x, int cellIndex_y, int cellIndex_z);
       int index(int i, int j, int k) {return k*(Nx*Ny) + j*Ny + i;}
       vector<T> tripletList;
       vector<T> newTripletList;
};

#endif // ST_FIELD_H_INCLUDED
