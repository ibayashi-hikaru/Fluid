#ifndef ST_FIELD_H_INCLUDED
#define ST_FIELD_H_INCLUDED
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Eigen;
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
       }
       int GridNum() const {return Nx;};
       double Dx() const {return dx;};
       void Init();
       void Advect(double dt);
       void AddForce(double dt);      
       void CG_Project(double dt);
       void SetForce(Vector3d force, Vector3d position);
       Vector3d GetVelocity(Vector3d position) const;
    private:
       const double rho = 1.0;
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
       double getVelocityX(Vector3d position) const;
       double getVelocityY(Vector3d position) const;
       double getVelocityZ(Vector3d position) const;
       void makeBoundary();
       bool isInside(Vector3d position) const;
       void setForceX(double fx, Vector3d position);
       void setForceY(double fy, Vector3d position);
       void setForceZ(double fz, Vector3d position);
       void clearForce();
       void initVelocity();
       void initPressure();
       Vector3d getLastPosition(Vector3d currentPosition, double dt);
};

#endif // ST_FIELD_H_INCLUDED
