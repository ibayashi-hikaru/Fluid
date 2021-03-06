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
            this->dx = 1.0;
            this->div = vector<vector<double>>(gridNum, vector<double>(gridNum));
            this->p = vector<vector<double>>(gridNum, vector<double>(gridNum));
            this->ux = vector<vector<double>>(gridNum + 1, vector<double>(gridNum));
            this->uy = vector<vector<double>>(gridNum, vector<double>(gridNum + 1));
            this->forcex = vector<vector<double>>(gridNum + 1, vector<double>(gridNum));
            this->forcey = vector<vector<double>>(gridNum, vector<double>(gridNum + 1));
       }
       int GridNum() const {return Nx;};
       double Dx() const {return dx;};
       void Init();
       void Advect(double dt);
       void AddForce(double dt);      
       void GS_Project(double dt);
       void CG_Project(double dt);
       void SetForce(Vector2d force, Vector2d position);
       Vector2d GetVelocity(Vector2d position) const;
       Vector2d TransformDisplayToField(Vector2d displayPosition, int width, int height) const;
       Vector2d TransformFieldToDisplay(Vector2d fieldPosition, int width, int height) const;
    private:
       const double rho = 1.0;
       unsigned long Nx, Ny;
       double dx;
       vector< vector<double>> div;
       vector< vector<double>> p;
       vector< vector<double>> ux;
       vector< vector<double>> uy;
       vector< vector<double>> forcex;
       vector< vector<double>> forcey;
       double getVelocityX(double x, double y) const;
       double getVelocityY(double x, double y) const;
       void makeBoundary();
       bool isInside(Vector2d position) const;
       void setForceX(double fx, Vector2d position);
       void setForceY(double fy, Vector2d position);
       void clearForce();
       void initVelocity();
       void initPressure();
       Vector2d getLastPosition(Vector2d currentPosition, double dt);
};

#endif // ST_FIELD_H_INCLUDED
