#ifndef ST_FIELD_H_INCLUDED
#define ST_FIELD_H_INCLUDED
#include <Eigen/Core>
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
       Field(unsigned long w, unsigned long h) {
            this->Nx = w;
            this->Ny = h;
            this->dx = 1.0;
            this->div = vector<vector<double>>(w, vector<double>(h));
            this->p = vector<vector<double>>(w, vector<double>(h));
            this->ux0 = vector<vector<double>>(w + 1, vector<double>(h));
            this->ux1 = vector<vector<double>>(w + 1, vector<double>(h));
            this->uy0 = vector<vector<double>>(w, vector<double>(h + 1));
            this->uy1 = vector<vector<double>>(w, vector<double>(h + 1));
            this->forcex = vector<vector<double>>(w + 1, vector<double>(h));
            this->forcey = vector<vector<double>>(w, vector<double>(h + 1));
       }
       void Init();
       void Advect(double dt);
       void AddForce(double dt);      
       void Project(double dt);
       void SetForce(Vector2d force, Vector2d position);
       Vector2d TransformDisplayToField(Vector2d displayPosition, int width, int height) const;
    private:
       const double rho = 1.0;
       unsigned long Nx, Ny;
       double dx;
       vector< vector<double>> div;
       vector< vector<double>> p;
       vector< vector<double>> ux0, ux1;
       vector< vector<double>> uy0, uy1;
       vector< vector<double>> forcex;
       vector< vector<double>> forcey;
       double getVelocityX(double x, double y) const;
       double getVelocityY(double x, double y) const;
       void makeBoundary();
       bool isInsideField(Vector2d position) const;
       void setForceX(double fx, Vector2d position);
       void setForceY(double fy, Vector2d position);
       void clearForce();
       void initVelocity();
};

#endif // ST_FIELD_H_INCLUDED
