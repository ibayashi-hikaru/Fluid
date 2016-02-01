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
            this->width = w;
            this->height = h;
            this->cellSize = 1.0;
            this->div = vector<vector<double>>(w, vector<double>(h));
            this->p = vector<vector<double>>(w, vector<double>(h));
            this->ux0 = vector<vector<double>>(w + 1, vector<double>(h));
            this->ux1 = vector<vector<double>>(w + 1, vector<double>(h));
            this->uy0 = vector<vector<double>>(w, vector<double>(h + 1));
            this->uy1 = vector<vector<double>>(w, vector<double>(h + 1));
       }
       void AddForce(double dt);      
    private:
       unsigned long width, height;
       double cellSize;
       vector< vector<double>> div;
       vector< vector<double>> p;
       vector< vector<double>> ux0, ux1;
       vector< vector<double>> uy0, uy1;
};

#endif // ST_FIELD_H_INCLUDED
