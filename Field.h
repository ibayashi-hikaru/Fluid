#ifndef ST_FIELD_H_INCLUDED
#define ST_FIELD_H_INCLUDED
#include <Eigen/Core>
#include <eigen3/unsupported/Eigen/FFT>
#include "Cell.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "Main.h"
#include "FieldUtility.h"
using namespace std;
using namespace Eigen;
/*
    Class's description
*/

class Field {
    public:
       unsigned long width, height;
       double cellSize;
       vector< vector<Cell> > cells;
       vector< vector<complex<double>>> ft_vx;
       vector< vector<complex<double>>> ft_vy;
       vector< vector<double>> div;
       Field(unsigned long w, unsigned long h) {
            this->width = w;
            this->height = h;
            this->cellSize = 1.0;
            this->cells = vector<vector<Cell>>{w, vector<Cell>{h, Cell{}}};
            for(int i = 0; i < w; i++) {
                for(int j = 0; j < h; j++) {
                    this->cells.at(i).at(j).position = Vector2d((i + 0.5) * cellSize, (j + 0.5) * cellSize);
                }
            }
            this->ft_vx = vector< vector<complex<double>>>(w, vector< complex<double>>(h, complex<double>(0.0, 0.0)));
            this->ft_vy = vector< vector<complex<double>>>(w, vector< complex<double>>(h, complex<double>(0.0, 0.0)));
            this->div = vector<vector<double>>(w, vector<double>(h));
       }
       void initVelocity();
       Vector2d getVelocity(Vector2d position) const; 
       void FFT2d();
       void invFFT2d();
       Vector2d traceParticle(Vector2d position, double dt) const;
       void addForce(double dt);
       void addTransport(double dt);
       void addDiffuse(double dt);
       void projectField();
       void swapVelocity();
       void makeSquareForceSource();
       void makeLineForceSource();
       void makeDualForceSource();
       void resetForceSource();
       void updateRot();
    private:
       Vector2d periodizedPosition(Vector2d position) const;
       Vector2i getNearestPointIndices(Vector2d position) const;
       Vector2i getIndicesOfDiscretePosition(Vector2d discretePosition) const;
       vector<Vector2d> getSurroundingVelocities(Vector2i index) const;
       bool isEdge(Vector2i positionIndices) const;
       bool isCorner(Vector2i positionIndices) const;
       bool isRightOrLeftSide(Vector2i positionIndices) const;
       bool isUpOrDownSide(Vector2i positionIndices) const;
};

#endif // ST_FIELD_H_INCLUDED
