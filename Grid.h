#ifndef ST_GRID_H_INCLUDED
#define ST_GRID_H_INCLUDED
#include <Eigen/Core>
#include <Eigen/unsupported/Eigen/FFT>
#include "Cell.h"
#include <vector>
using namespace std;
using namespace Eigen;
/*
    Class's description
*/
enum STEP{
    STEP0 = 0,
    STEP1 = 1,
    STEP2 = 2,
    STEP3 = 3,
    STEP4 = 4,
};

class Grid 
{
public:
   STEP CALC_STEP;
   unsigned long width, height;
   double cellSize;
   vector< vector<Cell> > cells;
   vector< vector<complex<double>>> ft_vx;
   vector< vector<complex<double>>> ft_vy;
   Grid(unsigned long w, unsigned long h){
        CALC_STEP = STEP0;
        this->width = w;
        this->height = h;
        this->cellSize = 1.0;
        this->cells = vector<vector<Cell>>{h, vector<Cell>{w, Cell{}}};
        for(int i = 0; i < h; i++){
            for(int j = 0; j < w; j++){
                this->cells.at(i).at(j).position = Vector2d(i + 0.5, j + 0.5);
            }
        }
        this->ft_vx = vector< vector<complex<double>>>(w, vector< complex<double>>(h, complex<double>(0.0, 0.0)));
        this->ft_vy = vector< vector<complex<double>>>(w, vector< complex<double>>(h, complex<double>(0.0, 0.0)));
   }
   Vector2d getVelocity(Vector2d position) const; 
   void FFT2d();
   void invFFT2d();
   Vector2d traceParticle(Vector2d position, double dt) const;
   void addForce(double dt);
   void addTransport(double dt);
   void addDiffuse(double dt);
   void projectField();
   void swapVelocity();
protected:
};

#endif // ST_GRID_H_INCLUDED

