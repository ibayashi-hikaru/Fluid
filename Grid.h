#ifndef ST_GRID_H_INCLUDED
#define ST_GRID_H_INCLUDED
#include "Cell.h"
#include "Vector2.h"
#include <vector>
using namespace std;
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
   unsigned long length, height;
   double cellSize;
   vector< vector<Cell> > cells;
   Grid(unsigned long l, unsigned long h){
        CALC_STEP = STEP0;
        this->length = l;
        this->height = h;
        this->cellSize = 1.0;
        this->cells = vector<vector<Cell>>{h, vector<Cell>{l, Cell{}}}; 
   }
   Vector2<double> getVelocity(Vector2<double> position); 
   double divergence(int x, int y) const;
   Vector2<double> traceParticle(const Vector2<double>& position, double dt) const;
   void addForce(double dt);
protected:
};

#endif // ST_GRID_H_INCLUDED

