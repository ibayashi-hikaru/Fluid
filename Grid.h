#ifndef ST_GRID_H_INCLUDED
#define ST_GRID_H_INCLUDED
#include "Cell.h"
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
   int length, height;
   double cellSize;
   vector< vector<Cell> > cells;
   Grid(int l, int h){
        CALC_STEP = STEP0;
        this->length = l;
        this->height = h;
        this->cells = vector<vector<Cell>>(h, vector<Cell>(l, Cell())); 
   }
   double divergence(int x, int y) const;
   void addForce(double dt);
protected:
};

#endif // ST_GRID_H_INCLUDED

