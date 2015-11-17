#ifndef ST_GRID_H_INCLUDED
#define ST_GRID_H_INCLUDED
#include "Cell.h"
#include <vector>
using namespace std;
/*
    Class's description
*/
class Grid 
{
public:
   int length, height;
   double cellSize;
   vector< vector<Cell> > cells;
   Grid(int l, int h){
        this->length = l;
        this->height = h;
        this->cells = vector<vector<Cell>>(h, vector<Cell>(l, Cell())); 
   }
   double div(int x, int y) const;
protected:
};

#endif // ST_GRID_H_INCLUDED

