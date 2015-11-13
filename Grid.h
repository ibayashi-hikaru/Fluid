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
   vector< vector<Cell> > cells;
protected:
};

#endif // ST_GRID_H_INCLUDED

