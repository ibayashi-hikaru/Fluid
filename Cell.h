#ifndef ST_CELL_H_INCLUDED
#define ST_CELL_H_INCLUDED
#include "Vector2.h"

/*
    Class's description
*/
class Cell
{
public:
    Vector2<double> u;
    double density;
    double size;
    Cell(){
        this->u.x = 0.0;
        this->u.y = 0.0;
        this->density = 1,0;
        this->size = 1,0;
    }
    Cell(double x, double y, double d, double s){
        this->u.x = x;
        this->u.y = y;
        this->density = d;
        this->size = s;
    }
protected:
};

#endif // ST_CELL_H_INCLUDED

