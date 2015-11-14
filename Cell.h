#ifndef ST_CELL_H_INCLUDED
#define ST_CELL_H_INCLUDED 

/*
    Class's description
*/
class Cell
{
public:
    double u0[2], u1[2];
    double density;
    double length;
    Cell(){
        this->u0[0] = 0.0;
        this->u0[1] = 0.0;
        this->density = 1,0;
        this->length = 1,0;
    }
    Cell(double x, double y, double d, double l){
        this->u0[0] = x;
        this->u0[1] = y;
        this->density = d;
        this->length = l;
    }
protected:
};

#endif // ST_CELL_H_INCLUDED

