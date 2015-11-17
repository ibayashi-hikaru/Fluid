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
    double size;
    Cell(){
        this->u0[0] = 0.0;
        this->u0[1] = 0.0;
        this->density = 1,0;
        this->size = 1,0;
    }
    Cell(double x, double y, double d, double s){
        this->u0[0] = x;
        this->u0[1] = y;
        this->density = d;
        this->size = s;
    }
protected:
};

#endif // ST_CELL_H_INCLUDED

