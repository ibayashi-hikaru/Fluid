#ifndef ST_CELL_H_INCLUDED
#define ST_CELL_H_INCLUDED
#include <Eigen/Core>
using namespace Eigen;
/*
    Class's description
*/
class Cell {
    public:
        Vector2d u0, u1;
        Vector2d position;
        double density;
        double size;
        Vector2d force;
        Cell(){
            this->u0 = Vector2d::Zero();
            this->u1 = Vector2d::Zero();
            this->density = 1.0;
            this->size = 1.0;
            this->force = Vector2d::Zero();
        }
        Cell(double x, double y, double d, double s, double fx, double fy ){
            this->u0.x() = x;
            this->u0.y() = y;
            this->u1.x() = x;
            this->u1.y() = y;
            this->density = d;
            this->size = s;
            this->force.x() = fx;
            this->force.y() = fy;
        }
};
#endif // ST_CELL_H_INCLUDED
