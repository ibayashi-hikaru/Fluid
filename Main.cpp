#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include "Main.h"
#include "Cell.h"
#include "Grid.h"

const double PI = 3.14159265358979323846;
const double NU = 0.01;
using namespace Eigen;
void export_velocity_field(Grid field);
int main(int argc, char** argv){
    Grid field(2, 2);
    for(int i=0; i<field.height; i++){
        for(int j=0; j<field.width; j++){
            field.cells[i][j].u0.x() = 1.0;
            field.cells[i][j].u0.y() = 0.0;
            field.cells[i][j].u1.x() = 1.0;
            field.cells[i][j].u1.y() = 0.0;
        } 
    }
    for(int i=0; i < field.width; i++){
        field.cells[i][field.height/2].force.y() = 4.0;
    }
    field.CALC_STEP = STEP0;
    field.addForce(0.1);
    field.addTransport(0.1);
    field.FFT2d();
    field.addDiffuse(0.1);
    field.projectField();
    field.invFFT2d();
    field.swapVelocity();
    export_velocity_field(field);
    return 0;
}

void export_velocity_field(Grid field){
    std::cout << "set xrange [" << 0 << ":" << field.width <<"]" << std::endl;
    std::cout << "set yrange [" << 0 << ":" << field.height <<"]" << std::endl;
    int id = 0;
    for(int i=0; i < field.width; i++){
        for(int j=0; j < field.height; j++){
            id++;   
            std::cout << "set arrow " << id << " from " << i << "," << j << " to " << i + field.cells[i][j].u0.x() << "," <<  j + field.cells[i][j].u0.y() << std::endl;
        } 
    }
    std::cout << "plot 0/1 notitle" << std::endl;
}
