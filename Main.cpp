#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include "Main.h"
#include "Cell.h"
#include "Field.h"

const double PI = 3.14159265358979323846;
const double NU = 0.01;
using namespace Eigen;
void export_velocity_field(Field field);
int main(int argc, char** argv) {
    Field field(2, 2);
    for(int i=0; i<field.height; i++) {
        for(int j=0; j<field.width; j++) {
            field.cells[i][j].u0.x() = 1.0;
            field.cells[i][j].u0.y() = 0.0;
            field.cells[i][j].u1.x() = 1.0;
            field.cells[i][j].u1.y() = 0.0;
        } 
    }
    for(int i=0; i < field.width; i++) {
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

void export_velocity_field(Field field) {
    std::cout << "set xrange [" << 0 << ":" << field.width * field.cellSize <<"]" << std::endl;
    std::cout << "set yrange [" << 0 << ":" << field.height * field.cellSize<<"]" << std::endl;
    int id = 0;
    for(int i=0; i < field.width; i++) {
        for(int j=0; j < field.height; j++) {
            id++;   
            std::cout << "set arrow " 
                      << id << " from " 
                      << (i + 0.5) * field.cellSize << "," 
                      << (j + 0.5) * field.cellSize << " to " 
                      << (i + 0.5) * field.cellSize + field.cells[i][j].force.x() << "," 
                      << (j + 0.5) * field.cellSize + field.cells[i][j].force.y() << std::endl;
        } 
    }
    std::cout << "plot 0/1 notitle" << std::endl;
}
