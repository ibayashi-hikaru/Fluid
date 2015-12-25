#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include "Main.h"
#include "Cell.h"
#include "Field.h"
#include "InterfaceUtility.h"

const double PI = M_PI;
const double NU = 0.01;
using namespace Eigen;
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
    InterfaceUtility::export_u0field_to_gnuplot(field);
    return 0;
}
