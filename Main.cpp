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
    bool gif_flag = false;
    bool plot_flag = false;
    for(int i = 0; i < argc; i++) {
        if(strncmp(argv[i], "-gif", 4) == 0) {
            gif_flag = true; 
        }
        if(strncmp(argv[i], "-plot", 5) == 0) {
            plot_flag = true; 
        } 
    } 
    Field field(64, 64);
    for(int i=0; i<field.height; i++) {
        for(int j=0; j<field.width; j++) {
            field.cells[i][j].u0.x() = 1.0;
            field.cells[i][j].u0.y() = 0.0;
            field.cells[i][j].u1.x() = 1.0;
            field.cells[i][j].u1.y() = 0.0;
        } 
    }
    field.makeLineForceSource();
    field.CALC_STEP = STEP0;
    if(gif_flag || plot_flag) InterfaceUtility::init_gnuplot(field);
    if(gif_flag) InterfaceUtility::init_gif(field);
    for(int i = 0; i < 50; i++) {
        field.addForce(0.1);
        field.addTransport(0.1);
        field.FFT2d();
        field.addDiffuse(0.1);
        field.projectField();
        field.invFFT2d();
        field.swapVelocity();
        if(gif_flag) InterfaceUtility::export_u0field_to_gnuplot(field); 
    }
    if(plot_flag) InterfaceUtility::export_u0field_to_gnuplot(field);
    return 0;
}
