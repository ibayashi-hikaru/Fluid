#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include <string>
#include "Main.h"
#include "Cell.h"
#include "Field.h"
#include "GnuplotUtility.h"

const double PI = M_PI;
const double NU = 0.01;
using namespace Eigen;
void fieldInit(Field& field);
void fieldInit(Field& field){
    field.initVelocity(); 
    field.makeLineForceSource();
}
int main(int argc, char** argv) {
    bool gif_flag = false;
    bool plot_flag = false;
    int time_cnt = 50;
    for(int i = 0; i < argc; i++) {
        if(strncmp(argv[i], "-gif", 4) == 0) {
            gif_flag = true; 
        }
        if(strncmp(argv[i], "-plot", 5) == 0) {
            plot_flag = true; 
        } 
        if(strncmp(argv[i], "-t", 2) == 0) {
            std::string a{argv[i]};   
            time_cnt = std::stoi(a.substr(2, -1));
        } 
    } 
    if(gif_flag || plot_flag) GnuplotUtility::init_gnuplot(field);
    if(gif_flag) GnuplotUtility::init_gif(field);
    for(int i = 0; i < time_cnt; i++) {
        field.addForce(0.1);
        field.addTransport(0.1);
        field.FFT2d();
        field.addDiffuse(0.1);
        field.projectField();
        field.invFFT2d();
        field.swapVelocity();
        if(gif_flag) GnuplotUtility::export_u0field_to_gnuplot(field); 
    Field field = Field(32, 32);
    fieldInit(field);
    }
    if(plot_flag) GnuplotUtility::export_u0field_to_gnuplot(field);
    return 0;
}
