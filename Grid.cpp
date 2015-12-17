#include <iostream>
#include <cmath>
#include "Grid.h"
#include "Main.h"

Vector2d
Grid::getVelocity(Vector2d position) const{
}

void
Grid::FFT2d(){
    vector< vector<double>> in_rows_x(height, vector<double>(width, 0.0));
    vector< vector<double>> in_rows_y(height, vector<double>(width, 0.0));
    for(int i=0; i < height; i++){
        vector<double> tmp_row_x(width), tmp_row_y(width);
        for(int j=0; j < width; j++){
           tmp_row_x.at(j) = cells[i][j].u1.x(); 
           tmp_row_y.at(j) = cells[i][j].u1.y(); 
        }
        in_rows_x.at(i) = tmp_row_x;
        in_rows_y.at(i) = tmp_row_y;
    }
    FFT<double> fft;
    vector< vector< complex<double>>> med_x(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_y(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    for (int i = 0; i < height; i++) {
        fft.fwd(med_x.at(i), in_rows_x.at(i));
        fft.fwd(med_y.at(i), in_rows_y.at(i));
    }
    vector< vector< complex<double>>> med_cols_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_cols_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for(int i=0; i < width; i++){
        for(int j=0; j < height; j++){
            med_cols_x.at(i).at(j) = med_x.at(j).at(i); 
            med_cols_y.at(i).at(j) = med_y.at(j).at(i); 
        }
    }
    vector< vector< complex<double>>> out_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> out_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for (int i = 0; i < width; i++) {
        fft.fwd(out_x.at(i), med_cols_x.at(i));
        fft.fwd(out_y.at(i), med_cols_y.at(i));
    }
    for(int i=0; i < height; i++){
        for(int j=0; j < width; j++){
            ft_vx.at(i).at(j) = out_x.at(j).at(i); 
            ft_vy.at(i).at(j) = out_y.at(j).at(i);
        }
    }
}
void
Grid::invFFT2d(){
    FFT<double> fft;
    vector< vector< complex<double>>> med_x(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_y(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    for(int i=0; i < height; i++){
        fft.inv(med_x.at(i), ft_vx.at(i));
        fft.inv(med_y.at(i), ft_vy.at(i));
    }
    vector< vector< complex<double>>> med_cols_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_cols_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for(int i=0; i < width; i++){
        for(int j=0; j < height; j++){
            med_cols_x.at(i).at(j) = med_x.at(j).at(i); 
            med_cols_y.at(i).at(j) = med_y.at(j).at(i); 
        }
    }
    vector< vector< complex<double>>> out_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> out_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for (int i = 0; i < width; i++) {
        fft.inv(out_x.at(i), med_cols_x.at(i));
        fft.inv(out_y.at(i), med_cols_y.at(i));
    }
    for(int i=0; i < height; i++){
        for(int j=0; j < width; j++){
            cells.at(i).at(j).u1.x() = out_x.at(j).at(i).real(); 
            cells.at(i).at(j).u1.y() = out_y.at(j).at(i).real();
        }
    }

}
// Runge-Kutta
Vector2d
Grid::traceParticle(Vector2d position, double dt) const{
    Vector2d k0 = dt * getVelocity(position);
    Vector2d k1 = dt * getVelocity(position - k0/2.0); 
    Vector2d k2 = dt * getVelocity(position - k1/2.0); 
    Vector2d k3 = dt * getVelocity(position - k2);
    return position - (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0; 
}

void
Grid::addForce(double dt){
    if(CALC_STEP == STEP0){
        for(int i=0; i<width; i++){
            for(int j=0; j<height; j++){
                cells[i][j].u1.x() = cells[i][j].u0.x() + dt*cells[i][j].force.x(); 
                cells[i][j].u1.y() = cells[i][j].u0.y() + dt*cells[i][j].force.y(); 
            } 
        }
        CALC_STEP = STEP1;
    } else {
        std::cout << "Force cannot be applied at this step." << std::endl; 
    }    
}

void
Grid::addTransport(double dt){
    if(CALC_STEP == STEP1){
        for(int i=0; i<width; i++){
            for(int j=0; j<height; j++){
                Vector2d current_position{i + 0.5, j + 0.5};
                Vector2d last_position = Grid::traceParticle(current_position, dt);
                cells[i][j].u1 += Grid::getVelocity(last_position) - cells[i][j].u0;
            } 
        }
        CALC_STEP = STEP2;
    } else {
        std::cout << "Transport cannot be applied at this step." << std::endl; 
    }
}

void
Grid::addDiffuse(double dt){
    if(CALC_STEP == STEP2){
        for(int i=0; i<width; i++){
            for(int j=0; j<height; j++){
                complex<double> ikx = complex<double>(0.0, (2.0*PI*i)/width); 
                complex<double> iky = complex<double>(0.0, (2.0*PI*j)/height); 
                ft_vx.at(i).at(j) =  ft_vx.at(i).at(j)/(1.0 - NU * dt * (ikx * ikx + iky * iky));
                ft_vy.at(i).at(j) =  ft_vy.at(i).at(j)/(1.0 - NU * dt * (ikx * ikx + iky * iky));
            } 
        }
        CALC_STEP = STEP3;
    } else {
        std::cout << "Diffuse cannot be applied at this step." << std::endl; 
    }
}

void
Grid::projectField(){
    if(CALC_STEP == STEP3){
        double inv_l = (double) 1.0/width;   
        double inv_h = (double) 1.0/height;   
        for(int i=0; i<width; i++){
            for(int j=0; j<height; j++){
                if(i == 0 && j == 0){
                    ft_vx.at(i).at(j) -= 0.0;
                    ft_vy.at(i).at(j) -= 0.0;
                } else {
                    complex<double> ikx = complex<double>(0.0, 2.0*PI*i*inv_l); 
                    complex<double> iky = complex<double>(0.0, 2.0*PI*j*inv_h);
                    double ik2 = -((2.0*PI*i*inv_l) * (2.0*PI*i*inv_l) + (2.0*PI*j*inv_h) * (2.0*PI*j*inv_h));
                    complex<double> ik_dot_w = ikx * ft_vx.at(i).at(j) + iky * ft_vy.at(i).at(j); // This variable name is based on the paper "stable fluid"
                    ft_vx.at(i).at(j) -= (1.0/ik2) * ik_dot_w * ikx;
                    ft_vy.at(i).at(j) -= (1.0/ik2) * ik_dot_w * iky;
                }
            } 
        }
        CALC_STEP = STEP4;
    } else {
        std::cout << "Project cannot be applied at this step." << std::endl; 
    }
}

void
Grid::swapVelocity(){
    for(int i=0; i<width; i++){
        for(int j=0; j<height; j++){
            cells[i][j].u0 = cells[i][j].u1;
        } 
    }
}
