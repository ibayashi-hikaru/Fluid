#include <iostream>
#include <cmath>
#include "Grid.h"
#include "Main.h"

Vector2d
Grid::getVelocity(Vector2d position) const{
    if(position.x() < 0.0 || length * cellSize < position.x() || position.y() < 0.0 || height * cellSize < position.y() ){
        std::cout << "Position is out of Grid." << std::endl;
        return Vector2d::Zero(); 
    } else {
        double dummy_x = position.x() - cellSize/2.0;
        double dummy_y = position.y() - cellSize/2.0;
        // Decide 4 points' velocity near to target point
        // f: floor, c: ceil, i: interporated
        Vector2d ffu, fcu, cfu, ccu;
        if(floor(dummy_x) == -1){
            ffu = cells[length][floor(dummy_y)].u0;
            fcu = cells[length][ceil(dummy_y)].u0;
            cfu = cells[ceil(dummy_x)][floor(dummy_y)].u0;
            ccu = cells[ceil(dummy_x)][ceil(dummy_y)].u0;
        } else if(floor(dummy_y) == -1){
            ffu = cells[floor(dummy_x)][height].u0;
            fcu = cells[floor(dummy_x)][ceil(dummy_y)].u0;
            cfu = cells[ceil(dummy_x)][height].u0;
            ccu = cells[ceil(dummy_x)][ceil(dummy_y)].u0;
        } else {
            ffu = cells[floor(dummy_x)][floor(dummy_y)].u0;
            fcu = cells[floor(dummy_x)][ceil(dummy_y)].u0;
            cfu = cells[ceil(dummy_x)][floor(dummy_y)].u0;
            ccu = cells[ceil(dummy_x)][ceil(dummy_y)].u0;
        }
        // Interporate velocity field.
        // First x direction.
        Vector2d ifu, icu; 
        ifu = ffu * (ceil(dummy_x) - dummy_x) + cfu * (dummy_x - floor(dummy_x));
        icu = fcu * (ceil(dummy_x) - dummy_x) + ccu * (dummy_x - floor(dummy_x));
        // Then y direction and return
        Vector2d iiu;
        iiu = ifu * (ceil(dummy_y) - dummy_y) + icu * (dummy_y - floor(dummy_y));
        return iiu; 
    }
}

double
Grid::divergence(int x, int y) const{
    if(x < 0 || x >= width || y < 0 || y >= height ){
        std::cout << "(" <<  x << ","  << y << ") is out of field" << std::endl;
        return 0.0; 
    }
    // Use backward difference
    double div_x, div_y;
    if(x == 0){
        div_x = (cells[0][y].u0.x() - cells[width-1][y].u0.x())/cellSize;
    } else {
        div_x = (cells[x][y].u0.x() - cells[x-1][y].u0.x())/cellSize;
    }
    if(y == 0){
        div_y = (cells[x][0].u0.y() - cells[x][height-1].u0.y())/cellSize;
    } else {
        div_y = (cells[x][y].u0.y() - cells[x][y-1].u0.y())/cellSize;
    }
    return div_x + div_y;
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
