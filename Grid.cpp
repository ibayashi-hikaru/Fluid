#include <iostream>
#include <cmath>
#include "Grid.h"

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
    if(x < 0 || x >= length || y < 0 || y >= height ){
        std::cout << "(" <<  x << ","  << y << ") is out of field" << std::endl;
        return 0.0; 
    }
    // Use backward difference
    double div_x, div_y;
    if(x == 0){
        div_x = (cells[0][y].u0.x() - cells[length-1][y].u0.x())/cellSize;
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
    vector< vector<double>> in_rows_x(height, vector<double>(length, 0.0));
    vector< vector<double>> in_rows_y(height, vector<double>(length, 0.0));
    for(int i=0; i < height; i++){
        vector<double> tmp_row_x(length), tmp_row_y(length);
        for(int j=0; j < length; j++){
           tmp_row_x.at(i) = cells[i][j].u1.x(); 
           tmp_row_y.at(i) = cells[i][j].u1.y(); 
        }
        in_rows_x.at(i) = tmp_row_x;
        in_rows_y.at(i) = tmp_row_y;
    }
    FFT<double> fft;
    vector< vector< complex<double>>> med_x(height, vector< complex<double>>(length, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_y(height, vector< complex<double>>(length, complex<double>(0.0, 0.0)));
    for (int i = 0; i < height; i++) {
        vector< complex<double>> tmp_med_x(length);
        vector< complex<double>> tmp_med_y(length);
        fft.fwd(tmp_med_x, in_rows_x.at(i));
        fft.fwd(tmp_med_y, in_rows_y.at(i));
        med_x.at(i) = tmp_med_x;
        med_y.at(i) = tmp_med_y;
    }
    vector< vector< complex<double>>> med_cols_x(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_cols_y(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for(int i=0; i < length; i++){
        for(int j=0; j < height; j++){
            med_cols_x.at(i).at(j) = med_x.at(j).at(i); 
            med_cols_y.at(i).at(j) = med_y.at(j).at(i); 
        }
    }
    vector< vector< complex<double>>> out_x(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> out_y(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for (int i = 0; i < length; i++) {
        vector< complex<double>> tmp_out_x(height);
        vector< complex<double>> tmp_out_y(height);
        fft.fwd(tmp_out_x, med_cols_x.at(i));
        fft.fwd(tmp_out_y, med_cols_y.at(i));
        out_x.at(i) = tmp_out_x;
        out_y.at(i) = tmp_out_y;
    }
    for(int i=0; i < height; i++){
        for(int j=0; j < length; j++){
            ft_vx.at(i).at(j) = out_x.at(j).at(i); 
            ft_vy.at(i).at(j) = out_y.at(j).at(i);
        }
    }
}
void
Grid::invFFT2d(){
    FFT<double> fft;
    vector< vector< complex<double>>> med_x(height, vector< complex<double>>(length, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_y(height, vector< complex<double>>(length, complex<double>(0.0, 0.0)));
    for(int i=0; i < height; i++){
        vector< complex<double>> tmp_med_x(length);
        vector< complex<double>> tmp_med_y(length);
        fft.inv(tmp_med_x, ft_vx.at(i));
        fft.inv(tmp_med_y, ft_vy.at(i));
        med_x.at(i) = tmp_med_x;
        med_y.at(i) = tmp_med_y;
    }
    vector< vector< complex<double>>> med_cols_x(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_cols_y(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for(int i=0; i < length; i++){
        for(int j=0; j < height; j++){
            med_cols_x.at(i).at(j) = med_x.at(j).at(i); 
            med_cols_y.at(i).at(j) = med_y.at(j).at(i); 
        }
    }
    vector< vector< complex<double>>> out_x(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> out_y(length, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for (int i = 0; i < length; i++) {
        vector< complex<double>> tmp_out_x(height);
        vector< complex<double>> tmp_out_y(height);
        fft.inv(tmp_out_x, med_cols_x.at(i));
        fft.inv(tmp_out_y, med_cols_y.at(i));
        out_x.at(i) = tmp_out_x;
        out_y.at(i) = tmp_out_y;
    }
    for(int i=0; i < height; i++){
        for(int j=0; j < length; j++){
            cells.at(i).at(j).u1.x() = out_x.at(j).at(i).real(); 
            cells.at(i).at(j).u1.y() = out_y.at(j).at(i).real();
        }
    }

}
// Runge-Kutta
Vector2d
Grid::traceParticle(Vector2d position, double dt) const{
    Vector2d k0 = dt * getVelocity(position);
    Vector2d k1 = dt * getVelocity(position + k0/2.0); 
    Vector2d k2 = dt * getVelocity(position + k1/2.0); 
    Vector2d k3 = dt * getVelocity(position + k2);
    return (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0; 
}

void
Grid::addForce(double dt){
    if(CALC_STEP == STEP0){
        for(int i=0; i<length; i++){
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
        for(int i=0; i<length; i++){
            for(int j=0; j<height; j++){
                Vector2d currentPosition{i + 0.5, j + 0.5};
                cells[i][j].u1 += Grid::traceParticle(currentPosition, dt) - cells[i][j].u0;
            } 
        }
        CALC_STEP = STEP2;
    } else {
        std::cout << "Transport cannot be applied at this step." << std::endl; 
    }
}
void
Grid::swapVelocity(){
    for(int i=0; i<length; i++){
        for(int j=0; j<height; j++){
            cells[i][j].u0 = cells[i][j].u1;
        } 
    }
}
