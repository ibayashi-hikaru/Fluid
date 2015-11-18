#include <iostream>
#include <cmath>
#include "Vector2.h"
#include "Grid.h"

Vector2<double>
Grid::getVelocity(const Vector2<double>& position) const{
    if(position.x < 0.0 || length * cellSize < position.x || position.y < 0.0 || height * cellSize < position.y ){
        std::cout << "Position is out of Grid." << std::endl;
        return Vector2<double>{}; 
    } else {
        double dummy_x = position.x - cellSize/2.0;
        double dummy_y = position.y - cellSize/2.0;
        // Decide 4 points' velocity near to target point
        // f: floor, c: ceil, i: interporated
        Vector2<double> ffu{}, fcu{}, cfu{}, ccu{};
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
        Vector2<double> ifu{}, icu{}; 
        ifu = ffu * (ceil(dummy_x) - dummy_x) + cfu * (dummy_x - floor(dummy_x));
        icu = fcu * (ceil(dummy_x) - dummy_x) + ccu * (dummy_x - floor(dummy_x));
        // Then y direction and return
        Vector2<double> iiu{};
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
        div_x = (cells[0][y].u0.x - cells[length-1][y].u0.x)/cellSize;
    } else {
        div_x = (cells[x][y].u0.x - cells[x-1][y].u0.x)/cellSize;
    }
    if(y == 0){
        div_y = (cells[x][0].u0.y - cells[x][height-1].u0.y)/cellSize;
    } else {
        div_y = (cells[x][y].u0.y - cells[x][y-1].u0.y)/cellSize;
    }
    return div_x + div_y;
}

// Runge-Kutta
Vector2<double>
Grid::traceParticle(const Vector2<double>& position, double dt) const{
   Vector2<double> k0 = dt * getVelocity(position);
   Vector2<double> k1 = dt * getVelocity(position + k0/2.0); 
   Vector2<double> k2 = dt * getVelocity(position + k1/2.0); 
   Vector2<double> k3 = dt * getVelocity(position + k2);
   return (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0; 
}

void
Grid::addForce(double dt){
    if(CALC_STEP == STEP0){
        for(int i=0; i<length; i++){
            for(int j=0; j<height; j++){
                cells[i][j].u1.x = cells[i][j].u0.x + dt*cells[i][j].force.x; 
                cells[i][j].u1.y = cells[i][j].u0.y + dt*cells[i][j].force.y; 
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
                Vector2<double> currentPosition{i + 0.5, j + 0.5};
                cells[i][j].u1 += Grid::traceParticle(currentPosition, dt) - cells[i][j].u0;
            } 
        }
        CALC_STEP = STEP2;
    } else {
        std::cout << "Transport cannot be applied at this step." << std::endl; 
    }
}
