#include <iostream>
#include "Grid.h"

double
Grid::divergence(int x, int y) const{
    if(x < 0 || x >= length || y < 0 || y >= height ){
        std::cout << "(" <<  x << ","  << y << ") is out of field" << std::endl;
        return 0.0; 
    }
    // Use backward difference
    double div_x, div_y;
    if(x == 0){
        div_x = (cells[0][y].u.x - cells[length-1][y].u.x)/cellSize;
    } else {
        div_x = (cells[x][y].u.x - cells[x-1][y].u.x)/cellSize;
    }
    if(y == 0){
        div_y = (cells[x][0].u.y - cells[x][height-1].u.y)/cellSize;
    } else {
        div_y = (cells[x][y].u.y - cells[x][y-1].u.y)/cellSize;
    }
    return div_x + div_y;
}

void
Grid::addForce(double dt){
    for(int i=0; i<length; i++){
        for(int j=0; j<height; j++){
            cells[i][j].u.x += dt*cells[i][j].force.x; 
            cells[i][j].u.y += dt*cells[i][j].force.y; 
        } 
    }    
}
