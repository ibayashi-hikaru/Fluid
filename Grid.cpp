#include <iostream>
#include "Grid.h"

double
Grid::div(int x, int y) const{
    if(x < 0 || x >= length || y < 0 || y >= height ){
        std::cout << "(" <<  x << ","  << y << ") is out of field" << std::endl;
        return 0.0; 
    }
    // Use backward difference
    double div_x, div_y;
    if(x == 0){
        div_x = (cells[1][y].u[0] - cells[0][y].u[0])/cellSize;
    } else {
        div_x = (cells[x][y].u[0] - cells[x-1][y].u[0])/cellSize;
    }
    if(y == 0){
        div_y = (cells[x][1].u[1] - cells[x][0].u[1])/cellSize;
    } else {
        div_y = (cells[x][y].u[1] - cells[x][y-1].u[1])/cellSize;
    }
    return div_x + div_y;
}
