#include <iostream>
#include <Eigen/Core>
#include "Cell.h"
#include "Grid.h"

using namespace Eigen;
int main(int argc, char** argv){
    Grid field(512, 512);
    for(int i=0; i<512; i++){
        for(int j=0; j<512; j++){
            field.cells[i][j].u0.x() = 1.0;
            field.cells[i][j].u0.y() = 0.0;
            field.cells[i][j].u1.x() = 1.0;
            field.cells[i][j].u1.y() = 0.0;
        } 
    }
    std::cout << field.getVelocity(Vector2d{36.6, 105.4}) << std::endl;
    return 0;
}
