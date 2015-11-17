#include <iostream>
#include "Cell.h"
#include "Grid.h"
#include "Vector2.h"

int main(int argc, char** argv){
    Grid field(512, 512);
    for(int i=0; i<512; i++){
        for(int j=0; j<512; j++){
            field.cells[i][j].u.x = 1.0;
            field.cells[i][j].u.y = 0.0;
        } 
    }
    std::cout << field.getVelocity(36.6, 105.4).y << std::endl;
    return 0;
}
