#include <iostream>
#include "Cell.h"
#include "Grid.h"

int main(int argc, char** argv){
    Grid field(512, 512);
    for(int i=0; i<512; i++){
        for(int j=0; j<512; j++){
            field.cells[i][j].u0[0] = 1.0;
            field.cells[i][j].u0[1] = 0.0;
        } 
    }
    std::cout << field.cells[100][100].u0[0] << std::endl;
    std::cout << "hello fluid!" << std::endl;
    return 0;
}
