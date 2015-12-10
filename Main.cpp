#include <iostream>
#include <Eigen/Core>
#include "Cell.h"
#include "Grid.h"

using namespace Eigen;
void export_velocity_field(Grid field);
int main(int argc, char** argv){
    Grid field(64, 64);
    for(int i=0; i<field.height; i++){
        for(int j=0; j<field.length; j++){
            field.cells[i][j].u0.x() = 1.0;
            field.cells[i][j].u0.y() = 0.0;
            field.cells[i][j].u1.x() = 1.0;
            field.cells[i][j].u1.y() = 0.0;
        } 
    }
    export_velocity_field(field);
    return 0;
}

void export_velocity_field(Grid field){
    std::cout << "set xrange [" << 0 << ":" << field.length <<"]" << std::endl;
    std::cout << "set yrange [" << 0 << ":" << field.height <<"]" << std::endl;
    int id = 0;
    for(int i=0; i < field.length; i++){
        for(int j=0; j < field.height; j++){
            id++;   
            std::cout << "set arrow " << id << " from " << i << "," << j << " to " << i + field.cells[i][j].u0.x() << "," <<  j + field.cells[i][j].u0.y() << std::endl;
        } 
    }
    std::cout << "plot 0/1 notitle" << std::endl;
}
