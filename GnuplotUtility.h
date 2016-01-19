#ifndef ST_GNUPLOTUTILITY_H_INCLUDED
#define ST_GNUPLOTUTILITY_H_INCLUDED
#include <iostream>
#include "Field.h"
using namespace std;

class GnuplotUtility {
    public:
        static void init_gnuplot(const Field& field) {
            cout << "set xrange [" << 0 << ":" << field.width * field.cellSize <<"]" << endl;
            cout << "set yrange [" << 0 << ":" << field.height * field.cellSize<<"]" << endl;
        }   
        static void init_gif(const Field& field) {
            cout << "set terminal gif animate optimize size 1024, 1024" << endl;
            cout << "set output 'fluid.gif'" << endl;
        }
        static void export_u0field_to_gnuplot(const Field& field) {
            int id = 0;
            for(int i=0; i < field.width; i++) {
                for(int j=0; j < field.height; j++) {
                    id++;   
                    cout << "set arrow " 
                         << id << " from " 
                         << (i + 0.5) * field.cellSize << "," 
                         << (j + 0.5) * field.cellSize << " to " 
                         << (i + 0.5) * field.cellSize + field.cells[i][j].u0.x() << "," 
                         << (j + 0.5) * field.cellSize + field.cells[i][j].u0.y() << endl;
                } 
            }
            cout << "plot 0/1 notitle" << endl;
        }
};

#endif // ST_GNUPLOTUTILITY_H_INCLUDED
