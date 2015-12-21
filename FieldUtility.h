#ifndef ST_FIELDUTILITY_H_INCLUDED
#define ST_FIELDUTILITY_H_INCLUDED
#include <Eigen/Core>
using namespace std;
using namespace Eigen;

class FieldUtility{
public:
    static Vector2d interpolate2d(Vector2d ff, Vector2d fc, Vector2d cf, Vector2d cc, Vector2d local_normalized_position){
        Vector2d fi = ff * (1.0 - local_normalized_position.y()) 
                    + fc * (local_normalized_position.y());
        Vector2d ci = cf * (1.0 - local_normalized_position.y()) 
                    + cc * (local_normalized_position.y());
        Vector2d ii = fi * (1.0 - local_normalized_position.x())
                    + ci * (local_normalized_position.x());
        return ii; 
    }
};
#endif // ST_FIELDUTILITY_H_INCLUDED
