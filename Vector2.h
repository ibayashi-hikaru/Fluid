#ifndef ST_VECTOR2_H_INCLUDED
#define ST_VECTOR2_H_INCLUDED

/*
    Class's description
*/
template <typename TYPE> class Vector2 
{
public:
     Vector2() : x(0), y(0){}
     Vector2(TYPE val_x, TYPE val_y) : x(val_x), y(val_y){}
     TYPE x, y;
     Vector2 operator+(const Vector2<TYPE>& v) const{
         return Vector2<TYPE> (x + v.x, y + v.y);
     }
     Vector2 operator-(const Vector2<TYPE>& v) const{
         
         return Vector2<TYPE> (x - v.x, y - v.y);
     }
     Vector2 operator*(TYPE a) const{
         return Vector2<TYPE> (x * a, y * a);
     }
protected:
};

#endif // ST_VECTOR2_H_INCLUDED

