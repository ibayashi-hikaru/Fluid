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
     Vector2<TYPE> operator+(const Vector2<TYPE>& v) const{
         return Vector2<TYPE> {x + v.x, y + v.y};
     }
     const Vector2<TYPE>& operator+=(const Vector2<TYPE>& v){
         x += v.x;
         y += v.y;
         return *this;
     }
     Vector2<TYPE> operator-(const Vector2<TYPE>& v) const{
         
         return Vector2<TYPE> {x - v.x, y - v.y};
     }
     const Vector2<TYPE>& operator-=(const Vector2<TYPE>& v){
         x -= v.x;
         y -= v.y;
         return *this;
     }
     Vector2<TYPE> operator*(TYPE a) const{
         return Vector2<TYPE> {x * a, y * a};
     }
     Vector2<TYPE> operator/(TYPE a) const{
         return Vector2<TYPE> {x / a, y / a};
     }
protected:
};
template <typename TYPE> 
Vector2<TYPE> operator*(TYPE a, const Vector2<TYPE>& v){
    return Vector2<TYPE> {a * v.x,  a * v.y};
}

#endif // ST_VECTOR2_H_INCLUDED

