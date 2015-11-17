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
protected:
};

#endif // ST_VECTOR2_H_INCLUDED

