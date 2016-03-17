#ifndef ST_INTERFACE_H_INCLUDED
#define ST_INTERFACE_H_INCLUDED
#include <vector>
#include <list>
#include "Field.h"
enum DrawMode {
    VELOCITY = 0,
    POINTS = 1,
};
class Interface {
    public:
        Interface() {
            lastTime = std::chrono::system_clock::now();
            field = Field(16);    
            lastPosition = Eigen::Vector2d::Zero();
            currentPosition = Eigen::Vector2d::Zero();
            forceSourcePosition = Eigen::Vector3d((field.GridNum() * field.Dx())/2.0,
                                                   field.GridNum() * field.Dx(),
                                                  (field.GridNum() * field.Dx())/2.0);
            saveFlg = false;
            noForceFlg = false;
            imageId = 0;
            elapsedTime = 0.0;
            azimuth = 0.0;
            elevation = 0.0;
            twist = 0.0;
            distance = 5.0;
            startFlg = false;
        }
        std::chrono::system_clock::time_point lastTime;
        double deltaTime; 
        Eigen::Vector2d lastPosition;
        Eigen::Vector2d currentPosition;
        Field field;
        std::vector< std::vector< std::vector<Eigen::Vector3d>>> points;
        DrawMode DRAW_MODE;
        Eigen::Vector3d forceSourcePosition;
        bool saveFlg;
        bool noForceFlg;
        int imageId;
        double elapsedTime;
        double azimuth;
        double elevation;
        double twist;
        double distance;
        int mButton;
        bool startFlg;
};
#endif // ST_INTERFACE_H_INCLUDED
