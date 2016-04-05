#ifndef ST_FIELD_H_INCLUDED
#define ST_FIELD_H_INCLUDED
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
using T = Eigen::Triplet<double>;
/*
    Class's description
*/

class Field {
    public:
        Field() {
        }
        Field(unsigned long gridNum);
        Field& operator=(const Field& other);
        unsigned long GridNum() const {return Nx;}
        double Dx() const {return dx;}
        void Init();
        void Advect(double dt);
        void AddForce(double dt);      
        void CG_Project(double dt);
        void CG_ProjectWithMarker(double dt);
        void UpdateMarkers(double dt);
        void SetForce(const Eigen::Vector3d& force, const Eigen::Vector3d& position);
        Eigen::Vector3d GetVelocity(const Eigen::Vector3d& position) const;
        std::vector< Eigen::Vector3d> sortedMarkersX;
        void CoutDiv();
        void Extrapolate();
    private:
        const double rho = 1.0;
        const double g = 0.0098;
        unsigned long Nx, Ny, Nz;
        double dx;
        std::vector< std::vector< std::vector< double>>> div;
        std::vector< std::vector< std::vector< double>>> p;
        std::vector< std::vector< std::vector< double>>> ux;
        std::vector< std::vector< std::vector< double>>> uy;
        std::vector< std::vector< std::vector< double>>> uz;
        std::vector< std::vector< std::vector< double>>> xSwap;
        std::vector< std::vector< std::vector< double>>> ySwap;
        std::vector< std::vector< std::vector< double>>> zSwap;
        std::vector< std::vector< std::vector< double>>> forcex;
        std::vector< std::vector< std::vector< double>>> forcey;
        std::vector< std::vector< std::vector< double>>> forcez;
        std::vector<int> allocator;
        std::vector<T> tripletList;
        std::vector<T> newTripletList;
        double getVelocityX(const Eigen::Vector3d& position) const;
        double getVelocityY(const Eigen::Vector3d& position) const;
        double getVelocityZ(const Eigen::Vector3d& position) const;
        void makeBoundary();
        bool isInside(const Eigen::Vector3d& position) const;
        void setForceX(double fx, const Eigen::Vector3d& position);
        void setForceY(double fy, const Eigen::Vector3d& position);
        void setForceZ(double fz, const Eigen::Vector3d& position);
        void clearForce();
        void initVelocity();
        void initPressure();
        void initMarkers();
        Eigen::Vector3d getLastPosition(const Eigen::Vector3d& currentPosition, double dt);
        void addGravityForce(double dt);
        void sortMarkers();
        bool existsMarker(size_t cellIndex_x, size_t cellIndex_y, size_t cellIndex_z);
        unsigned long index(unsigned long i, unsigned long j, unsigned long k) const {return k*(Nx*Ny) + j*Ny + i;}
        Eigen::Vector3d centerPosition(size_t i, size_t j, size_t k) const {
            Eigen::Vector3d centerPosition((i + 0.5) * dx, (j + 0.5) * dx, (k + 0.5) * dx);
            return centerPosition;
        }
        void waterDrop(double x, double y, double z, double radius);
        void storeWater(double rate);
        double getAveVelocityX(int i, int j, int k) const;
        double getAveVelocityY(int i, int j, int k) const;
        double getAveVelocityZ(int i, int j, int k) const;
        void updateVelocityBySwap();
        void copyVelocity();
};

#endif // ST_FIELD_H_INCLUDED
