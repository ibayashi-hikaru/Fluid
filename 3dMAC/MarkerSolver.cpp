#include "Field.h"

void
Field::UpdateMarkers(double dt) {
    for(size_t i = 0; i < Nx; i++) {
        for(size_t j = 0; j < Ny; j++) {
            for(size_t k = 0; k < Nz; k++) {
                sortedMarkersX[index(i, j, k)] += dt * GetVelocity(sortedMarkersX[index(i, j, k)]);
            }
        }
    }
    sortMarkers();
}

void
Field::initMarkers() {
    sortedMarkersX.clear();
    sortedMarkersX.shrink_to_fit();
    waterDrop(Nx*(5.0/6.0), Ny/2.0, Nz/2.0, Nx/6.0);
    storeWater(1.0/3.0);
    sortMarkers();
}

void
Field::sortMarkers() {
    sort(sortedMarkersX.begin(),
         sortedMarkersX.end(),
         [](const Eigen::Vector3d& a, const Eigen::Vector3d& b){return a.x() < b.x();}
        );
}

bool
Field::existsMarker(size_t cellIndex_x, size_t cellIndex_y, size_t cellIndex_z) {
    if(sortedMarkersX.begin()->x() > (cellIndex_x + 1) * dx) return false;
    if((sortedMarkersX.end() - 1)->x() < cellIndex_x * dx) return false;
    // Returns an iterator pointing to the first element in the range [first,last) which does not compare less than val.
    auto lower_it = lower_bound(sortedMarkersX.begin(),
                                sortedMarkersX.end(),
                                Eigen::Vector3d(cellIndex_x * dx, 0, 0),
                                [](const Eigen::Vector3d& a, const Eigen::Vector3d& b){return a.x() < b.x();}
                                );
    // Returns an iterator pointing to the first element in the range [first,last) which compares greater than val.
    auto upper_it = upper_bound(sortedMarkersX.begin(),
                                sortedMarkersX.end(),
                                Eigen::Vector3d((cellIndex_x + 1) * dx, 0, 0),
                                [](const Eigen::Vector3d& a, const Eigen::Vector3d& b){return a.x() <= b.x();}
                               );
    if(upper_it == lower_it) return false;
    for(auto it = lower_it; it < upper_it; it++) {
        if(   (cellIndex_y) * dx <= it->y() && it->y() <= (cellIndex_y + 1) * dx 
           && (cellIndex_z) * dx <= it->z() && it->z() <= (cellIndex_z + 1) * dx) return true;
    }
    return false;
}
