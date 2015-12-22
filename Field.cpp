#include "Field.h"
#include "FieldUtility.h"
// フィールドの大きさに合わせて、変換した位置を代入すること。
Vector2d
Field::getVelocity(Vector2d position) const {
    Vector2d nearestDiscrete = getNearestDiscretePosition(position);
    Vector2i nearestPositionIndices = getIndicesOfDiscretePosition(nearestDiscrete); 
    if(!isEdge(nearestPositionIndices)) {
        vector<Vector2d> velocities;
        velocities = getSurroundingVelocities(nearestDiscrete);
        Vector2d local_normalized_position;
        // セルの中心までの距離なのでcellSize/2.0を足す。
        local_normalized_position.x() = (position.x() - nearestDiscrete.x() + cellSize/2.0)/cellSize;
        local_normalized_position.y() = (position.y() - nearestDiscrete.y() + cellSize/2.0)/cellSize;
        return FieldUtility::interpolate2d(velocities.at(0), velocities.at(1), velocities.at(2), velocities.at(3),
                                           local_normalized_position);
    } else {
        vector<Vector2d> velocities;
        if(isCorner(nearestPositionIndices)) {
            velocities.push_back(cells.at(width - 1).at(height - 1).u0);        
            velocities.push_back(cells.at(width - 1).at(0).u0);        
            velocities.push_back(cells.at(0).at(height - 1).u0);        
            velocities.push_back(cells.at(0).at(0).u0);        
        }else if(isRightOrLeftSide(nearestPositionIndices)) {
            velocities.push_back(cells.at(width - 1).at(nearestPositionIndices.y() - 1).u0);        
            velocities.push_back(cells.at(width - 1).at(nearestPositionIndices.y()).u0);        
            velocities.push_back(cells.at(0).at(nearestPositionIndices.y() - 1).u0);        
            velocities.push_back(cells.at(0).at(nearestPositionIndices.y()).u0);        
        }else if(isUpOrDownSide(nearestPositionIndices)) {
            velocities.push_back(cells.at(nearestPositionIndices.x() - 1).at(height - 1).u0);        
            velocities.push_back(cells.at(nearestPositionIndices.x() - 1).at(0).u0);        
            velocities.push_back(cells.at(nearestPositionIndices.x()).at(height - 1).u0);        
            velocities.push_back(cells.at(nearestPositionIndices.x()).at(0).u0);        
        }
        Vector2d local_normalized_position;
        // セルの中心までの距離なのでcellSize/2.0を足す。
        local_normalized_position.x() = (position.x() - nearestDiscrete.x() + cellSize/2.0)/cellSize;
        local_normalized_position.y() = (position.y() - nearestDiscrete.y() + cellSize/2.0)/cellSize;
        return FieldUtility::interpolate2d(velocities.at(0), velocities.at(1), velocities.at(2), velocities.at(3),
                                           local_normalized_position);
    }    
    std:: cout << "Field::getVelocity error" << std::endl;
    return Vector2d::Zero();
}

// 最も近い格子点の位置を取得する(doubleなので正確ではない)
Vector2d
Field::getNearestDiscretePosition(Vector2d position) const {
    Vector2d discretePosition;
    Vector2d positionSurplus;
    positionSurplus.x() = fmod(position.x(), cellSize);
    positionSurplus.y() = fmod(position.y(), cellSize);
    if(positionSurplus.x() < cellSize) {
        discretePosition.x() = position.x() - positionSurplus.x(); 
    } else {
        discretePosition.x() = position.x() + (cellSize - positionSurplus.x()); 
    }
    if(positionSurplus.y() < cellSize) {
        discretePosition.y() = position.y() - positionSurplus.y(); 
    } else {
        discretePosition.y() = position.y() + (cellSize - positionSurplus.y()); 
    }
    return discretePosition;
}

Vector2i
Field::getIndicesOfDiscretePosition(Vector2d discretePosition) const {
    Vector2i indices;
    // 流石に雑過ぎかも。。。
    indices.x() = int((discretePosition.x() + cellSize/2.0)/cellSize);
    indices.y() = int((discretePosition.y() + cellSize/2.0)/cellSize);
    return indices;
}

// 格子点の位置からその周りにある定義されてた４つの速度を返す。
vector<Vector2d>
Field::getSurroundingVelocities(Vector2d discretePosition) const {
    Vector2i index(getIndicesOfDiscretePosition(discretePosition));
    vector<Vector2d> velocities;
    velocities.push_back(cells.at(index.x() - 1).at(index.y() - 1).u0);
    velocities.push_back(cells.at(index.x() - 0).at(index.y() - 1).u0);
    velocities.push_back(cells.at(index.x() - 1).at(index.y() - 0).u0);
    velocities.push_back(cells.at(index.x() - 0).at(index.y() - 0).u0);
    return velocities;
}

bool
Field::isEdge(Vector2i positionIndices) const {
    return positionIndices.x() == 0 ||
           positionIndices.y() == 0 ||
           positionIndices.x() == width ||
           positionIndices.y() == height;
}

bool
Field::isCorner(Vector2i positionIndices) const {
    return (positionIndices.x() == 0 && positionIndices.y() == 0) ||
           (positionIndices.x() == 0 && positionIndices.y() == height) ||
           (positionIndices.x() == width && positionIndices.y() == 0) ||
           (positionIndices.x() == width && positionIndices.y() == height);
}

bool
Field::isRightOrLeftSide(Vector2i positionIndices) const {
    return (positionIndices.x() == 0 || positionIndices.x() == width)
           && positionIndices.y() != 0
           && positionIndices.y() != height;
}
bool
Field::isUpOrDownSide(Vector2i positionIndices) const {
    return (positionIndices.y() == 0 || positionIndices.y() == height)
           && positionIndices.x() != 0
           && positionIndices.x() != width;
}

Vector2d
Field::periodizePosition(Vector2d position) const {
    Vector2d periodizedPosition;
    periodizedPosition.x() = fmod(position.x(), width * cellSize);
    periodizedPosition.y() = fmod(position.y(), height * cellSize);
    return periodizedPosition;
}

void
Field::FFT2d() {
    vector< vector<double>> in_rows_x(height, vector<double>(width, 0.0));
    vector< vector<double>> in_rows_y(height, vector<double>(width, 0.0));
    for(int i=0; i < height; i++) {
        vector<double> tmp_row_x(width), tmp_row_y(width);
        for(int j=0; j < width; j++) {
           tmp_row_x.at(j) = cells[i][j].u1.x(); 
           tmp_row_y.at(j) = cells[i][j].u1.y(); 
        }
        in_rows_x.at(i) = tmp_row_x;
        in_rows_y.at(i) = tmp_row_y;
    }
    FFT<double> fft;
    vector< vector< complex<double>>> med_x(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_y(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    for (int i = 0; i < height; i++) {
        fft.fwd(med_x.at(i), in_rows_x.at(i));
        fft.fwd(med_y.at(i), in_rows_y.at(i));
    }
    vector< vector< complex<double>>> med_cols_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_cols_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for(int i=0; i < width; i++) {
        for(int j=0; j < height; j++) {
            med_cols_x.at(i).at(j) = med_x.at(j).at(i); 
            med_cols_y.at(i).at(j) = med_y.at(j).at(i); 
        }
    }
    vector< vector< complex<double>>> out_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> out_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for (int i = 0; i < width; i++) {
        fft.fwd(out_x.at(i), med_cols_x.at(i));
        fft.fwd(out_y.at(i), med_cols_y.at(i));
    }
    for(int i=0; i < height; i++) {
        for(int j=0; j < width; j++) {
            ft_vx.at(i).at(j) = out_x.at(j).at(i); 
            ft_vy.at(i).at(j) = out_y.at(j).at(i);
        }
    }
}
void
Field::invFFT2d() {
    FFT<double> fft;
    vector< vector< complex<double>>> med_x(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_y(height, vector< complex<double>>(width, complex<double>(0.0, 0.0)));
    for(int i=0; i < height; i++) {
        fft.inv(med_x.at(i), ft_vx.at(i));
        fft.inv(med_y.at(i), ft_vy.at(i));
    }
    vector< vector< complex<double>>> med_cols_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> med_cols_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for(int i=0; i < width; i++) {
        for(int j=0; j < height; j++) {
            med_cols_x.at(i).at(j) = med_x.at(j).at(i); 
            med_cols_y.at(i).at(j) = med_y.at(j).at(i); 
        }
    }
    vector< vector< complex<double>>> out_x(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    vector< vector< complex<double>>> out_y(width, vector< complex<double>>(height, complex<double>(0.0, 0.0)));
    for (int i = 0; i < width; i++) {
        fft.inv(out_x.at(i), med_cols_x.at(i));
        fft.inv(out_y.at(i), med_cols_y.at(i));
    }
    for(int i=0; i < height; i++) {
        for(int j=0; j < width; j++) {
            cells.at(i).at(j).u1.x() = out_x.at(j).at(i).real(); 
            cells.at(i).at(j).u1.y() = out_y.at(j).at(i).real();
        }
    }

}
// Runge-Kutta
Vector2d
Field::traceParticle(Vector2d position, double dt) const {
    Vector2d k0 = dt * getVelocity(periodizePosition(position));
    Vector2d k1 = dt * getVelocity(periodizePosition(position - k0/2.0)); 
    Vector2d k2 = dt * getVelocity(periodizePosition(position - k1/2.0)); 
    Vector2d k3 = dt * getVelocity(periodizePosition(position - k2));
    return position - (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0; 
}

void
Field::addForce(double dt) {
    if(CALC_STEP == STEP0) {
        for(int i=0; i<width; i++) {
            for(int j=0; j<height; j++) {
                cells[i][j].u1.x() = cells[i][j].u0.x() + dt*cells[i][j].force.x(); 
                cells[i][j].u1.y() = cells[i][j].u0.y() + dt*cells[i][j].force.y(); 
            } 
        }
        CALC_STEP = STEP1;
    } else {
        std::cout << "Force cannot be applied at this step." << std::endl; 
    }    
}

void
Field::addTransport(double dt) {
    if(CALC_STEP == STEP1) {
        for(int i=0; i<width; i++) {
            for(int j=0; j<height; j++) {
                Vector2d current_position{(i + 0.5) * cellSize, (j + 0.5) * cellSize};
                Vector2d last_position = periodizePosition(traceParticle(current_position, dt));
                cells[i][j].u1 += getVelocity(last_position) - cells[i][j].u0;
            } 
        }
        CALC_STEP = STEP2;
    } else {
        std::cout << "Transport cannot be applied at this step." << std::endl; 
    }
}

void
Field::addDiffuse(double dt) {
    if(CALC_STEP == STEP2) {
        for(int i=0; i<width; i++) {
            for(int j=0; j<height; j++) {
                complex<double> ikx = complex<double>(0.0, (2.0*PI*i)/width); 
                complex<double> iky = complex<double>(0.0, (2.0*PI*j)/height); 
                ft_vx.at(i).at(j) =  ft_vx.at(i).at(j)/(1.0 - NU * dt * (ikx * ikx + iky * iky));
                ft_vy.at(i).at(j) =  ft_vy.at(i).at(j)/(1.0 - NU * dt * (ikx * ikx + iky * iky));
            } 
        }
        CALC_STEP = STEP3;
    } else {
        std::cout << "Diffuse cannot be applied at this step." << std::endl; 
    }
}

void
Field::projectField() {
    if(CALC_STEP == STEP3) {
        double inv_w = (double) 1.0/width;   
        double inv_h = (double) 1.0/height;   
        for(int i=0; i<width; i++) {
            for(int j=0; j<height; j++) {
                if(i == 0 && j == 0) {
                    ft_vx.at(i).at(j) -= 0.0;
                    ft_vy.at(i).at(j) -= 0.0;
                } else {
                    complex<double> ikx = complex<double>(0.0, 2.0*PI*i*inv_w); 
                    complex<double> iky = complex<double>(0.0, 2.0*PI*j*inv_h);
                    double ik2 = -((2.0*PI*i*inv_w) * (2.0*PI*i*inv_w) + (2.0*PI*j*inv_h) * (2.0*PI*j*inv_h));
                    // This variable name is based on the paper "stable fluid"
                    complex<double> ik_dot_w = ikx * ft_vx.at(i).at(j) + iky * ft_vy.at(i).at(j); 
                    ft_vx.at(i).at(j) -= (1.0/ik2) * ik_dot_w * ikx;
                    ft_vy.at(i).at(j) -= (1.0/ik2) * ik_dot_w * iky;
                }
            } 
        }
        CALC_STEP = STEP4;
    } else {
        std::cout << "Project cannot be applied at this step." << std::endl; 
    }
}

void
Field::swapVelocity() {
    for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
            cells[i][j].u0 = cells[i][j].u1;
        } 
    }
}
