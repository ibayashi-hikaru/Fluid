#include "Field.h"
#include "FieldUtility.h"
void
Field::initVelocity() {
    for(int i=0; i < height; i++) {
        for(int j=0; j < width; j++) {
            cells[i][j].u0.x() = 0.0;
            cells[i][j].u0.y() = 0.0;
            cells[i][j].u1.x() = 0.0;
            cells[i][j].u1.y() = 0.0;
        } 
    }
}
// フィールドの大きさに合わせて、変換した位置を代入すること。
Vector2d
Field::getVelocity(Vector2d position) const {
    Vector2i nearestIndices = getNearestPointIndices(position);
    if(!isEdge(nearestIndices)) {
        vector<Vector2d> velocities;
        velocities = getSurroundingVelocities(nearestIndices);
        Vector2d local_normalized_position;
        // セルの中心までの距離なので0.5が必要。
        local_normalized_position.x() = (position.x() - (nearestIndices.x() - 0.5) * cellSize) / cellSize;
        local_normalized_position.y() = (position.y() - (nearestIndices.y() - 0.5) * cellSize) / cellSize;
        return FieldUtility::interpolate2d(velocities.at(0), velocities.at(1), velocities.at(2), velocities.at(3),
                                           local_normalized_position);
    } else {
        vector<Vector2d> velocities;
        if(isCorner(nearestIndices)) {
            velocities.push_back(cells.at(width - 1).at(height - 1).u1);        
            velocities.push_back(cells.at(width - 1).at(0).u1);        
            velocities.push_back(cells.at(0).at(height - 1).u1);        
            velocities.push_back(cells.at(0).at(0).u1);        
        }else if(isRightOrLeftSide(nearestIndices)) {
            velocities.push_back(cells.at(width - 1).at(nearestIndices.y() - 1).u1);        
            velocities.push_back(cells.at(width - 1).at(nearestIndices.y()).u1);        
            velocities.push_back(cells.at(0).at(nearestIndices.y() - 1).u1);        
            velocities.push_back(cells.at(0).at(nearestIndices.y()).u1);        
        }else if(isUpOrDownSide(nearestIndices)) {
            velocities.push_back(cells.at(nearestIndices.x() - 1).at(height - 1).u1);        
            velocities.push_back(cells.at(nearestIndices.x() - 1).at(0).u1);        
            velocities.push_back(cells.at(nearestIndices.x()).at(height - 1).u1);        
            velocities.push_back(cells.at(nearestIndices.x()).at(0).u1);        
        }
        Vector2d local_normalized_position;
        // セルの中心までの距離なので0.5が必要。
        local_normalized_position.x() = (position.x() - (nearestIndices.x() - 0.5) * cellSize) / cellSize;
        local_normalized_position.y() = (position.y() - (nearestIndices.y() - 0.5) * cellSize) / cellSize;
        return FieldUtility::interpolate2d(velocities.at(0), velocities.at(1), velocities.at(2), velocities.at(3),
                                           local_normalized_position);
    }    
    std:: cout << "Field::getVelocity error" << std::endl;
    return Vector2d::Zero();
}

// 最も近い格子点を取得する
Vector2i
Field::getNearestPointIndices(Vector2d position) const {
    Vector2d discretePosition;
    Vector2d positionSurplus;
    positionSurplus.x() = fmod(position.x(), cellSize);
    positionSurplus.y() = fmod(position.y(), cellSize);
    if(positionSurplus.x() < cellSize / 2.0) {
        discretePosition.x() = position.x() - positionSurplus.x(); 
    } else {
        discretePosition.x() = position.x() + (cellSize - positionSurplus.x()); 
    }
    if(positionSurplus.y() < cellSize / 2.0) {
        discretePosition.y() = position.y() - positionSurplus.y(); 
    } else {
        discretePosition.y() = position.y() + (cellSize - positionSurplus.y()); 
    }
    return getIndicesOfDiscretePosition(discretePosition);
}

Vector2i
Field::getIndicesOfDiscretePosition(Vector2d discretePosition) const {
    Vector2i indices;
    // double の誤差で僅かに小さくなる時があるので、cellSize/2.0を足してから割り算する。 
    indices.x() = int((discretePosition.x() + cellSize/2.0)/cellSize);
    indices.y() = int((discretePosition.y() + cellSize/2.0)/cellSize);
    return indices;
}

// 格子点からその周りにある定義されてた４つの速度を返す。
vector<Vector2d>
Field::getSurroundingVelocities(Vector2i index) const {
    vector<Vector2d> velocities;
    velocities.push_back(cells.at(index.x() - 1).at(index.y() - 1).u0);
    velocities.push_back(cells.at(index.x() - 1).at(index.y() - 0).u0);
    velocities.push_back(cells.at(index.x() - 0).at(index.y() - 1).u0);
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
Field::periodizedPosition(Vector2d position) const {
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
           tmp_row_x.at(j) = cells.at(i).at(j).u1.x(); 
           tmp_row_y.at(j) = cells.at(i).at(j).u1.y(); 
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
    // Vector2d k0 = getVelocity(periodizedPosition(position));
    // Vector2d k1 = getVelocity(periodizedPosition(position + (-dt/2.0) * k0)); 
    // Vector2d k2 = getVelocity(periodizedPosition(position + (-dt/2.0) * k1); 
    // Vector2d k3 = getVelocity(periodizedPosition(position + (-dt) * k2));
    // return position + (-dt/6.0) * (k0 + 2.0*k1 + 2.0*k2 + k3); 
    return position + (-dt) * getVelocity(periodizedPosition(position));
}

void
Field::addForce(double dt) {
    for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
            cells.at(i).at(j).u1.x() = cells.at(i).at(j).u0.x() + dt*cells.at(i).at(j).force.x(); 
            cells.at(i).at(j).u1.y() = cells.at(i).at(j).u0.y() + dt*cells.at(i).at(j).force.y(); 
        } 
    }
}

void
Field::addTransport(double dt) {
    for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
            Vector2d current_position{(i + 0.5) * cellSize, (j + 0.5) * cellSize};
            Vector2d last_position = traceParticle(current_position, dt);
            cells.at(i).at(j).u1 = getVelocity(periodizedPosition(last_position));
        } 
    }
}

void
Field::addDiffuse(double dt) {
    for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
            complex<double> ikx = complex<double>(0.0, (2.0*PI*i)/width); 
            complex<double> iky = complex<double>(0.0, (2.0*PI*j)/height); 
            ft_vx.at(i).at(j) =  ft_vx.at(i).at(j)/(1.0 - NU * dt * (ikx * ikx + iky * iky));
            ft_vy.at(i).at(j) =  ft_vy.at(i).at(j)/(1.0 - NU * dt * (ikx * ikx + iky * iky));
        } 
    }
}

void
Field::projectField() {
    double inv_w = (double) 1.0/width;   
    double inv_h = (double) 1.0/height;   
    for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
            if(i == 0 && j == 0) {
                ft_vx.at(i).at(j) -= complex<double>(0.0, 0.0);
                ft_vy.at(i).at(j) -= complex<double>(0.0, 0.0);
            } else {
                complex<double> kx = complex<double>(0.0, 2.0*PI*i*inv_w); 
                complex<double> ky = complex<double>(0.0, 2.0*PI*j*inv_h);
                // This variable name is based on the paper "stable fluid"
                complex<double> k_dot_w = kx * ft_vx.at(i).at(j) + ky * ft_vy.at(i).at(j); 
                ft_vx.at(i).at(j) -= (1.0/(kx * kx + ky * ky)) * k_dot_w * kx;
                ft_vy.at(i).at(j) -= (1.0/(kx * kx + ky * ky)) * k_dot_w * ky;
            }
        } 
    }
}

void
Field::swapVelocity() {
    for(int i=0; i<width; i++) {
        for(int j=0; j<height; j++) {
            cells.at(i).at(j).u0 = cells.at(i).at(j).u1;
        } 
    }
}
void
Field::makeSquareForceSource() {
    for(int i = (width*7)/16 + 1; i < (width*9)/16; i++) {
        cells.at(i).at((height*7)/16).force.y() = -5.0;
        cells.at(i).at((height*9)/16).force.y() = 5.0; 
    }
    for(int i = (height*7)/16 + 1; i < (height*9)/16; i++) {
        cells.at((width*7)/16).at(i).force.x() = -5.0;
        cells.at((width*9)/16).at(i).force.x() = 5.0;
    }
    cells.at((width*7)/16).at((height*7)/16).force.x() = -5.0;
    cells.at((width*7)/16).at((height*7)/16).force.y() = -5.0;
    cells.at((width*7)/16).at((height*9)/16).force.x() = -5.0;
    cells.at((width*7)/16).at((height*9)/16).force.y() = 5.0;
    cells.at((width*9)/16).at((height*7)/16).force.x() = 5.0;
    cells.at((width*9)/16).at((height*7)/16).force.y() = -5.0;
    cells.at((width*9)/16).at((height*9)/16).force.x() = 5.0;
    cells.at((width*9)/16).at((height*9)/16).force.y() = 5.0;
}
void
Field::makeLineForceSource() {
    for(int i = (width*7)/16 + 1; i < (width*9)/16; i++) {
        cells.at(i).at(height/2).force.y() = 5.0;
    }
}
void
Field::makeDualForceSource() {
    for(int i = (width*3)/16 + 1; i < (width*5)/16; i++) {
        cells.at(i).at(0).force.y() = 5.0;
    }
    for(int i = (width*6)/16 + 1; i < (width*8)/16; i++) {
        cells.at(i).at(height - 1).force.y() = -5.0;
    }
}
void
Field::resetForceSource() {
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            cells.at(i).at(j).force.x() = 0.0;
            cells.at(i).at(j).force.y() = 0.0;
        }
    }
}

void
Field::updateRot() {
    for(int i = 0; i < width; i++) {
        for(int j = 0; j < height; j++) {
            if(i != 0 && i != width -1 && j != 0 && j != height - 1) {
                div.at(i).at(j) = (cells.at(i + 1).at(j).u1.x() - cells.at(i - 1).at(j).u1.x())/(2.0 * cellSize)
                                    + (cells.at(i).at(j + 1).u1.y() - cells.at(i).at(j - 1).u1.y())/(2.0 * cellSize);
            } else {
                if(i == 0) {
                    if(j == 0) {
                        div.at(0).at(0) = (cells.at(1).at(0).u1.x() - cells.at(width - 1).at(0).u1.x())/(2.0 * cellSize)
                                        + (cells.at(0).at(j + 1).u1.y() - cells.at(0).at(height - 1).u1.y())/(2.0 * cellSize);
                         
                    } else if(j == height - 1) {
                        div.at(0).at(height - 1) = (cells.at(1).at(height - 1).u1.x() - cells.at(width - 1).at(height - 1).u1.x())/(2.0 * cellSize)
                                                 + (cells.at(0).at(0).u1.y() - cells.at(0).at(height - 2).u1.y())/(2.0 * cellSize);
                    } else {
                        div.at(0).at(j) = (cells.at(1).at(j).u1.x() - cells.at(width - 1).at(j).u1.x())/(2.0 * cellSize)
                                       + (cells.at(0).at(j + 1).u1.y() - cells.at(0).at(j - 1).u1.y())/(2.0 * cellSize);
                    }
                }
                if(i == width - 1) {
                    if(j == 0) {
                        div.at(width - 1).at(0) = (cells.at(0).at(0).u1.x() - cells.at(width - 2).at(0).u1.x())/(2.0 * cellSize)
                                    + (cells.at(width - 1).at(1).u1.y() - cells.at(width - 1).at(height - 1).u1.y())/(2.0 * cellSize);
                    } else if(j == height - 1) {
                        div.at(width - 1).at(height - 1) = (cells.at(0).at(height - 1).u1.x() - cells.at(width - 2).at(height - 1).u1.x())/(2.0 * cellSize)
                                                        + (cells.at(width - 1).at(0).u1.y() - cells.at(width - 1).at(height - 2).u1.y())/(2.0 * cellSize);
                    } else {
                        div.at(width - 1).at(j) = (cells.at(0).at(j).u1.x() - cells.at(width - 2).at(j).u1.x())/(2.0 * cellSize)
                                                + (cells.at(width - 1).at(j + 1).u1.y() - cells.at(width - 1).at(j - 1).u1.y())/(2.0 * cellSize);
                    }
                }
                if(j == 0) {
                    if(i != 0 && i != width - 1) {
                        div.at(i).at(0) = (cells.at(i + 1).at(0).u1.x() - cells.at(i - 1).at(0).u1.x())/(2.0 * cellSize)
                                            + (cells.at(i).at(1).u1.y() - cells.at(i).at(height - 1).u1.y())/(2.0 * cellSize);
                    } 
                }
                if(j == height - 1) {
                    if(i != 0 && i != width - 1) {
                        div.at(i).at(height - 1) = (cells.at(i + 1).at(height - 1).u1.x() - cells.at(i - 1).at(height - 1).u1.x())/(2.0 * cellSize)
                                            + (cells.at(i).at(0).u1.y() - cells.at(i).at(height - 2).u1.y())/(2.0 * cellSize);
                    } 
                }
            }
        }
    }
}