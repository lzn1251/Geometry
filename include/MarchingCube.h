#ifndef GEOMETRY_MARCHING_CUBE_H
#define GEOMETRY_MARCHING_CUBE_H

#include "Basic.h"  

namespace Geometry { 

class MarchingCube {
public:
    MarchingCube();

    std::vector<Triangle> run(int gridSize, float isovalue, float (*func)(Vector3f pos));

private:
    Vector3f LinearInterpolate(const Vector3f& p0, const Vector3f& p1, 
        float value0, float value1, float isolevel);

    static const int edgeTable[256];

    static const int triTable[256][16];

    static const int edgeConnections[12][2];
};

}

#endif