#ifndef GEOMTRY_BASIC_H
#define GEOMTRY_BASIC_H 

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <iostream>

typedef Eigen::Vector3f Vector3f;

struct Triangle {
    Vector3f v0;
    Vector3f v1;
    Vector3f v2;

    Triangle() {
        v0 = v1 = v2 = Vector3f(0.0f);
    }

    Triangle(const Vector3f& p0, const Vector3f& p1, const Vector3f& p2)
      : v0(p0), v1(p1), v2(p2) {}
};


#endif