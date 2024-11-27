#ifndef __frame_h__
#define __frame_h__

#include <vector>

#include "point3D.h"

using std::vector;

class Frame {
public:
    vector<Point3D> points; // A, B, C, D;
    double square;
    Point3D norm;
    Point3D center;
    bool triangle;

    Frame() {};

    Frame(const Point3D& _A, const Point3D& _B, const Point3D& _C, const Point3D& _D);

    Frame(const vector<Point3D>& points);
};

#endif