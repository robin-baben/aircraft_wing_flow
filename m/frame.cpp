#include <vector>
#include <algorithm>
#include <math.h>


#include "point3D.h"
#include "frame.h"

using std::vector;

Frame::Frame(const Point3D& _A, const Point3D& _B, const Point3D& _C, const Point3D& _D) {
    Point3D d1 = _A - _C;
    Point3D d2 = _B - _D;

    square = abs(VecProd_Point(d1, d2)) / 2;
    norm = VecProd_Point(d1, d2) / (2*square);


    center = (_A + _B + _C + _D) / 4;
    if (_A == _B || _A == _C || _A == _D)
        center = (_B + _C + _D) / 3;
    else if (_B == _C || _B == _D)
        center = (_A + _C + _D) / 3;
    else if (_C == _D)
        center = (_A + _C + _B) / 3;

    this->points.push_back(_A);
    this->points.push_back(_B);
    this->points.push_back(_C);
    this->points.push_back(_D);
}

Frame::Frame(const vector<Point3D>& p) {
    Point3D d1 = p[0] - p[2];
    Point3D d2 = p[1] - p[3];
    this->triangle = false;

    this->square = abs(VecProd_Point(d1, d2)) / 2;
    this->norm = VecProd_Point(d1, d2) / (2*square);

    this->points.push_back(p[0]);
    for (int i=1; i < 4; ++i) {
        if (p[i] == p[i-1]) {
            this->triangle = true;
        } else {
            this->points.push_back(p[i]);
        }
    }
    if (p[3] == p[0]) {
        this->triangle = true;
    }

    if (this->triangle) {
        this->center = (points[0] + points[1] + points[2]) / 3;
    } else {
        this->center = (points[0] + points[1] + points[2] + points[3]) / 4;
    }
}

double max_z(Frame const& f) {
    vector<double> zs;
    for (Point3D p : f.points) {
        zs.push_back(p.z);
    }
    return *std::max_element(zs.begin(), zs.end());
}

Point3D w_sigma(const Point3D x, Frame sigma) {
    if (!sigma.triangle) {
        return -1 * (Bio_Savar(x, sigma.points[0], sigma.points[1]) +
                    Bio_Savar(x, sigma.points[1], sigma.points[2]) +
                    Bio_Savar(x, sigma.points[2], sigma.points[3]) +
                    Bio_Savar(x, sigma.points[3], sigma.points[0]))  / (4 * M_PI);
    } else {
        return -1 * (Bio_Savar(x, sigma.points[0], sigma.points[1]) +
                    Bio_Savar(x, sigma.points[1], sigma.points[2]) +
                    Bio_Savar(x, sigma.points[2], sigma.points[0]))  / (4 * M_PI);
    }
}