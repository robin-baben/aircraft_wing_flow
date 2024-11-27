#include <vector>


#include "point3D.h"
#include "frame.h"

using std::vector;

Frame::Frame(const Point3D& _A, const Point3D& _B, const Point3D& _C, const Point3D& _D) {
        this->four_points.push_back(_A);
        this->four_points.push_back(_B);
        this->four_points.push_back(_C);
        this->four_points.push_back(_D);

        Point3D d1 = _A - _C;
        Point3D d2 = _B - _D;

        Point3D center;

        square = abs(VecProd_Point(d1, d2)) / 2;
        norm = VecProd_Point(d1, d2) / (2*square);
        center = (_A + _B + _C + _D) / 4;
        if (_A == _B || _A == _C || _A == _D)
            center = (_B + _C + _D) / 3;
        else if (_B == _C || _B == _D)
            center = (_A + _C + _D) / 3;
        else if (_C == _D)
            center = (_A + _C + _B) / 3;
    }