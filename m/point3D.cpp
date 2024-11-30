#include <iostream>
#include <cmath>

#include <vector>
#include <algorithm>


#include "point3D.h"


using std::vector;


Point3D operator*(Point3D const& p, double num) {
    return Point3D(p.x * num, p.y * num, p.z * num);
}

Point3D operator*(double num, Point3D const& p) {
    return Point3D(p.x * num, p.y * num, p.z * num);
}

Point3D operator/(Point3D const& p, double num) {
    return Point3D(p.x / num, p.y / num, p.z / num);
}

std::ostream & operator<<(std::ostream & s, const Point3D & P){
    s << '(' << P.x << ", " << P.y << ", " << P.z << ')';
    return s;
}

Point3D operator+(Point3D const& P1, Point3D const& P2) {
    return Point3D(P1.x + P2.x, P1.y + P2.y, P1.z + P2.z);
}

void operator+=(Point3D& P1, const Point3D& P2) {
    P1.x += P2.x;
    P1.y += P2.y;
    P1.z += P2.z;
}

Point3D operator-(const Point3D& P1, const Point3D& P2) {
        return Point3D(P1.x - P2.x, P1.y - P2.y, P1.z - P2.z);
}

double abs(const Point3D& P) {
    return sqrt(P.x*P.x + P.y*P.y + P.z*P.z);
}

bool operator==(const Point3D& P1, const Point3D& P2) {
    return ( fabs(P1.x-P2.x) < 1e-10 && fabs(P1.y - P2.y) < 1e-10 && fabs(P1.z - P2.z) < 1e-10);
}

bool operator!=(const Point3D& P1, const Point3D& P2) {
    return ( fabs(P1.x-P2.x) > 1e-6 || fabs(P1.y - P2.y) > 1e-6 || fabs(P1.z - P2.z) > 1e-6);
}

double DotProd_Point(const Point3D& P1, const Point3D& P2) { return P1.x * P2.x + P1.y * P2.y + P1.z * P2.z; }

Point3D VecProd_Point(Point3D const & P1, Point3D const& P2)
{
    Point3D res;
    res.x = P1.y * P2.z - P1.z * P2.y;
    res.y = P1.z * P2.x - P1.x * P2.z;
    res.z = P1.x * P2.y - P1.y * P2.x;
    return res;
}

Point3D Bio_Savar(const Point3D& X, const Point3D& P1, const Point3D& P2)
{
    Point3D f;
    Point3D f1 = VecProd_Point(P2 - P1, X - P1) / (DotProd_Point(P2 - P1, P2 - P1) * DotProd_Point(X - P1, X - P1) - DotProd_Point(P2 - P1, X - P1) * DotProd_Point(P2 - P1, X - P1));
    double f2 = DotProd_Point(P2 - P1, X - P2) / abs(X - P2);
    double f3 = DotProd_Point(P2 - P1, X - P1) / abs(X - P1);

    f = f1 * (f2 - f3);
    return f;
}