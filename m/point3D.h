#ifndef __point3D_h__
#define __point3D_h__

#include <iostream>


#include <vector>

class Point3D {
public:
    double x, y, z;

    Point3D() {};

    Point3D(const double & _x, const double & _y,  const double & _z):
    x(_x), y(_y), z(_z) {}

    Point3D & operator=(const Point3D & P) {
        this->x=P.x;
        this->y=P.y;
        this->z=P.z;
        return *this;
    }
};

std::ostream & operator<<(std::ostream & s, const Point3D & P);

Point3D operator*(Point3D const& p, double num);

Point3D operator/(Point3D const& p, double num);

Point3D operator+(const Point3D& P1, const Point3D& P2);

Point3D operator-(const Point3D& P1, const Point3D& P2);
double abs(const Point3D& P);

bool operator==(const Point3D& P1, const Point3D& P2);

bool operator!=(const Point3D& P1, const Point3D& P2);

double DotProd_Point(const Point3D& P1, const Point3D& P2);

Point3D VecProd_Point(const Point3D& P1, const Point3D& P2);

Point3D Bio_Savar(const Point3D& X, const Point3D& P1, const Point3D& P2);

double max_z(Point3D const& P);

#endif