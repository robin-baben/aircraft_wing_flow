#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;
using std::string;
using std::vector;


extern "C"
{
    extern void zgesv_(int* N, int* Nrhs, std::complex<double>* A, int* ldA, int* Ipvt, std::complex<double> *B, int* ldB, int* info);
}

void print_matrix(std::complex<double> *matrix, int n) {

    for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++){
            printf("%.3lf%+.3lfi ", real(matrix[i*n+j]), imag(matrix[i*n+j]));
        }
       printf("\n");
    }
}

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

    Point3D operator-(Point3D P2) const {
        return Point3D(x - P2.x, y - P2.y, z - P2.z);
    }

    Point3D operator/(double num) const {
        return Point3D(x / num, y / num, z / num);
    }

    Point3D operator*(double num) const {
        return Point3D(x * num, y * num, z * num);
    }
};

double abs(const Point3D P) {
    return sqrt(P.x*P.x + P.y*P.y + P.z*P.z);
}

std::ostream & operator<<(std::ostream & s, const Point3D & P){
    s << '(' << P.x << ", " << P.y << ", " << P.z << ')' << std::endl;
    return s;
}

Point3D operator+(const Point3D& P1, const Point3D& P2) {
        return Point3D(P1.x + P2.x, P1.y + P2.y, P1.z + P2.z);
}

double DotProd_Point(Point3D& P1, Point3D& P2)
{
    double res;
    res = P1.x * P2.x + P1.y * P2.y + P1.z * P2.z;
    return res;
}

Point3D VecProd_Point(Point3D& P1, Point3D& P2)
{
    Point3D res;
    res.x = P1.y * P2.z - P1.z * P2.y;
    res.y = P1.z * P2.x - P1.x * P2.z;
    res.z = P1.x * P2.y - P1.y * P2.x;
    return res;
}



class Frame {
    public:
    Point3D A, B, C, D;
    double square;
    Point3D norm;

    Frame(const Point3D& _A, const Point3D& _B, const Point3D& _C, const Point3D& _D) {
        this->A = _A;
        this->B = _B;
        this->C = _C;
        this->D = _D;

        Point3D d1 = _A - _C;
        Point3D d2 = _B - _D;

        square = abs(VecProd_Point(d1, d2)) / 2;
        norm = VecProd_Point(d1, d2) / (2*square);
    }
};

void init(std::string path, vector<Frame>& frames, vector<int>& up_ind, vector<int>& down_ind) {

    std::string line;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
 
    std::ifstream in(path); // окрываем файл для чтения
    if (in.is_open())
    {
        // while (std::getline(in, line))
        // {
        //     string s;
        //     std::stringstream ss(line);
        //     while (getline(ss, s, " ")) {
        //         // store token string in the vector
        //         
        //     }
        // }
        while (std::getline(in, line)) {
            std::istringstream ss(line);
            ss >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3 >> x4 >> y4 >> z4;
            Point3D A(x1, y1, z1);
            Point3D B(x2, y2, z2);
            Point3D C(x3, y3, z3);
            Point3D D(x4, y4, z4);


            frames.push_back(Frame(A, B, C, D));

            double xs[4] = {x1, x2, x3, x4};
            int count = 0;
            for (int i = 0; i < 4; ++i) {
                if (abs(xs[i] - 1.0) < 1e-6) {
                    count++;
                }
            }

            if (count == 2) {

            }

        }
    }
    in.close();     // закрываем файл
}

    

int main() {
    // srand(time(0));

    // int n = 1000;
    // int Nrhs = 1;
    // int k = 0;

    // while (k < 5) {
    //     std::complex<double>* A =(std::complex<double>*) calloc(n*n ,sizeof(std::complex<double>));
    //     std::complex<double>* A2 =(std::complex<double>*) calloc(n*n ,sizeof(std::complex<double>));
    //     std::complex<double>* b =(std::complex<double>*) calloc(n ,sizeof(std::complex<double>));
    //     std::complex<double>* b2 =(std::complex<double>*) calloc(n ,sizeof(std::complex<double>));
    //     int* Ipvt = (int*) calloc(n, sizeof(int));
    //     int info;

    //     for (int i = 0; i < n; i++) {
    //         b[i] = std::complex<double>((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
    //         b2[i] = b[i];
    //         for (int j = 0; j < n; j++) {
    //             A[j*n+i] = std::complex<double>((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
    //             A2[j*n+i] = A[j*n+i];
    //         }
    //     }

    //     double t = omp_get_wtime();
    //     zgesv_(&n, &Nrhs, A, &n, Ipvt, b, &n, &info);
    //     t = omp_get_wtime() - t;
    //     double err = 0;

    //     for (int i = 0; i < n; ++i) {
    //         std::complex<double> rhs = 0;
    //         for (int j = 0; j < n; ++j) {
    //             rhs += A2[i + j*n] * b[j];
    //         }
    //         err += abs(b2[i] - rhs);
    //     }

    //     std::cout <<"size: " << n << " time: " << t << " error: " << err << std::endl;

    //     n *= 2;
    //     k++;

    //     free(A);
    //     free(A2);
    //     free(b);
    //     free(b2);
    //     free(Ipvt);
    // }
    //Point3D P1(1.5, 1.5, 10.3);
    //Point3D P2(0.5, 3.0, 11.0);
    //cout << P1 - P2;

    vector<Frame> frames;
    vector<int> up_ind;
    vector<int> down_ind;
    init("wing_10_20.dat", frames, up_ind, down_ind);

    return 0;
}
    