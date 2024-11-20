#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

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
};
std::ostream & operator<<(std::ostream & s, const Point3D & P){
    s << '(' << P.x << ", " << P.y << ", " << P.z << ')' << std::endl;
    return s;
}

class Frame {
    public:
    Point3D A, B, C, D;

    Frame(const Point3D& _A, const Point3D& _B, const Point3D& _C, const Point3D& _D) {
        this->A = _A;
        this->B = _B;
        this->C = _C;
        this->D = _D;
    }
};

vector<Frame> init(std::string path) {
    vector<Frame> frames;

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
            frames.push_back(Frame(Point3D(x1, y1, z1), Point3D(x2, y2, z2), Point3D(x3, y3, z3), Point3D(x4, y4, z4)));
        }
    }
    in.close();     // закрываем файл
    return frames;
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
    // Point3D P1(1.5, 1.5, 10.3);
    // cout << P1;
    vector<Frame> frames = init("wing_10_20.dat");

    return 0;
}
    