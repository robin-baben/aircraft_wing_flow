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

#include "point3D.h"
#include "frame.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

void init(
    const std::string path, // путь к файлу
    vector<Frame>& frames, // массив ячеек крыла для заполнения
    vector<Frame>& sled, // массив ячеек вихревого следа для заполнения
    double par, // длина ячейки в вихревом следе
    Point3D& W_st // направление скорости потока
) {

    std::string line;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

    double alpha = 10;
    double beta = 10;

    W_st.x = cos(alpha) * cos(beta);
    W_st.y = sin(alpha) * cos(beta);
    W_st.z = sin(beta);
    
    std::ifstream in(path); // окрываем файл для чтения
    if (in.is_open())
    {
        while (std::getline(in, line)) {
            std::istringstream ss(line);
            ss >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3 >> x4 >> y4 >> z4;
            Point3D A(x1, y1, z1);
            Point3D B(x2, y2, z2);
            Point3D C(x3, y3, z3);
            Point3D D(x4, y4, z4);
            vector<Point3D> vec_points{A, B, C, D};
            // cout << A << B << C << D;
            frames.push_back(Frame(A, B, C, D));

            vector<Point3D> frame_on_trace;
            bool flag = true;
            int counter = 0;
            for (int i = 0; i < 4; ++i) {
                if (abs(vec_points[i].x - 1.0) < 1e-6) {
                    if(flag) { // тут флаг, чтобы по направлению обхода закинуть точки
                        frame_on_trace.push_back(Point3D(vec_points[i].x + par, vec_points[i].y, vec_points[i].z));
                        frame_on_trace.push_back(vec_points[i]);
                    } else {
                        frame_on_trace.push_back(vec_points[i]);
                        frame_on_trace.push_back(Point3D(vec_points[i].x + par, vec_points[i].y, vec_points[i].z));
                    }
                    counter++;
                    flag = true;
                }
            }

            if (counter== 2) {
                // проверка, что это не треугольник == что [0] и [2] точки не совпадают
                if (frame_on_trace[0] != frame_on_trace[2]) {
                    sled.push_back(Frame(frame_on_trace[0], frame_on_trace[1], frame_on_trace[2], frame_on_trace[3]));
                }
           
            }
        }

        // for (int i=0; i<sled.size(); ++i){
        //     for(int j=0; j<4; ++j){
        //         cout << sled[i].four_points[j];
        //     }
        //     cout << endl;
        // }

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
    // Point3D P1(1.0, 0.0, 0.0);
    // Point3D P2(1.0, 0.0, 1.0);
    // Point3D P3(0.0, 0.0, 1.0);
    // Point3D P4(0.0, 0.0, 0.0);
    
    // Frame frame1 (P1, P2, P3, P4);

    vector<Frame> frames;
    vector<Frame> sled;
    Point3D W_st;
    
    double par = 10.0;
    init("wing_10_20.dat", frames, sled, par, W_st);


    cout << frames[0].norm << endl;
    cout << frames[0].square << endl;




    // vector<Frame> frames;
    // vector<int> up_ind;
    // vector<int> down_ind;
    // init("wing_10_20.dat", frames, up_ind, down_ind);

    return 0;
}