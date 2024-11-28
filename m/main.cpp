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
#include <math.h>
#include <iostream>

#include "point3D.h"
#include "frame.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

void init(
    const std::string path, // путь к файлу
    vector<Frame>& frames, // массив ячеек крыла для заполнения
    vector<Frame>& tr_neib_up, // массив верхних ячеек соседствующих со следом
    vector<Frame>& tr_neib_down, // массив нижних ячеек соседствующих со следом
    //vector<int>& sled_indexes, //массив для хранения индексов ячеек содержащих границу крыла прилегающую к следу
    Point3D& W_st // направление скорости потока
) {

    std::string line;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;

    double alpha = 10 * M_PI / 180;
    double beta = 10 * M_PI / 180;

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
            
            Frame v(vec_points);
            v.ind = size(frames);
            frames.push_back(v);

            if (!frames.back().triangle) { //два раза добавляем одно и то же, надо переделать
                for (Point3D p : vec_points) {
                    if (fabs(p.x - 1.0) < 1e-10) {
                        if (v.norm.y > 0) {
                            tr_neib_up.push_back(Frame(vec_points));
                            tr_neib_up.back().ind = size(frames);
                        } else {
                            tr_neib_down.push_back(Frame(vec_points));
                            tr_neib_down.back().ind = size(frames);
                        }
                        break;
                    }
                }
                
                // vector<Point3D> frame_on_trace;
                // bool flag = false;
                // for (int i = 0; i < 4; ++i) {
                //     if (fabs(vec_points[i].x - 1.0) < 1e-10) {
                //         if(flag) { // тут флаг, чтобы по направлению обхода закинуть точки
                //             frame_on_trace.push_back(Point3D(vec_points[i].x + par, vec_points[i].y, vec_points[i].z));
                //             frame_on_trace.push_back(vec_points[i]);
                //         } else {
                //             frame_on_trace.push_back(vec_points[i]);
                //             frame_on_trace.push_back(Point3D(vec_points[i].x + par, vec_points[i].y, vec_points[i].z));
                //         }
                //         flag = true;
                //     }
                // }
                // if (flag) {
                //     sled.push_back(Frame(frame_on_trace)); // закидываем созданную ячейку в след
                //     sled_indexes.push_back(size(frames)-1); // закидываем индекс соседствующей ячейки из сетки на крыле
                // }
            }
        }

        // for (int i=0; i<sled.size(); ++i){
        //     for(int j=0; j<4; ++j){
        //         cout << sled[i].points[j];
        //     }
        //     cout << endl;
        // }

    }
    in.close();     // закрываем файл
}

void init_trace(
    vector<Frame>& tr_neib_up, 
    vector<Frame>& tr_neib_down,
    vector<Frame>& trace,
    const double par // длина ячейки в вихревом следе
) {
    //отсортируем ячейки по координате z
    std::sort(tr_neib_up.begin(), tr_neib_up.end(), [](const Frame& f1, const Frame& f2) {return max_z(f1) > max_z(f2);});
    std::sort(tr_neib_down.begin(), tr_neib_down.end(), [](const Frame& f1, const Frame& f2) {return max_z(f1) > max_z(f2);});

    for (Frame f : tr_neib_up) {
        vector<Point3D> frame_on_trace;
        bool flag = false;
        for (int i = 0; i < 4; ++i) {
            if (fabs(f.points[i].x - 1.0) < 1e-10) {
                if(flag) { // тут флаг, чтобы по направлению обхода закинуть точки
                    frame_on_trace.push_back(Point3D(f.points[i].x + par, f.points[i].y, f.points[i].z));
                    frame_on_trace.push_back(f.points[i]);
                } else {
                    frame_on_trace.push_back(f.points[i]);
                    frame_on_trace.push_back(Point3D(f.points[i].x + par, f.points[i].y, f.points[i].z));
                }
                flag = true;
            }
        }
        if (flag) {
            trace.push_back(Frame(frame_on_trace)); // закидываем созданную ячейку в след
        }
    }
}

void fill_matrix(
    double* A, //матрица для заполнения
    double* b, // вектор для заполнения
    Point3D& W_st, 
    vector<Frame>& frames, 
    vector<Frame>& trace,
    const vector<Frame>& tr_neib_up, //ячейки соседние сверху
    const vector<Frame>& tr_neib_down) // ячейки соседние снизу
{
    
    int full_size = frames.size() + trace.size() + 1;
    //заполним всю матрицу и весь вектор нулями, можно при инициализации сделать, наверное

    for (int i = 0; i < frames.size(); i++) {
        for (int j = 0; j < frames.size(); j++) {
            if (!frames[j].triangle) {
                A[i + full_size *j] = DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[3]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[3], frames[j].points[0]), frames[i].norm);
            } else {
                A[i + full_size * j] = DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[0]), frames[i].norm);
            }
        }
    }
    cout << 333 << endl;
    

    for (int i = 0; i < trace.size(); i++) {
        //пробегаемся по каждой строке, бурем индекс из up и down
        A[i + frames.size() + full_size * tr_neib_up[i].ind] = -1.0; // на линии отрыва сверху
        A[i + frames.size() + full_size * tr_neib_down[i].ind] = 1.0; // на линии отрыва снизу
    }
    cout << 444 << endl;

    for (int i = 0; i < frames.size(); i++) { // в следе нет треугольников, поэтому без проверки
        for (int j = 0; j < trace.size(); j++) {
            A[i + full_size * (j + frames.size())] = DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[0], trace[j].points[1]), frames[i].norm) +
                DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[1], trace[j].points[2]), frames[i].norm) +
                DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[2], trace[j].points[3]), frames[i].norm) +
                DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[3], trace[j].points[0]), frames[i].norm);
        }

    }
    cout << 555 << endl;

    for (int i = frames.size(); i < full_size - 1; i++) {
        for (int j = frames.size(); j < full_size - 1; j++) {
            if (i == j)
                A[i + full_size * j] = 1.0;
        }
    }
    cout << 666 << endl;

    for (int j = 0; j < frames.size(); j++) { // поменять, для правильного заполнения строк, и правильными ли значениями заполняем
        A[(full_size - 1) + full_size * j] = frames[j].square;
    }
    cout << 777 << endl;

    for (int i = 0; i < frames.size(); i++) { //самый правый столбец матрицы
        A[i + full_size * (full_size - 1)] = 1.0;
    }
    cout << 888 << endl;

    for (int i = 0; i < frames.size(); i++)
        b[i] = DotProd_Point(W_st, frames[i].norm);

    cout << 999 << endl;

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
    Point3D P1(1.0, 0.0, 0.0);
    Point3D P2(1.0, 0.0, 1.0);
    Point3D P3(0.0, 0.0, 1.0);
    Point3D P4(0.0, 0.0, 0.0);
    vector<Point3D> p = {P1, P2, P3, P4}; 

    
    
    // Frame frame1(p);
    // cout << frame1.triangle << endl;

    vector<Frame> frames;
    vector<Frame> tr_neib_up; // массив верхник ячеек соседствующих со следом, его надо дополнительно обработать
    vector<Frame> tr_neib_down; // массив ячеек соседствующих со следом, его надо дополнительно обработать
    vector<Frame> trace;
    //vector<int> trace_neighbours_ind;
    Point3D W_st;
    
    double par = 10.0;

    //
    init("wing_10_20.dat", frames, tr_neib_up, tr_neib_down, W_st);

    init_trace(tr_neib_up, tr_neib_down, trace, par);

    cout << size(frames) << endl;
    cout << size(tr_neib_up) << endl;
    cout << size(tr_neib_down) << endl;

    // for (auto f : trace) {
    //     for (Point3D p : f.points) {
    //         cout << p << ' ';
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    // for (auto f : tr_neib_down) {
    //     cout << f.ind << endl;
    // }

    int full_size = frames.size() + trace.size() + 1;
    double* A = (double*) calloc(full_size*full_size, sizeof(double));
    double* b = (double*) calloc(full_size, sizeof(double));

    fill_matrix(A, b, W_st, frames, trace, tr_neib_up, tr_neib_down);

    // vector<Frame> frames;
    // vector<int> up_ind;
    // vector<int> down_ind;
    // init("wing_10_20.dat", frames, up_ind, down_ind);

    return 0;
}