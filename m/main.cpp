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

extern "C"
{
    extern void dgesv_(int* N, int* Nrhs, double* A, int* ldA, int* Ipvt, double* B, int* ldB, int* info);
}

void init(
    const std::string path, // путь к файлу
    vector<Frame>& frames, // массив ячеек крыла для заполнения
    vector<Frame>& tr_neib_up, // массив верхних ячеек соседствующих со следом
    vector<Frame>& tr_neib_down // массив нижних ячеек соседствующих со следом
) {

    std::string line;
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;


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
            v.ind = size(frames)-1;
            frames.push_back(v);

            if (!frames.back().triangle) { //находим верхние и нижние ячейкми на линии отрыва
                for (Point3D p : vec_points) {
                    if (fabs(p.x - 1.0) < 1e-10) {
                        if (v.norm.y > 0) {
                            tr_neib_up.push_back(Frame(vec_points));
                            tr_neib_up.back().ind = size(frames)-1;
                        } else {
                            tr_neib_down.push_back(Frame(vec_points));
                            tr_neib_down.back().ind = size(frames)-1;
                        }
                        break;
                    }
                }
            }
        }
    }
    in.close();     // закрываем файл
}


void init_trace(
    vector<Frame>& tr_neib_up, //массив верхних ячеек
    vector<Frame>& tr_neib_down, // массив нижних ячеек на линии отрыва
    vector<Frame>& trace, //ячейки вихревого следа
    const double par // длина ячейки в вихревом следе
) {
    //отсортируем ячейки по координате z, тогда одному индексу будет соответствовать верхняя и нижняя ячейка, смежные по линии отрыва
    std::sort(tr_neib_up.begin(), tr_neib_up.end(), [](const Frame& f1, const Frame& f2) {return max_z(f1) > max_z(f2);});
    std::sort(tr_neib_down.begin(), tr_neib_down.end(), [](const Frame& f1, const Frame& f2) {return max_z(f1) > max_z(f2);});

    for (Frame f : tr_neib_up) { // так как пробегаемся по верхним ячейка, то нормаль у вихревого следа должна смотреть вверх
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
    Point3D& W_st, //вектор набегающей скорости
    vector<Frame>& frames, 
    vector<Frame>& trace,
    const vector<Frame>& tr_neib_up, //ячейки соседние сверху
    const vector<Frame>& tr_neib_down // ячейки соседние снизу
) {
    int full_size = frames.size() + trace.size() + 1;
    //заполним всю матрицу и весь вектор нулями, можно при инициализации сделать, наверное

    for (int i = 0; i < frames.size(); i++) {
        for (int j = 0; j < frames.size(); j++) {
            if (!frames[j].triangle) {
                A[i + full_size *j] = -1 * (DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[3]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[3], frames[j].points[0]), frames[i].norm)) / (4 * M_PI);
            } else {
                A[i + full_size * j] = -1* (DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[0]), frames[i].norm)) / (4 * M_PI);
            }
        }
    }
    

    for (int i = 0; i < trace.size(); i++) {
        //пробегаемся по каждой строке, бурем индекс из up и down
        A[i + frames.size() + full_size * tr_neib_up[i].ind] = -1.0; // на линии отрыва сверху
        A[i + frames.size() + full_size * tr_neib_down[i].ind] = 1.0; // на линии отрыва снизу
    }

    for (int i = 0; i < frames.size(); i++) { // в следе нет треугольников, поэтому без проверки
        for (int j = 0; j < trace.size(); j++) {
            A[i + full_size * (j + frames.size())] = -1 *  (DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[0], trace[j].points[1]), frames[i].norm) +
                DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[1], trace[j].points[2]), frames[i].norm) +
                DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[2], trace[j].points[3]), frames[i].norm) +
                DotProd_Point(Bio_Savar(frames[i].center, trace[j].points[3], trace[j].points[0]), frames[i].norm)) / (4 * M_PI);
        }

    }

    for (int i = frames.size(); i < full_size - 1; i++) {
        A[i + full_size * i] = 1.0;
    }

    for (int j = 0; j < frames.size(); j++) {
        A[(full_size - 1) + full_size * j] = frames[j].square;
    }

    for (int i = 0; i < frames.size(); i++) { //самый правый столбец матрицы
        A[i + full_size * (full_size - 1)] = 1.0;
    }

    for (int i = 0; i < frames.size(); i++)
        b[i] = -DotProd_Point(W_st, frames[i].norm);
}

void fill_matrix_potential_flow(
    double* A, //матрица для заполнения
    double* b, // вектор для заполнения
    Point3D& W_st, //вектор набегающей скорости
    vector<Frame>& frames
) {
    int full_size = frames.size() + 1;

    for (int i = 0; i < frames.size(); i++) {
        for (int j = 0; j < frames.size(); j++) {
            if (!frames[j].triangle) {
                A[i + full_size *j] = -1 * (DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[3]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[3], frames[j].points[0]), frames[i].norm)) / (4 * M_PI);
            } else {
                A[i + full_size * j] = -1* (DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]), frames[i].norm) +
                    DotProd_Point(Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[0]), frames[i].norm)) / (4 * M_PI);
            }
        }
    }

    for (int i = 0; i < frames.size(); i++) { //самый правый столбец матрицы
        A[i + full_size * (full_size - 1)] = 1.0;
    }

    for (int j = 0; j < frames.size(); j++) { //последняя строка
        A[(full_size - 1) + full_size * j] = frames[j].square;
    }

    for (int i = 0; i < frames.size(); i++)
        b[i] = -DotProd_Point(W_st, frames[i].norm);
}

double square_surface(const std::vector<Frame> &frames) {
    double square = 0.0;
    for (Frame f : frames) {
            square += f.square;
    }
    return square;
} 

Point3D lift_force(
    vector<Frame> frames, // ячейки крыла
    vector<double> g, //коэффициенты прибилжения
    const double rho
) {
    Point3D F(0.0, 0.0, 0.0);

    for (int i = 0; i < frames.size(); ++i) {
        Point3D W;

        if (!frames[i].triangle) {
            W = -1 * (Bio_Savar(frames[i].center, frames[i].points[0], frames[i].points[1])+
                Bio_Savar(frames[i].center, frames[i].points[1], frames[i].points[2]) +
                Bio_Savar(frames[i].center, frames[i].points[2], frames[i].points[3]) +
                Bio_Savar(frames[i].center, frames[i].points[3], frames[i].points[0])) / (4 * M_PI);
        } else {
            W = -1* (Bio_Savar(frames[i].center, frames[i].points[0], frames[i].points[1]) +
                Bio_Savar(frames[i].center, frames[i].points[1], frames[i].points[2]) +
                Bio_Savar(frames[i].center, frames[i].points[2], frames[i].points[0])) / (4 * M_PI);
        }

        cout << frames[i].square * frames[i].norm * 2 * g[i] * g[i] * DotProd_Point(W, W) * rho << endl;
        F += frames[i].square * frames[i].norm * 2 * g[i] * g[i] * DotProd_Point(W, W) * rho;
    }

    return F;
}

void write_frames_to_file(
    const std::string path,
    const vector<Frame> &frames
) {
    std::ofstream out;          // поток для записи
    out.open(path, std::ios::app);      // открываем файл для записи
    if (out.is_open())
    {
        out << endl;
        for (Frame f : frames) {
            for (Point3D p : f.points) {
                out << p << " ";
            }
            out << endl;
        }
    }
    out.close();
}

void write_answer_to_file(
    const std::string path,
    const double* b,
    int size
) {
    std::ofstream out;          // поток для записи
    out.open(path);      // открываем файл для записи
    if (out.is_open())
    {   
        out << size << endl;
        for (int i = 0; i < size; ++i) {
            out << b[i] << std::endl;
        }
    }
    out.close();
}

void write_answer_to_file_point(
    const std::string path,
    const vector<Point3D> velocity) 
{
    std::ofstream out;          // поток для записи
    out.open(path);      // открываем файл для записи
    if (out.is_open())
    {   
        // out << size << endl;
        for (int i = 0; i < velocity.size(); ++i) 
        {
            out << velocity[i] << std::endl;
        }
    }
    out.close();
}

void write_answer_to_file_double_vec(
    const std::string path,
    const vector<double> doubles_vec) 
{
    std::ofstream out;          // поток для записи
    out.open(path);      // открываем файл для записи
    if (out.is_open())
    {   
        out << doubles_vec.size() << endl;
        for (int i = 0; i < doubles_vec.size(); ++i)
            out << doubles_vec[i] << std::endl;
    }
    out.close();
}

void print_matrix(
    const double* A,
    int size
) {
    for(int i=0; i<size; ++i)
    {
        for(int j=0; j< size; ++j)
        {
            if (A[i + j*size] >= 0.0)
                printf(" %.3lf ", A[i + j*size]);
            else
                printf("%.3lf ", A[i + j*size]);
        }
        cout << endl;
    }
}

// Фнукция поиска вектора скоростей в точках коллокациях
vector<Point3D> Velocity_surface(const Point3D& W_start, const vector<double>& g, const vector<Frame>& frames, const vector<Frame>& trace)
{
    vector<Point3D> W;
    
    for(int i=0; i<frames.size(); ++i)
    {
        Point3D part_sum{0.0, 0.0, 0.0};
        for(int j=0; j<frames.size(); ++j)
        {
            if (!frames[j].triangle)
            {
                part_sum += (-1 * (Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]) +
                                   Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]) +
                                   Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[3]) +
                                   Bio_Savar(frames[i].center, frames[j].points[3], frames[j].points[0]))  / (4 * M_PI)) * g[j];
            } 
            else 
            {
                part_sum += (-1 * (Bio_Savar(frames[i].center, frames[j].points[0], frames[j].points[1]) +
                                   Bio_Savar(frames[i].center, frames[j].points[1], frames[j].points[2]) +
                                   Bio_Savar(frames[i].center, frames[j].points[2], frames[j].points[0]))  / (4 * M_PI)) * g[j];
            }
        }
        for(int j=0; j<trace.size(); ++j)
        {
            if (!frames[j].triangle)
            {
                part_sum += (-1 * (Bio_Savar(frames[i].center, trace[j].points[0], trace[j].points[1]) +
                                   Bio_Savar(frames[i].center, trace[j].points[1], trace[j].points[2]) +
                                   Bio_Savar(frames[i].center, trace[j].points[2], trace[j].points[3]) +
                                   Bio_Savar(frames[i].center, trace[j].points[3], trace[j].points[0]))  / (4 * M_PI)) * g[frames.size() + j];
            } 
            else 
            {
                part_sum += (-1 * (Bio_Savar(frames[i].center, trace[j].points[0], trace[j].points[1]) +
                                   Bio_Savar(frames[i].center, trace[j].points[1], trace[j].points[2]) +
                                   Bio_Savar(frames[i].center, trace[j].points[2], trace[j].points[0]))  / (4 * M_PI)) * g[frames.size() + j];
            }
        }
        W.push_back(W_start + part_sum);
    }
    return W;
}

vector<double> Pressure_coeff(const Point3D& W_start, const vector<Point3D>& W)
{
    vector<double> pressure_coeff_vec;

    for(int i=0; i<W.size(); ++i)
        pressure_coeff_vec.push_back(1 - (abs(W[i]) * abs(W[i]))/(abs(W_start) * abs(W_start)) );

    return pressure_coeff_vec;
}




int main() {
    vector<Frame> frames;
    vector<Frame> tr_neib_up; // массив верхник ячеек соседствующих со следом, его надо дополнительно обработать
    vector<Frame> tr_neib_down; // массив ячеек соседствующих со следом, его надо дополнительно обработать
    vector<Frame> trace;
    Point3D W_st;

    double alpha = 10.0 * M_PI / 180;
    double beta = 0.0;
    double rho = 1.0;

    // у вектора набегающей скорости еще модуль должен быть
    W_st.x = cos(alpha) * cos(beta);
    W_st.y = sin(alpha) * cos(beta);
    W_st.z = sin(beta);
    
    double par = 10.0;

    init("wing_10_20.dat", frames, tr_neib_up, tr_neib_down); //перенести расчет W в йункцию инит, чтобы вычисленные W можно было использовать

    init_trace(tr_neib_up, tr_neib_down, trace, par);

    // cout << "squares of triangles:" << endl;
    // for (Frame f : frames) {
    //     if (f.triangle) {
    //         cout << f.square << endl;
    //     }
    // }

    // cout << "indxs frames with zero square:" << endl;
    // for (Frame f : tr_neib_up) {
    //     if (f.square < 1e-10) {
    //         cout << f.ind << endl;
    //     }
    // }

    // for (int i = 430; i < 460; ++i) {
    //     cout << frames[i].norm << endl;
    // }

    // cout << endl;
    // for (int i = 400; i < 430; ++i) {
    //     cout << frames[i].norm << endl;
    // }

    
    // for (Frame f : tr_neib_down) {
    //     cout << f.ind << ' ';
    // }
    // cout << endl;
    // for (Frame f : tr_neib_up) {
    //     cout << f.ind << ' ';
    // }
    // cout << endl;

    // cout << size(frames) << endl;
    // cout << size(tr_neib_up) << endl;
    // cout << size(tr_neib_down) << endl;
    // cout << trace.size() << endl;

    int full_size = frames.size() + trace.size() + 1;
    double* A = (double*) calloc(full_size*full_size, sizeof(double));
    double* b = (double*) calloc(full_size, sizeof(double));
    

    fill_matrix(A, b, W_st, frames, trace, tr_neib_up, tr_neib_down);
    // fill_matrix_potential_flow(A, b, W_st, frames);
    int* Ipvt = (int*) calloc(full_size, sizeof(int));
    int info;
    int Nrhs = 1;

    //print_matrix(A, full_size);
    // for (int i = 0; i < full_size; ++i) {
    //     cout << i << ": " << A[i + (full_size-2) * full_size] << endl;
    // }

    // cout << Bio_Savar(frames[31].center, trace[0].points[0], trace[0].points[1]) << endl;
    // cout << Bio_Savar(frames[31].center, trace[0].points[1], trace[0].points[2]) << endl;
    // cout << Bio_Savar(frames[31].center, trace[0].points[2], trace[0].points[3]) << endl;
    // cout << Bio_Savar(frames[31].center, trace[0].points[3], trace[0].points[0]) << endl;

    // cout << endl;

    // cout << frames[32].center << endl;
    // cout << frames[31].center << endl;
    // cout << trace[0].points[2] << endl;
    // cout << trace[0].points[3] << endl;

    // cout << DotProd_Point(trace[0].points[3] - trace[0].points[2], trace[0].points[3] - trace[0].points[2]) * DotProd_Point(frames[31].center - trace[0].points[2], frames[31].center - trace[0].points[2]) << endl;
    // cout << DotProd_Point(trace[0].points[3] - trace[0].points[2], frames[31].center - trace[0].points[2]) * DotProd_Point(trace[0].points[3] - trace[0].points[2], frames[31].center - trace[0].points[2]) << endl;

    // cout << square_surface(frames) << endl;

    


    dgesv_(&full_size, &Nrhs, A, &full_size, Ipvt, b, &full_size, &info);

    vector<double> b_vec;
    for(int i=0; i<full_size; ++i)
        b_vec.push_back(b[i]); 

    vector<Point3D> W = Velocity_surface(W_st, b_vec, frames, trace);


// Ищем поле скоростей на крыле и заносим в файл.
    // for(int i=0; i < W.size(); ++i)
    //     cout << W[i] << endl;
    // write_answer_to_file_point("velocity_file.gv", W);


// Ищем коэф-ты давления на крыле и заносим в файл
    vector<double> Cp = Pressure_coeff(W_st, W);
    for(int i=0; i < Cp.size(); ++i)
        cout << Cp[i] << endl;
    write_answer_to_file_double_vec("press_coeff.txt", Cp);


    //write_answer_to_file("hello.txt", b, full_size-1);

    // vector<double> g;
    // for (int i = 0; i < frames.size(); ++i) {
    //     g.push_back(b[i]);
    // }

    //Point3D F = lift_force(frames, g, rho);
    //cout << F.y << endl;
    // int n = 100;

    // double* A2 =(double*) calloc(n*n ,sizeof(double));
    // double* A3 =(double*) calloc(n*n ,sizeof(double));
    // double* b2 =(double*) calloc(n ,sizeof(double));
    // double* b3 =(double*) calloc(n ,sizeof(double));

    // for (int i = 0; i < n; i++) {
    //     b2[i] = 1.0;
    //     b3[i] = b2[i];
    //     for (int j = 0; j < n; j++) {
    //         A2[j*n+i] = 1.0*rand()/RAND_MAX;
    //         A3[j*n+i] = A2[j*n+i];
    //     }
    // }


    //dgesv_(&n, &Nrhs, A2, &n, Ipvt, b2, &n, &info);
    // double err = 0;

    // for (int i = 0; i < n; ++i) {
    //     double rhs = 0;
    //     for (int j = 0; j < n; ++j) {
    //         rhs += A3[i + j*n] * b2[j];
    //     }
    //     err += abs(b3[i] - rhs);
    // }
    //cout << info << endl;
    

    return 0;
}