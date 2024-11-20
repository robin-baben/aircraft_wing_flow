#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <iostream>
#include <omp.h>


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

int main() {
    srand(time(0));

    int n = 1000;
    int Nrhs = 1;
    int k = 0;

    while (k < 5) {
        std::complex<double>* A =(std::complex<double>*) calloc(n*n ,sizeof(std::complex<double>));
        std::complex<double>* A2 =(std::complex<double>*) calloc(n*n ,sizeof(std::complex<double>));
        std::complex<double>* b =(std::complex<double>*) calloc(n ,sizeof(std::complex<double>));
        std::complex<double>* b2 =(std::complex<double>*) calloc(n ,sizeof(std::complex<double>));
        int* Ipvt = (int*) calloc(n, sizeof(int));
        int info;

        for (int i = 0; i < n; i++) {
            b[i] = std::complex<double>((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
            b2[i] = b[i];
            for (int j = 0; j < n; j++) {
                A[j*n+i] = std::complex<double>((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
                A2[j*n+i] = A[j*n+i];
            }
        }

        double t = omp_get_wtime();
        zgesv_(&n, &Nrhs, A, &n, Ipvt, b, &n, &info);
        t = omp_get_wtime() - t;
        double err = 0;

        for (int i = 0; i < n; ++i) {
            std::complex<double> rhs = 0;
            for (int j = 0; j < n; ++j) {
                rhs += A2[i + j*n] * b[j];
            }
            err += abs(b2[i] - rhs);
        }

        std::cout <<"size: " << n << " time: " << t << " error: " << err << std::endl;

        n *= 2;
        k++;

        free(A);
        free(A2);
        free(b);
        free(b2);
        free(Ipvt);
    }
    return 0;
}