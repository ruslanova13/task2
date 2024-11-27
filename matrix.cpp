#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <thread>


void multiply(double * A, double * B, double * X, int n)
{
    int n_workers = std::thread::hardware_concurrency();
    std::vector <std::thread> workers;
    int block_size = n / n_workers + 1;

    for(int start = 0; start < n; start += block_size)
    {
        workers.push_back(std::thread([=]{
            for(int i = start; i < n && i < start + block_size; ++i)
            {
                for(int j = 0; j < n; j++)
                {
                    double Xij = 0;
                    Xij = 0;
                    for(int p = 0; p < n; p++)
                    {
                        Xij += A[ i * n + p] * B[p * n + j];
                    }
                    X[i * n + j] = Xij;
                }
            }
        }));
    }
    for(int wid = 0; wid < workers.size(); ++wid)
    {
        workers[wid].join();
    }
}

double norm(double * A, int n)
{
    double s = 0;
    for(int i = 0; i < n*n; i++)
    {
        s += (A[i] * A[i]);
    }

    return sqrt(s);
}

double drift(double * A, double * X, int n)
{   
    double *B = (double *)malloc(sizeof(double) * n * n); // присоединенная матрица для обратной 
    multiply(A, X, B, n); // A -> A, X -> A^-1 B -> A * A^-1

    for(int i = 0; i < n; i++)
    {
        B[i * n + i] -= 1.0;
    }

    return norm(B, n);
}