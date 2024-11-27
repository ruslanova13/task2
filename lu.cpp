#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <thread>
#include <vector>


#define MAX(a, b) ((a) > (b) ? (a) : (b)) //максимум из двух элементов
#define MAX3(x, y, eps) (MAX (MAX (x, y), eps)) //максимум из трех элементов
#define ALMOST_EQUAL(x, y, eps) (fabs ((x) - (y)) < ((eps) * MAX3 (x, y, 1.0))) //сравнить на равенство два числа с точностью eps

void decomposition(double * A, double * L, double * U, int n, int *p)
{
    *p = 0;

    for(int i = 0; i < n * n; i++)
    {
        L[i] = 0;
        U[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        L[ i * n + i] = 1;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if( i <= j)
            {
                double summa = 0;
                for(int k = 0; k < i; k++)
                {
                    summa += L[i * n + k] * U[k * n + j];
                }
                U[i * n + j] = A[i * n + j] - summa;
            }

            if (i > j)
            { 
                double summa = 0;
                for(int k = 0; k < j; k++)
                {
                    summa += L[i * n + k] * U[k * n + j];
                }
                if(ALMOST_EQUAL(U[j * n + j], 0, 1e-10))
                {
                    printf("no exist LU-decomposition\n");
                    *p = 1;
                    return;
                } else {
                L[i * n + j] = (A[i * n + j] - summa) / U[j * n + j];
                }
            }     
        }   
    }
}