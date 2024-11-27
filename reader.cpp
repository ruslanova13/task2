#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>

#define MAX(x, y) ((x) > (y) ? (x) : (y))
 
double *read_square_matrix(const char * filename, int n)
{
    FILE *f = fopen(filename, "r");
    if(!f) {
        return NULL;
    }

    double *data = (double *)malloc(sizeof(double) * n * n);

    assert(data != NULL);

    for (int i = 0; i < n*n; i++)
    {
        if (fscanf(f, "%lf", data + i) != 1 )
        {
            printf("Не удалось прочитать элемент матрицы\n");
            free(data);
            return NULL;
        }
    }

    fclose(f);

    return data;
}


double f(int n, int k, int i, int j)
{
    switch(k)
    {
        case 1: return n -  MAX(i, j) + 1;
        case 2: return MAX(i, j);
        case 3: return fabs(i - j);
        case 4: return 1.0 / (i + j - 1);
        default: assert(0);
    }
}

double * make_matrix(int n, int k)
{
    double *data = (double *)malloc(sizeof(double) * n * n);
    assert(data != NULL);

    for ( int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            data[i * n + j] = f(n, k, i + 1, j + 1);
        }
    }
    return data;
}