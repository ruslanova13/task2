#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <thread>
#include <vector>


#define MAX(a, b) ((a) > (b) ? (a) : (b)) //максимум из двух элементов
#define MAX3(x, y, eps) (MAX (MAX (x, y), eps)) //максимум из трех элементов
#define ALMOST_EQUAL(x, y, eps) (fabs ((x) - (y)) < ((eps) * MAX3 (x, y, 1.0))) //сравнить на равенство два числа с точностью eps


void print_matrix(double *data, int n, int k, FILE *f = stdout);
void swap_rows(double * A, int n, int i1, int i2) // перестановка строк
{
    for(int j = 0; j < n; j++)
    {
        double temp = A[i1 * n + j];
        A[i1 * n + j] = A[i2 * n + j];
        A[i2 * n + j] = temp;  
    }
}

void swap_column(double * A, int n, int i1, int i2) // перестановка столбцов
{
    for(int i = 0; i < n; i++)
    {
        double temp = A[i * n + i1];
        A[i * n + i1] = A[i * n + i2];
        A[i * n + i2] = temp;  
    }
}
void scale_row(double * A, int n, int i, double alpha) // домножение i-ой строки на скаляр
{
    for ( int j = 0; j < n; j++) 
    {
        A[i * n + j] *= alpha;
    }
}

void substract(double * A, int n, int i1, int i2, double k) // вычитаем из i2-ой строчки i1-ую строчку, умноженную на k
{
    for(int j = 0; j < n; j++)
    {
        A[i2 * n + j] -= A[i1 * n + j] * k;
    }
}


void substruct_below(double * A, double *X, int n, int j)
{
    int n_workers = std::thread::hardware_concurrency();
    std::vector <std::thread> workers;
    int block_size = (n - (j + 1) + 1) / n_workers + 1;
    if(block_size < 100)
        block_size = 100;

    for(int i = j + 1; i < n; i += block_size)
    {
        workers.push_back(std::thread(
            [=]{
            for(int i2 = i; i2 < n && i2 < i + block_size; ++i2)
            {
                double k = A[i2 * n + j];

                substract(A, n, j, i2, k);
                substract(X, n, j, i2, k);
            }
        }));
    }
    for(int wid = 0; wid < workers.size(); ++wid)
    {
        workers[wid].join();
    }
}

void substruct_above(double * A, double *X, int n, int j)
{
    int n_workers = std::thread::hardware_concurrency();
    std::vector <std::thread> workers;
    int block_size = (j - 1) / n_workers + 1;
    if(block_size < 100)
        block_size = 100;

    for(int i = j - 1; i >= 0; i -= block_size)
    {
        workers.push_back(std::thread([=]{
            for(int i2 = i; i2 >= 0 && i2 > i - block_size; --i2)
            {
                double k = A[i2 * n + j];

                substract(A, n, j, i2, k);
                substract(X, n, j, i2, k);
            }
        }));
    }
    for(int wid = 0; wid < workers.size(); ++wid)
    {
        workers[wid].join();
    }
}


// Вычисляет позицию максимального по модулю значения в столбцах [i, n) строки i.
int argmax(double *A, int n, int i) {
    int n_workers = std::thread::hardware_concurrency();
    std::vector <std::thread> workers;
    std::vector<int> partial(n_workers);
    int block_size = (n - i) / n_workers + 1;
    if(block_size < 100)
        block_size = 100;

    for(int start = i, widx = 0; start < n; start += block_size, ++widx)
    {
        partial.push_back(n);
        workers.push_back(std::thread([=, &partial]{
            double max = 0;
            int argmax_block = start;
            for(int j = start; j < start + block_size && j < n; ++j)
            {
                if(j == start || fabs(A[i * n + j]) > max) {
                    max = fabs(A[i * n + j]);
                    argmax_block = j;
                }
            }
            partial[widx] = argmax_block;
        }));
    }
    for(int wid = 0; wid < workers.size(); ++wid)
    {
        workers[wid].join();
    }
    // Максимум максимумов
    int argmax = partial[0];
    double max = fabs(A[i * n + argmax]);
    for(int w = 1; w < n_workers; w++)
    {
        if(partial[w] < n && fabs(A[i * n + partial[w]]) > max)
        {
            max = fabs(A[partial[w]]);
            argmax = partial[w];
        }
    }
    return argmax;
}


int invers_matrix(int n, double * A, double * X, int *changes) // нахождение обратной
{
    // заполнение единичной матрицы
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            X[i * n + j] = 0;
        }
        X[i * n + i] = 1;
    }

    for( int i = 0; i < n; i++)
    {
        int p = argmax(A, n, i);
        double s = A[i * n + p];
        if( ALMOST_EQUAL(s, 0, 1e-10))
        {
                printf("\n не существует обратной, нулевая строчка\n\n");
                assert(0);
                return -1;
        }

        changes[i] = p;
        if ( p != i)
        {
            swap_column(A, n, p, i);
            
        }
        
        p = i;

        double coef = 1.0 / A[i * n + i];
        scale_row(A, n, i, coef);
        scale_row(X, n, i, coef);
        
        substruct_below(A, X, n, i);

    }

    for ( int i2 = n - 1; i2 > 0; i2--)
    {
        substruct_above(A, X, n, i2);
    }

    for(int i = n - 1; i >= 0; --i)
    {
        if(changes[i] != i) {
            swap_rows(X, n, i, changes[i]);
        }
    }
    return 0;
}
