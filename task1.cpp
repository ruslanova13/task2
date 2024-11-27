#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <cassert>


double *read_square_matrix(const char * filename, int n); // чтение из файла
double * make_matrix(int n, int k); // создание матрицы по формуле
int invers_matrix(int n, double * A, double * X, int *changes); // поиск обратной
void multiply(double * A, double * B, double * X, int n); // умножение матриц
double drift(double * A, double * X, int n);
void decomposition(double * A, double * L, double * U, int n, int *p);
double norm(double * A, int n);

double *copy_matrix(double *m, int n)
{
    double *matrix = (double *)malloc(sizeof(double) * n * n);
    assert(matrix);
    memcpy(matrix, m, sizeof(double) * n * n);
    return matrix;
}

void print_matrix(double *data, int n, int k, FILE *f = stdout) // печать
{

    fprintf(f, "\n");
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            fprintf(f, "%10.3e ", data[i * n + j]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");

}


int main (int argc, const char* argv[] )
{
    int n, m, k;
    const char* filename;
    double *matrix = NULL, *matrix_orig = NULL;
    struct timespec start, finish;
    double elapsed;
    int p = 0;

    if (argc < 4 || argc > 5)
    {
        printf("\nNot right count arguments, should be 4 or 5\n\n"); // если неверное число аргументов 
        return -1;
    }

    n = atoi(argv[1]); // размер матрицы
    m = atoi(argv[2]); // количество выводимых данных
    k = atoi(argv[3]); // номер формулы

    if(k > 0) {
        matrix = make_matrix(n, k);
        matrix_orig = copy_matrix(matrix, n);
    }
    else if (k == 0 && argc != 5)
    {
        printf("\nNo filname, please give me one\n\n"); // если нет файла
        return -1;
    }
    else {
        filename = argv[4];
        matrix = read_square_matrix(filename, n);
        matrix_orig = copy_matrix(matrix, n);
        if (matrix == NULL)
        {
            printf("Матрица не прочиталась\n");
            return -1;
        }
    }
    print_matrix(matrix, n, m, stdout); 

    double *X = (double *)malloc(sizeof(double) * n * n); // присоединенная матрица для обратной 
    int *changes = (int *)malloc(sizeof(int) * n); // вспомогательный массив для сохранения перестановок

    clock_gettime(CLOCK_MONOTONIC, &start);
    invers_matrix(n, matrix, X, changes);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    
    printf("обратная\n");
    print_matrix(X, n, m, stdout); // обратная

    printf("{%d, %f},\n", n, elapsed); //  печать время     
    printf("\nnorm = %10.3e\n", drift(matrix_orig, X, n));
    
    double *L = (double *)malloc(sizeof(double) * n * n);
    assert(L != NULL);
    double *U = (double *)malloc(sizeof(double) * n * n);
    assert(U != NULL);

    clock_gettime(CLOCK_MONOTONIC, &start);
    decomposition(matrix_orig, L, U, n, &p);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    if (p == 1)
    {
        return -1;
    }
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("L\n");
    print_matrix(L, n, m, stdout);
    printf("U\n");
    print_matrix(U, n, m, stdout);
    multiply(L, U, matrix, n);

    for(int i = 0; i < n * n; i++)
    {
        matrix_orig[i] -= matrix[i];
    }
    printf("\nnorm = %10.3e\n", norm(matrix_orig, n));
    printf("{%d, %f},\n", n, elapsed); //  печать время
    
    free(matrix);
    free(matrix_orig);
    free(X);
    free(changes);
    //free(L);
    //free(U);

    return 0;
}