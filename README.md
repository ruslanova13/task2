для запуска нужно ввести команду 

./a.out n_workers n m k 

если k = 0, то дополнительно указывается имя файла, из которого считывается матрица

* n_workers - число потоков
* n - размер матрицы
* m - размер выводимого минора
* k - номер формулы

например

./a.out 10 2000 3 2
