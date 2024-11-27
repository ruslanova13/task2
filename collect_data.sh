#!/bin/bash



# Аня: ./collect_data.sh  ./a.out >anya.dat
# Ульяна: ./collect_data.sh ../../MurMur_UA/task1/a.out >ua.dat

PROG=$1  # Это аргумент командной строки (имя выполняемого файла для тестирования)
MAT_SIZES="4, 8, 16, 32, 64, 512, 1024, 1500, 1700, 2048, 2100, 2500, 3000, 3500, 3700"
MINOR_SIZE=4
FORMULA=3

echo "{"
for size in ${MAT_SIZES}
do
    ${PROG} ${size} ${MINOR_SIZE} ${FORMULA} | grep "{"
done
echo "}"