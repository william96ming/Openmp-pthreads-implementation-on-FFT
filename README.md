# Openmp-pthreads-implementation-on-FFT
Open MP: 
srun --pty -p debug -N 1 -n 1 --cpus-per-task=8 --mem=2Gb /bin/bash 

FFT_DIF_omp.c:  gcc -fopenmp -lm FFT_DIF_omp.c -o FFT_DIF_omp
FFT_DIT_omp.c:  gcc -fopenmp -lm FFT_DIT_omp.c -o FFT_DIF_omp

Pthread: 
srun --partition=short --constraint=cascadelake --cpus-per-task=16 --nodes=1 --pty /bin/bash

FFT_DIT_pthread.c:  gcc -pthread -lm FFT_DIT_pthread.c -o FFT_DIT_pthread
