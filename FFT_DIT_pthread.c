#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define POWER 15
#define TRSIZ 32768
#define isign -1
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define NUM_threads 16

/* global variable */
double* data;
double wr[15][32768] = { 0 };
double wi[15][32768] = { 0 };
int n = 2 * TRSIZ;

//pthreads input data for FFT compute
struct thread_data {
	int thread_id;
	int mmax;
	int istep;
	int m;
	int q;
};

struct thread_data thread_data_array[NUM_threads];

//function: calculate the parameter wr & wi
void *calculate_parameter(void *thread_id) {
	
	int m = 0;
	long q = (long) thread_id;
	int mmax = pow(2, q + 1);
	double theta = isign * (6.28318530717959 / mmax);
	double wtemp = sin(0.5 * theta);
	double wpr = -2.0 * wtemp * wtemp;
	double wpi = sin(theta);
	wr[q][0] = 1.0;
	wi[q][0] = 0.0;
	for (m = 1; m < mmax / 2; m++)
	{
		wr[q][m] = wr[q][m - 1] + (wr[q][m - 1] * wpr - wi[q][m - 1] * wpi);
		wi[q][m] = wi[q][m - 1] + (wr[q][m - 1] * wpi + wi[q][m - 1] * wpr);

	}
}

void FFT_compute(void *threadarg) {

	struct thread_data* my_data;
	my_data = (struct thread_data*)threadarg;
	int mmax = my_data->mmax;
	int istep = my_data->istep;
	int q = my_data->q;
	int m = my_data->m;
	int i, j;
	double tempr, tempi;

	for (i = 2 * m - 1; i <= n; i += istep)
	{
		j = i + mmax;
		tempr = wr[q][m - 1] * data[j] - wi[q][m - 1] * data[j + 1];
		tempi = wr[q][m - 1] * data[j + 1] + wi[q][m - 1] * data[j];
		data[j] = data[i] - tempr;
		data[j + 1] = data[i + 1] - tempi;
		data[i] = data[i] + tempr;
		data[i + 1] = data[i + 1] + tempi;
	}

}

int main()
{	
	/* local variable */
	double data1[2 * TRSIZ];
	int num_loop = POWER - 1;
	double wpr, wpi;
	double tempr, tempi;
	int i = 0, j = 0, m = 0, q = 0, t = 0, mmax, istep;
	data = &data1[0] - 1;
	j = 1;
	long taskid[NUM_threads];

	//initial the data
	for (i = 0; i < TRSIZ; i++)
	{
		data1[2 * i] = i;
		data1[2 * i + 1] = 0;
	}

	// Start measuring time
	struct timeval begin, end;
	gettimeofday(&begin, 0);

	// do the bit-reversal
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j + 1], data[i + 1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	// calculate the FFT

	//initial pthreads
	pthread_t threads[NUM_threads];

	//calculate parameter
	//create pthreads
	for (t = 0; t <= num_loop; t++)
	{
		taskid[t] = t;
		pthread_create(&threads[t], NULL, calculate_parameter, (void*)taskid[t]);
	}

	//synchronization
	for (t = 0; t <= num_loop; t++)
		pthread_join(threads[t], NULL);

	for (q = 0; q <= num_loop; q++)
	{
		mmax = pow(2, q + 1);
		int count = mmax / 2;
		istep = mmax << 1;
		for (t = 0; t < NUM_threads; t++) {
			thread_data_array[t].mmax = mmax;
			thread_data_array[t].istep = istep;
			thread_data_array[t].q = q;
		}
		while (count > NUM_threads)
		{
			for (t = 0; t < NUM_threads; t++) {
				thread_data_array[t].thread_id = t;
				thread_data_array[t].m = count - t;
				pthread_create(&threads[t], NULL, FFT_compute, (void*)&thread_data_array[t]);
			}
			for (t = 0; t < NUM_threads; t++)
				pthread_join(threads[t], NULL);
			count -= NUM_threads;
		}
		for (t = 0; t < count; t++) {
			thread_data_array[t].thread_id = t;
			thread_data_array[t].m = count - t;
			pthread_create(&threads[t], NULL, FFT_compute, (void*)&thread_data_array[t]);
		}			
		for (t = 0; t < count; t++)
			pthread_join(threads[t], NULL);

	}

	// Stop measuring time and calculate the elapsed time
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	printf("The runtime is %d seconds, %d micorseconds\n", seconds, microseconds);

	// print the results
	printf("\nFourier components from the DIT algorithm:");
	for (t = 0; t < 20; t += 2)
		printf("\n%f %f", data[t + 1], data[t + 2]);
} // end of dittt()
