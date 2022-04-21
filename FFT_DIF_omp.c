#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define TRSIZ 32768
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define NUM_THREADS 5

double wr[15][32768] = { 0 };
double wi[15][32768] = { 0 };

int main()
{
	double wtemp, wpr, wpi, theta;
	double tempr,tempi;
	int N = TRSIZ;
	int i = 0, n = 0, j = 0, x = 0, k = 0, q = 0, m = 0, isign = -1, istep, mmax;
	double data1[2 * TRSIZ];
	double* data;
	data = &data1[0] - 1;
	n = N * 2;
//	mmax = n / 2;

	//initial the data
	for (i = 0; i < TRSIZ; i++)
	{
		data1[2 * i] = i;
		data1[2 * i + 1] = 0;
	}

	// Start measuring time
	struct timeval begin, end;
	gettimeofday(&begin, 0);

	omp_set_num_threads(NUM_THREADS);

	// calculate the FFT
	// calculate the paramter wi, wr
#pragma omp parallel for private (q,mmax,theta,wtemp,wpr,wpi,m)
	for (q = 14; q >= 0; q--) {
		mmax = pow(2, q+1);
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr[q][0] = 1.0;
		wi[q][0] = 0.0;
		for (m = 1; m < mmax / 2; m++)
		{
			wtemp = wr[q][m - 1];
			wr[q][m] = wr[q][m - 1] + wtemp * wpr - wi[q][m - 1] * wpi;
			wi[q][m] = wi[q][m - 1] + wtemp * wpi + wi[q][m - 1] * wpr;
		}
	}
	
	for (q = 14; q >= 0; q--) {	
		mmax = pow(2, q+1);
		istep = mmax << 1;
#pragma omp parallel for private (i,m,tempi,j,tempr)		
		for (m = 1; m <= mmax / 2; m++)
		{
			for (i = 2 * m - 1; i <= n; i += istep) {
				j = i + mmax;
				tempr = data[i];
				tempi = data[i + 1];
				data[i] = data[i] + data[j];
				data[i + 1] = data[i + 1] + data[j + 1];
				data[j] = tempr - data[j];
				data[j + 1] = tempi - data[j + 1];
				tempr = data[j];
				data[j] = wr[q][m - 1] * data[j] - wi[q][m - 1] * data[j + 1];
				data[j + 1] = wr[q][m - 1] * data[j + 1] + wi[q][m - 1] * tempr;

			}
		}
	}

	// do the bit-reversal
	x = 1;
	for (i = 1; i < n; i += 2) {
		if (x > i) {
			SWAP(data[x], data[i]);
			SWAP(data[x + 1], data[i + 1]);
		}
		m = n >> 1;
		while (m >= 2 && x > m) {
			x -= m;
			m >>= 1;
		}
		x += m;
	}

	// Stop measuring time and calculate the elapsed time
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	printf("The runtime is %d seconds, %d micorseconds\n", seconds, microseconds);

	// print the results
	printf("\nFourier components from the DIT algorithm:");
	for (k = 0; k < 20; k += 2)
		printf("\n%f %f", data[k + 1], data[k + 2]);
} // end of diftt()