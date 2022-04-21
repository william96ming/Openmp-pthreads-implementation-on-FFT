#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

#define POWER 15
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#define NUM_THREADS 4

double wr[15][32768] = { 0 };
double wi[15][32768] = { 0 };



int main()
{	
	int TRSIZ=pow(2,POWER);
	double wtemp, wpr, wpi, theta;
	double tempr, tempi;
	int num_loop = POWER-1;
	int N = TRSIZ;
	int i = 0, q = 0, j = 0, n = 0, k = 0, m = 0, isign = -1, istep, mmax;
	double data1[2 * TRSIZ];
	double* data;
	data = &data1[0] - 1;
	n = N * 2;
	j = 1;

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
	for (i = 1; i < n; i += 2) 
	{	
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


	//calculate the FFT
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for private (q,mmax,theta,wtemp,wpr,wpi,m)
	for (q = 0; q <= num_loop; q++)
	{
		mmax = pow(2, q + 1);
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr[q][0] = 1.0;
		wi[q][0] = 0.0;
		for (m = 1; m < mmax/2; m++)
		{
			wr[q][m] = wr[q][m - 1] + (wr[q][m - 1] * wpr - wi[q][m - 1] * wpi);
			wi[q][m] = wi[q][m - 1] + (wr[q][m - 1] * wpi + wi[q][m - 1] * wpr);
		}
	}



//#pragma omp parallel for private (mmax,istep,m,i,j,tempr,tempi)
	for (q = 0; q <= num_loop; q++) 
	{
		mmax = pow(2, q + 1);
		istep = mmax << 1;
#pragma omp parallel for private (m,i,j,tempr,tempi)
		for (m = 1; m <= mmax / 2; m++)
		{
			for (i = 2*m-1; i <= n; i += istep) 
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

	}
	


	// Stop measuring time and calculate the elapsed time
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	printf("The runtime of FFT_DIT_omp is %d seconds, %d micorseconds\n", seconds, microseconds);


	// print the results
	printf("\nFourier components from the DIT algorithm:");
	printf("\nPrint the first 10 output:");
	for (k = 0; k < 20; k += 2)
		printf("\n%f %f", data[k + 1], data[k + 2]);
		printf("\n");
} // end of dittt()
