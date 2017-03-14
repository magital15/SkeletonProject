#include "kernel.h"
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define TPB 32

// DEVICE WORK FUNCTIONS
__device__
int getRemainderGPU(int coefficient, int mod) {
	coefficient = coefficient % mod;
	coefficient += mod;
	return coefficient % mod;
}

__device__
int getMultGPU(int a, int b) {
	return a * b;
}

__global__
void append(int* d_out, int* d_a, int* d_b, int lenA, int lenB) {
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i == 0)	{
		d_out[0] = lenA;
		return;
	}
	else if (i == lenA*NUMPRIMES + 1) {
		d_out[i] = lenB;
		return;
	}
	if (i > (lenA + lenB)*NUMPRIMES) return;

	const int x = d_a[i];
	const int y = d_b[i];
	d_out[i] = x;
	d_out[i + lenA*NUMPRIMES + 1] = y;
}

__global__
void addMods2(int* d_out, int* d_a, int* d_b, int* primes, int len) {
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i == 0)	{
		d_out[0] = len;
		return;
	}
	if (i > len*NUMPRIMES) return;
	const int r = (i - 1) / len;
	const int mod = primes[r];
	const int x = d_a[i];
	const int y = d_b[i];
	d_out[i] = getRemainderGPU(x + y, mod);
}

__global__
void scalarMultMods2(int* d_out, int* d_in, int scalar, int* primes, int len) {
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i == 0)	{
		d_out[0] = len;
		return;
	}
	if (i > len*NUMPRIMES) return;
	const int r = (i - 1) / len;
	const int mod = primes[r];
	const int x = d_in[i];
	d_out[i] = getRemainderGPU(x*scalar, mod);
}

__global__
void multMods2(int* d_out, int* d_a, int* d_b, int* primes, int lenA, int lenB) {
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const int j = blockIdx.y*blockDim.y + threadIdx.y;
	if (i - j > 3 or i - j < -3) return;
	else if (i == 0 or j == 0) {
		d_out[0] = lenA + lenB - 1;
		return;
	}
	else if (i + j > (lenA + lenB - 1)*NUMPRIMES + 2) {
		return;
	}
	const int r = (i + j) % ((lenA + lenB - 1)*NUMPRIMES + 1);
	const int mod = primes[r];
	const int x = d_a[i];
	const int y = d_b[j];
	printf("[%i][%i] made it\n",i, j);
	int out = getRemainderGPU(getMultGPU(x,y), 3);
	printf("%i is out for [%i][%i]\n", out, i, j);
	d_out[i+j] += getRemainderGPU(x*y, 3);
}

// KERNEL WRAPPER FUNCTIONS
int* addGPU(int* d_a, int* d_b, int* d_primes) {
	// Grab the length
	int lenA = lenGrabber(d_a);
	int lenB = lenGrabber(d_b);
	int len = lenA >= lenB ? lenA : lenB;
	int d_len = len*NUMPRIMES + 1;
	
	// Store data for d_out
	int* d_out = 0;
	cudaMalloc(&d_out, d_len*sizeof(int));

	// Copy smaller GPU Poly into bigger if  needed, then perform Kernel
	if (lenA > lenB) {
		int* d_b2 = 0;
		cudaMalloc(&d_b2, d_len*sizeof(int));
		d_b2 = copyGPUGivenLen(d_b, len);
		addMods2<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b2, d_primes, len);
		cudaFree(d_b2);
	}
	else if (lenA < lenB) {
		int* d_a2 = 0;
		cudaMalloc(&d_a2, d_len*sizeof(int));
		d_a2 = copyGPUGivenLen(d_a, len);
		addMods2<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a2, d_b, d_primes, len);
		cudaFree(d_a2);
	}
	else {
		addMods2<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b, d_primes, len);
	}	

	return d_out;
}

int* scalarModsGPU(int* d_in, int scalar, int* d_primes) {
	// Grab the length
	int len = lenGrabber(d_in);
	int d_len = len*NUMPRIMES + 1;
	
	// Store data for d_out
	int* d_out = 0;
	cudaMalloc(&d_out, d_len*sizeof(int));

	scalarMultMods2<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in, scalar, d_primes, len);

	return d_out;
}

// Faster than making a new function that combines features
int* subtractGPU(int* d_a, int* d_b, int* d_primes) {
	int* d_out = 0;
	d_out = addGPU(d_a, scalarModsGPU(d_b, -1, d_primes), d_primes);
	return d_out;
}

int* multGPU(int* d_a, int* d_b, int* d_primes) {
	// Grab the length
	int lenA = lenGrabber(d_a);
	int lenB = lenGrabber(d_b);
	int len = lenA + lenB - 1;
	int d_len = len*NUMPRIMES + 1;

	// Store data for d_out
	int* d_out = 0;
	cudaMalloc(&d_out, d_len*sizeof(int));

	dim3 dimBlock(32, 32);
	dim3 dimGrid(1, 1);
	multMods2<<<dimGrid, dimBlock>>>(d_out, d_a, d_b, d_primes, lenA, lenB);
	return d_out;	
}

int* appendGPU(int* d_a, int* d_b) {
	int lenA = lenGrabber(d_a);
	int lenB = lenGrabber(d_b);
	int len = lenA + lenB;
	int d_len = len*NUMPRIMES + 2;

	// Set gpu poly memory
	int* d_out = 0;
	cudaMalloc(&d_out, d_len*sizeof(int));
	
	append<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b, lenA, lenB);

	return d_out;
}

