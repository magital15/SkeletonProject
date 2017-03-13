#include "kernel.h"
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define TPB 32

// DEVICE WORK FUNCTIONS
__device__
int getRemainderGPU(int coefficient, int mod)
{
	coefficient = coefficient % mod;
	coefficient += mod;
	return coefficient % mod;
}

// WORKS
int* makeGPUPoly(Poly a) {
	// Create Pointers to arrays on CPU and GPU
	int* GPUresult = 0;
	int* CPUresult = 0;

	// Set length of members and total device length
	int len = a.length;
	int d_len = len*NUMPRIMES + 1;

	// Set memory allocation for CPU and GPU results
	CPUresult = (int*)calloc(d_len, sizeof(int));
	cudaMalloc(&GPUresult, (d_len)*sizeof(int));

	// Set the first value in device array to the length
	CPUresult[0] = len;

	int pos = 1;
	// Copy the original data into the CPUresult
	for (int i = 1; i < NUMPRIMES + 1; i++)
	{
		for (int j = 0; j < len; j++)
		{
			CPUresult[pos] = a.members[i].coeffs[j];
			pos++;
		}
	}

	// Copy the CPUresult into the GPUresult
	cudaMemcpy(GPUresult, CPUresult, d_len*sizeof(int), cudaMemcpyHostToDevice);

	free(CPUresult);
	return GPUresult;
}

// WORKS
Poly getGPUPoly(int* d_in) {
	// Grab the length from the device poly
	int len = lenGrabber(d_in);
	int d_len = len*NUMPRIMES + 1;
	
	// Make a new poly based on the length
	Poly result = makePolyGivenLength(len);
	
	// Copy d_in into a CPUresult
	int* CPUresult = (int*)calloc(d_len, sizeof(int));
	cudaMemcpy(CPUresult, d_in, d_len*sizeof(int), cudaMemcpyDeviceToHost);

	int pos = 1;
	for (int i = 1; i <= NUMPRIMES; i++) {
		for (int j = 0; j < len; j++) {
			result.members[i].coeffs[j] = CPUresult[pos];
			pos++;
		}
	}

	cudaFree(d_in);
	return result;
}

// WORKS
int* makeGPUPrimes(int* primes) {
	// Create Pointers to arrays on CPU and GPU
	int* GPUresult = 0;

	// Set memory allocation for CPU and GPU results
	cudaMalloc(&GPUresult, NUMPRIMES*sizeof(int));

	// Copy the primes into the GPUresult
	cudaMemcpy(GPUresult, primes, NUMPRIMES*sizeof(int), cudaMemcpyHostToDevice);

	return GPUresult;
}

// WORKS
__global__
void addMods2(int* d_out, int* d_a, int* d_b, int* primes, int len)
{
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

// WORKS
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

// WORKS
// Grabs the length from a GPU poly
int lenGrabber(int* d_in) {
	int* holder = (int*)calloc(1, sizeof(int));
	cudaMemcpy(holder, d_in, sizeof(int), cudaMemcpyDeviceToHost);
	return holder[0];
}

// WORKS
// Copys a GPU poly into a bigger one
int* copyGPUGivenLen(int* d_in, int len) {
	// Get lengths
	int origLen = lenGrabber(d_in);
	int d_origLen = origLen*NUMPRIMES + 1;
	int d_len = len*NUMPRIMES + 1;

	// Allocate memory for the result
	int* result = 0;
	cudaMalloc(&result, d_len*sizeof(int));

	// Allocate memory for the Host holder
	int* CPUholder = (int*)calloc(d_origLen, sizeof(int));

	// Allocate memory for updated Host
	int* CPUresult = (int*)calloc(d_len, sizeof(int));

	// Copy Device data into CPU holder
	cudaMemcpy(CPUholder, d_in, d_origLen*sizeof(int), cudaMemcpyDeviceToHost);

	// Iterate to get the new CPUresult
	int counter = 0;
	int pos = 1;
	CPUresult[0] = len;
	for (int i = 1; i < d_origLen; i++) {
		if (counter == origLen) {
			for (int j = 0; j < (len - origLen); j++) {
				CPUresult[pos] = 0;
				pos++;
			}
			counter = 0;
		}
		CPUresult[pos] = CPUholder[i];
		pos++;
		counter++;
	}
	// Copy CPU result into GPU result
	cudaMemcpy(result, CPUresult, d_len*sizeof(int), cudaMemcpyHostToDevice);

	free(CPUholder);
	free(CPUresult);
	return result;
}


















