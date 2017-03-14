#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Function that sets the primes
int* setPrimes(int primes[]) {
	int* result = (int*)calloc(NUMPRIMES, sizeof(int));
	for (int i = 0; i < NUMPRIMES; i++) {
		result[i] = primes[i];
	}
	return result;
}

// Function that sets the device primes
int* makeGPUPrimes(int* primes) {
	// Create Pointers to arrays on CPU and GPU
	int* GPUresult = 0;

	// Set memory allocation for CPU and GPU results
	cudaMalloc(&GPUresult, NUMPRIMES*sizeof(int));

	// Copy the primes into the GPUresult
	cudaMemcpy(GPUresult, primes, NUMPRIMES*sizeof(int), cudaMemcpyHostToDevice);

	return GPUresult;
}

// Function to create an initial Poly
Poly makeNewPoly(int coeffArray[], int len, int* primeArray) {
	Poly result;
	result.length = len;
	// Allocate memory
	for (int i = 0; i < NUMPRIMES + 1; i++) {
		result.members[i].coeffs = (int*)calloc(len, sizeof(int));
	}
	// Copy first row of coefficients
	result.members[0].coeffs = coeffArray;
	// Fill the rest of the members
	getMods(result, primeArray);

	return result;
}

// Allocates memory to a Poly with a given length
Poly makePolyGivenLength(int len) {
	Poly result;
	result.length = len;
	// Allocate memory
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(len, sizeof(int));
	}

	return result;
}



// Makes a Poly on the GPU from an existing Poly
// LOOK AT GETTING RID OF POLY ENTIRELY, AND ONLY WORK WITH INPUT ARRAYS
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

// Puts a GPU Poly into a Poly
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

// Copy's an existing Poly into a longer one
Poly copyIntoBigger(Poly a, int len) {
	Poly result;
	result.length  = len;
	// Copy the original data into the new poly
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(len, sizeof(int));
		memcpy(result.members[i].coeffs, a.members[i].coeffs, a.length*sizeof(int));
	}

	return result;
}

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

// Grabs the length from a GPU poly
int lenGrabber(int* d_in) {
	int* holder = (int*)calloc(1, sizeof(int));
	cudaMemcpy(holder, d_in, sizeof(int), cudaMemcpyDeviceToHost);
	return holder[0];
}

// Prints in a format that can imported to Mathematica
void printForReconstruction(Poly g, int* primeArray) {
	printf("mp = {");
	for (int j = 1; j <= NUMPRIMES; j++)
	{
		printf("{%d, ", primeArray[j - 1]);
		for (int i = 0; i < g.length - 1; i++)
		{
			printf("%i,", g.members[j].coeffs[i]);
		}
		printf("%i}", g.members[j].coeffs[g.length - 1]);
		if (j != NUMPRIMES)
			printf(",");
		else
			printf("};");
		printf("\n");
	}
}
