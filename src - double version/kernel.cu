#include "kernel.h"
#include <stdio.h>
#include <math.h>
#define TPB 32

__device__
double getRemainder(double coefficient, double mod)
{
	coefficient = fmod(coefficient, mod);
	coefficient += mod;
	return fmod(coefficient, mod);
}

__device__
double getMult(double a, double b)
{
	return a * b;
}

__global__
void takeMod(double* d_out, double* d_in, double mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getRemainder(x, mod);
}

__global__
void addMods(double* d_out, double* d_a, double* d_b, double mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_a[i];
    const float y = d_b[i];
	d_out[i] = getRemainder(x + y, mod);
}

__global__
void scalarMultMods(double* d_out, double* d_in, double scalar, double mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getRemainder(x*scalar, mod);
}

__global__
void minusMult(double* d_out, double* d_in, double scalar, double mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getMult(x, scalar);
}

__global__
void mult(double* d_out, double* d_in, double scalar, double mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getRemainder(getMult(x, scalar), mod);
}

void getMods(Poly in, double* primes)
{
	double len = in.length;

	// Declare pointers to device arrays
	double *d_in = 0;
	double *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(double));
	cudaMalloc(&d_out, len*sizeof(double));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[0].coeffs, len*sizeof(double), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		takeMod<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in, 
											  primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(in.members[i].coeffs, d_out, len*sizeof(double), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(d_out);
}

void addPolys(Poly a, Poly b, Poly c, double* primes)
{
	double len;
	if (a.length > b.length)
	{
		len = a.length;
	}
	else
	{
		len = b.length;
	}

	// Declare pointers to device arrays
	double *d_a = 0;
	double *d_b = 0;
	double *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_a, len*sizeof(double));
	cudaMalloc(&d_b, len*sizeof(double));
	cudaMalloc(&d_out, len*sizeof(double));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_a, a.members[i].coeffs, len*sizeof(double), 
				   cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, len*sizeof(double), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b, 
											  primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(double), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
}

void scalarMultPoly(Poly in, Poly out, double scalar, double* primes)
{
	double len = in.length;

	// Declare pointers to device arrays
	double *d_in = 0;
	double *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(double));
	cudaMalloc(&d_out, len*sizeof(double));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[i].coeffs, len*sizeof(double), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		scalarMultMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in,
											scalar, primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(out.members[i].coeffs, d_out, len*sizeof(double), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(d_out);
}

void subtractPolys(Poly a, Poly b, Poly c, double* primes)
{
	scalarMultPoly(b, c, -1, primes);
	addPolys(a, c, c, primes);

//	Above is same as Below?

/*
	int len = a.length;

	// Declare pointers to device arrays
	int *d_a = 0;
	int *d_b = 0;
	int *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_a, len*sizeof(int));
	cudaMalloc(&d_b, len*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_a, a.members[i].coeffs, len*sizeof(int), 
				   cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, len*sizeof(int), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		minusMult<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_b, -1, 
											    primes[i-1]);
		addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_out, 
											  primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(int), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
*/
}

void multiplyPolys(Poly a, Poly b, Poly c, double* primes)
{
	double len = c.length;

	// Declare pointers to device arrays
	double *d_a = 0;
	double *d_b = 0;
	double *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_a, len*sizeof(double));
	cudaMalloc(&d_b, len*sizeof(double));
	cudaMalloc(&d_out, len*sizeof(double));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_a, a.members[i].coeffs, len*sizeof(double), 
				   cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, len*sizeof(double), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		minusMult<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_b, -1, 
											    primes[i-1]);
		addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_out, 
											  primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(double), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
}
