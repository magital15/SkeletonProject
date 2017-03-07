#include "kernel.h"
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define TPB 32

__device__
int getRemainder(int coefficient, int mod)
{
	coefficient = coefficient % mod;
	coefficient += mod;
	return coefficient % mod;
}

__device__
int getMult(int a, int b)
{
	return a * b;
}

__global__
void takeMod(int* d_out, int* d_in, int mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getRemainder(x, mod);
}

__global__
void addMods(int* d_out, int* d_a, int* d_b, int mod, int size)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i > size) return;
	const float x = d_a[i];
    const float y = d_b[i];
	d_out[i] = getRemainder(x + y, mod);

	//printf("d_out[%d] = %d \n", i, d_out[i]);
}

__global__
void scalarMultMods(int* d_out, int* d_in, int scalar, int mod, int size)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i > size) return;
	const float x = d_in[i];
	d_out[i] = getRemainder(x*scalar, mod);
}

// expects d_out to be large enough to hold i+monomial elements
// not sure if works for negative scalars
__global__
void monomialScalarMultMods(int* d_out, int* d_in, int scalar, int monomial, int mod, int size)
{
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
		if (i > size) return;
        const float x = d_in[i];
        d_out[i + monomial] = getRemainder(x*scalar, mod);
		
		//printf("d_temp[%d] = %d \n", i, d_out[i]);
}

__global__
void minusMult(int* d_out, int* d_in, int scalar, int mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getMult(x, scalar);
}

void getMods(Poly in, int* primes)
{
	int len = in.length;

	// Declare pointers to device arrays
	int *d_in = 0;
	int *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[0].coeffs, len*sizeof(int), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		takeMod<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in, 
											  primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(in.members[i].coeffs, d_out, len*sizeof(int), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(d_out);
}

void addPolys(Poly a, Poly b, Poly c, int* primes)
{
	int len;
	if (a.length > b.length)
	{
		len = a.length;
		b = copyIntoBigger(b, a.length);
	}
	else if (a.length == b.length)
	{
		len = a.length;
	}
	else
	{
		len = b.length;
		a = copyIntoBigger(a, b.length);
	}
	

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
		addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b, 
											  primes[i-1], len);

		// reconstruct answer

		// Copy results from device to host
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(int), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
}

void scalarMultPoly(Poly in, Poly out, int scalar, int* primes)
{
	int len = in.length;

	// Declare pointers to device arrays
	int *d_in = 0;
	int *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES+1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[i].coeffs, len*sizeof(int), 
				   cudaMemcpyHostToDevice);
  
		// Launch kernel to compute and store modded polynomial values
		scalarMultMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in,
											scalar, primes[i-1], len);

		// Copy results from device to host
		cudaMemcpy(out.members[i].coeffs, d_out, len*sizeof(int), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(d_out);
}

void subtractPolys(Poly a, Poly b, Poly c, int* primes)
{
	scalarMultPoly(b, c, -1, primes);
	addPolys(a, c, c, primes);
}

void multiplyPolys(Poly a, Poly b, Poly c, int* primes)
{
	int len = c.length;
	Poly shorterPoly = a.length <= b.length ? a : b;
	Poly longerPoly = a.length <= b.length ? b : a;

	// Declare pointers to device arrays
	int *d_a = 0;
	int *d_b = 0;
	int *d_out = 0;
	int *d_temp = 0;

	//Allocate memory for device arrays

	cudaMalloc(&d_a, longerPoly.length*sizeof(int));
	cudaMalloc(&d_b, shorterPoly.length*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));
	cudaMalloc(&d_temp, len*sizeof(int));
	cudaMemset(d_temp, 0, len*sizeof(int));

	// Do this for all polys in Polyset
	for (int i = 1; i <= NUMPRIMES; i++) {
		/*
		if (i == 5) {
			int k = 0;
		}
		*/
		cudaMemset(d_out, 0, len*sizeof(int));

		// Copy input data from host to device
		cudaMemcpy(d_a, longerPoly.members[i].coeffs, longerPoly.length*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, shorterPoly.members[i].coeffs, shorterPoly.length*sizeof(int), cudaMemcpyHostToDevice);
  		
		for (int j = 0; j < shorterPoly.length; j++) {
			cudaMemset(d_temp, 0, len*sizeof(int));

			// Launch kernel to compute and store modded polynomial values
			monomialScalarMultMods<<<(len + TPB - 1) / TPB, TPB>>>(d_temp, d_a, shorterPoly.members[i].coeffs[j], j, primes[i], longerPoly.length);
			addMods << <(len + TPB - 1) / TPB, TPB >> >(d_out, d_temp, d_out, primes[i], longerPoly.length);
		}
		// Copy results from device to host
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(int), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
    cudaFree(d_temp);
}

int2 LCM(int2 a, int2 b){
	// by default, return the product of the two numbers and the larger of the two powers
	int2 c = {a.x*b.x, a.y >= b.y ? a.y : b.y};
	return c;
}

void sPoly(Poly a, Poly b, Poly c, int* primes)
{
	int len = c.length; // should be 1 less than longerPoly.length
	
	// Declare pointers to device arrays
	int *d_a = 0;
	int *d_b = 0;
	int *d_out = 0;
	int *d_tempA = 0;
	int *d_tempB = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_a, a.length*sizeof(int));
	cudaMalloc(&d_b, b.length*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));
	cudaMalloc(&d_tempA, len*sizeof(int));
	cudaMemset(d_tempA, 0, len*sizeof(int));
	cudaMalloc(&d_tempB, len*sizeof(int));
	cudaMemset(d_tempB, 0, len*sizeof(int));
	
	// Do this for all polys in Polyset
	for (int i = 0; i < NUMPRIMES; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_a, a.members[i].coeffs, a.length*sizeof(int), 
				   cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, b.length*sizeof(int), 
				   cudaMemcpyHostToDevice);
  		
		// find LCM of highest power monomial in a and b.
		int2 lastA = {a.members[i].coeffs[a.length-1], a.length-1};
		int2 lastB = {b.members[i].coeffs[b.length-1], b.length-1};
		int2 lcm = LCM(lastA, lastB);
		
		
		int aScalar = lcm.x / lastA.x;
		int aMonomial = lcm.y - lastA.y;
		int bScalar = lcm.x / lastB.x;
		int bMonomial = lcm.y = lastB.y;
		
		// multiply both a and b by the requisite scale factor to get the last term to equal the LSM
		//monomialScalarMultMods<<<(len + TPB - 1)/TPB, TPB>>>(d_tempA, d_a, aScalar, aMonomial, primes[i], a.length);
		//monomialScalarMultMods<<<(len + TPB - 1)/TPB, TPB>>>(d_tempB, d_b, bScalar, bMonomial, primes[i]), b.length;
				
		// return a-b 
		// does this work, allowing negative mods??
		scalarMultMods<<<(len + TPB - 1)/TPB, TPB>>>(d_tempB, d_tempB, -1, primes[i], a.length); 
		addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_tempA, d_tempB, d_out, primes[i], len);

		// Copy results from device to host
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(int), 
				   cudaMemcpyDeviceToHost);
	}
  
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
        cudaFree(d_tempA);
	cudaFree(d_tempB);
}
