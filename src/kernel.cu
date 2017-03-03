#include "kernel.h"
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define TPB 32

// Commented out because this is currently working in the main.
/*
#pragma region gcd / chinese remainder thm helper methods, reconstruct kernel and launcher


// C function for extended Euclidean Algorithm
// code authors: GeeksForGeeks
__device__
int gcdExtended(int a, int b, int *x, int *y)
{
	// Base Case
	if (a == 0)
	{
		*x = 0, *y = 1;
		return b;
	}

	int x1, y1; // To store results of recursive call
	int gcd = gcdExtended(b%a, a, &x1, &y1);

	// Update x and y using results of recursive
	// call
	*x = y1 - (b / a) * x1;
	*y = x1;

	return gcd;
}

// Function to find modulo inverse of a
// PRECONDITION: a and m are coprime
// code authors: GeeksForGeeks
__device__
int2 modInverse(int a, int m)
{
	int x, y;
	int g = gcdExtended(a, m, &x, &y);
	if (g != 1) {
		printf("Inverse doesn't exist for %d, %d", a, m);
		return int2{ 0, 0 };
	}
	else
	{
		// m is added to handle negative x
		int res = (x%m + m) % m;
		int otherRes = (y%a + a) % a;
		

		//printf("%d * %d = 1 (mod %d)\n", a, res, m);
		//printf("So %d is the multiplicative inverse of %d (mod %d)\n", res, a, m);
		//printf("%d * %d = 1 (mod %d)\n", m, otherRes, a);
		//printf("So %d is the multiplicative inverse of %d (mod %d)\n", otherRes, m, a);

		return int2{ res, otherRes };
	}
}

__device__
int reconstruct(Poly a, int col, int *primeArray) {

	int nextMember = 1;
	int nextPrime = 0;

	int a1 = a.members[nextMember++].coeffs[col];
	int a2 = a.members[nextMember].coeffs[col];

	int p1 = primeArray[nextPrime++];
	int p2 = primeArray[nextPrime];

	int a12 = 0;
	int prevAnswer = -1;

	while (nextPrime < NUMPRIMES && a12 != prevAnswer) {

		// only enter this on iterations beyond the first one
		if (nextPrime != 1) {
			a1 = a12;
			a2 = a.members[nextMember].coeffs[col];
			p1 = p1*p2;
			p2 = primeArray[nextPrime];
		}

		int2 multiplicativeInverses = modInverse(p1, p2);
		int k1modp2 = multiplicativeInverses.x;
		int k2modp1 = multiplicativeInverses.y;

		prevAnswer = a12;
		a12 = (a2*k1modp2*p1 + a1*k2modp1*p2) % (p1*p2);

		nextMember++;
		nextPrime++;
	}

	return a12;
}

void reconstructKernel(Poly a, int *primes, int size)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= size) return;
	a.members[0].coeffs[i] = reconstruct(a, i, primes); // POINTER ERROR
}

void reconstructPoly(Poly in, int* primes)
{
	int len = in.length;

	// Declare pointers to device arrays
	int *d_in = 0;
	int *d_primes = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(int));
	cudaMalloc(&d_primes, NUMPRIMES*sizeof(int));

	cudaMemcpy(d_primes, primes, NUMPRIMES*sizeof(int),cudaMemcpyHostToDevice);
	
	// HOW TO ACTUALLY MOVE ALL DATA OVER?

	// Do this for all polys in Polyset
	for (int i = 1; i < NUMPRIMES + 1; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[i].coeffs, len*sizeof(int),
			cudaMemcpyHostToDevice);
	}

	//paralelize by coefficient number
	for (int i = 0; i < len; i++) {
		// Launch kernel to compute and store modded polynomial values
		reconstructKernel <<<(len + TPB - 1) / TPB, TPB >> >(d_in, primes, NUMPRIMES);
	}
	
	// Copy results from device to host
	cudaMemcpy(in.members[0].coeffs, d_in, len*sizeof(int),	cudaMemcpyDeviceToHost);
	

	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(primes);
}
#pragma endregion 
*/

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
void addMods(int* d_out, int* d_a, int* d_b, int mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_a[i];
    const float y = d_b[i];
	d_out[i] = getRemainder(x + y, mod);
}

__global__
void scalarMultMods(int* d_out, int* d_in, int scalar, int mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getRemainder(x*scalar, mod);
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
	}
	else
	{
		len = b.length;
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
											  primes[i-1]);

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
											scalar, primes[i-1]);

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

void multiplyPolys(Poly a, Poly b, Poly c, int* primes)
{
	int len = c.length;

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
}
