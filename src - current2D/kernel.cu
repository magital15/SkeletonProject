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
void addMods(int* d_out, int* d_a, int* d_b, int* primes, int size)
{
	const int c = blockIdx.x*blockDim.x + threadIdx.x;
	const int r = blockIdx.y*blockDim.y + threadIdx.y;
//	if (i > size) return;
	int mod = primes[r];
	const float x = d_a.members[r].coeffs[c];
	const float y = d_b.members[r].coeffs[c];
	d_c.members[r].coeffs[c] = getRemainder(x + y, mod);
	printf("Hello?");
}

__global__
void scalarMultMods(int* d_out, int* d_in, int scalar, int mod, int size)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= size) return;
	const float x = d_in[i];
	d_out[i] = getRemainder(x*scalar, mod);

	//printf("negative [%d] = %d \n", i, d_out[i]);
}

// expects d_out to be large enough to hold i+monomial elements
// not sure if works for negative scalars
__global__
void monomialScalarMultMods(int* d_out, int* d_in, int scalar, int monomial, int mod, int size)
{
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
		if (i >= size) return;
        const float x = d_in[i];
        d_out[i + monomial] = getRemainder(x*scalar, mod);
		
		//printf("scaled [%d] = %d \n", i, d_out[i]);
}

__global__
void minusMult(int* d_out, int* d_in, int scalar, int mod)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const float x = d_in[i];
	d_out[i] = getMult(x, scalar);


	//printf("negative [%d] = %d \n", i, d_out[i]);
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
	cudaMalloc(&d_a, len*NUMPRIMES*sizeof(int));
	cudaMalloc(&d_b, len*NUMPRIMES*sizeof(int));
	cudaMalloc(&d_out, len*NUMPRIMES*sizeof(int));

	int pos = 0;
	// Copy input data from host to device
	for (int i = 1; i < NUMPRIMES+1; i++) {
		cudaMemcpy(d_a, a.members[i].coeffs, len*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, len*sizeof(int), cudaMemcpyHostToDevice);
		pos += len;
	}

	// Launch kernel to compute and store modded polynomial values
	addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b, 
										  primes, len);
	
	// Copy input data from host to device
	for (int i = 1; i < NUMPRIMES+1; i++) {
		cudaMemcpy(c.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
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
	
	shorterPoly = copyIntoBigger(shorterPoly, len);
	longerPoly = copyIntoBigger(longerPoly, len);
	
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
		cudaMemset(d_out, 0, len*sizeof(int));

		// Copy input data from host to device
		cudaMemcpy(d_a, longerPoly.members[i].coeffs, longerPoly.length*sizeof(int), cudaMemcpyHostToDevice);

		// unless we paralelize and run all the shorterPoly monomial multiplications at once, we aren't  currently
		// making use of having d_b over on the device
		// cudaMemcpy(d_b, shorterPoly.members[i].coeffs, shorterPoly.length*sizeof(int), cudaMemcpyHostToDevice);
  		
		for (int j = 0; j < shorterPoly.length; j++) {
			cudaMemset(d_temp, 0, len*sizeof(int));

			// Launch kernel to compute and store modded polynomial values
			monomialScalarMultMods<<<(len + TPB - 1) / TPB, TPB>>>(d_temp, d_a, shorterPoly.members[i].coeffs[j], j, primes[i-1], longerPoly.length);
			addMods<<<(len + TPB - 1) / TPB, TPB>>>(d_out, d_temp, d_out, primes[i-1], len);
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
	cudaMalloc(&d_a, (len+1)*sizeof(int));
	cudaMalloc(&d_b, (len + 1)*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));
	cudaMalloc(&d_tempA, (len + 1)*sizeof(int));
	cudaMalloc(&d_tempB, (len + 1)*sizeof(int));
	
	// Do this for all polys in Polyset
	for (int i = 1; i <= NUMPRIMES; i++) {
		cudaMemset(d_a, 0, (len + 1)*sizeof(int));
		cudaMemset(d_b, 0, (len + 1)*sizeof(int));
		cudaMemset(d_out, 0, len*sizeof(int));
		cudaMemset(d_tempA, 0, (len + 1)*sizeof(int));
		cudaMemset(d_tempB, 0, (len + 1)*sizeof(int));

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
		int bMonomial = lcm.y - lastB.y;
		/*
		printf("\n-------------------\n");
		printf("A has highest number %dx^%d and B has highest number %dx^%d \n", lastA.x, lastA.y, lastB.x, lastB.y);
		printf("So their LCM is %dx^%d \n", lcm.x, lcm.y);
		
		printf("\nlenth of a %d\n", a.length);
		printf("lenth of b %d\n\n", b.length);


		printf("aScalar: %d\n", aScalar);
		printf("aMonomial: %d\n", aMonomial);
		printf("bScalar: %d\n", bScalar);
		printf("bMonomial: %d\n", bMonomial);
		*/

		int currPrime = primes[i - 1];

		// multiply both a and b by the requisite scale factor to get the last term to equal the LSM
		monomialScalarMultMods<<<(a.length + TPB - 1)/TPB, TPB>>>(d_tempA, d_a, aScalar, aMonomial, currPrime, len+1);
		monomialScalarMultMods << <(b.length + TPB - 1) / TPB, TPB >> >(d_tempB, d_b, bScalar, bMonomial, currPrime, len + 1);
						
		// return a-b 
		scalarMultMods << <(len + TPB - 1) / TPB, TPB >> >(d_tempB, d_tempB, -1, currPrime, len);
		addMods << <(len + TPB - 1) / TPB, TPB >> >(d_out, d_tempB, d_tempA, currPrime, len);

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

Poly exponentiate(Poly a, int exp, int* primeArray) {
	int len = exp*(a.length - 1) + 1;
	
	Poly result = copyIntoBigger(a, len);
	// ideally we won't have to move the intermediate results back to the CPU at each step

	for(int i = 1; i < exp; i++) {
		multiplyPolys(result, a, result, primeArray);
	}
	return result;
}
