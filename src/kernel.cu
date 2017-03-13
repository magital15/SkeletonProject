#include "kernel.h"
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define TPB 32

// DEVICE WORK FUNCTIONS
__device__
int getRemainder(int coefficient, int mod)
{
	coefficient = coefficient % mod;
	coefficient += mod;
	return coefficient % mod;
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
}

__global__
void scalarMultMods(int* d_out, int* d_in, int scalar, int mod, int size)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= size) return;
	const float x = d_in[i];
	d_out[i] = getRemainder(x*scalar, mod);
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
}

// KERNEL WRAPPER FUNCTIONS
void getMods(Poly in, int* primes) {
	int len = in.length;

	// Declare pointers to device arrays
	int *d_in = 0;
	int *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));

	// Take Mod of each member
	for (int i = 1; i <= NUMPRIMES; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[0].coeffs, len*sizeof(int), cudaMemcpyHostToDevice);

		// Launch kernel to compute and store modded polynomial values
		takeMod<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in, primes[i-1]);

		// Copy results from device to host
		cudaMemcpy(in.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
	}
	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(d_out);
}

Poly addPolys(Poly a, Poly b, int* primes) {
	// Deal with lengths
	int len = a.length >= b.length ? a.length : b.length;
	a = a.length < b.length ? copyIntoBigger(a, b.length) : a;
	b = a.length > b.length ? copyIntoBigger(b, a.length) : b;
	Poly result = makePolyGivenLength(len);

	// Declare pointers to device arrays
	int *d_a = 0;
	int *d_b = 0;
	int *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_a, len*sizeof(int));
	cudaMalloc(&d_b, len*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));

	// Do this for all members
	for (int i = 1; i <= NUMPRIMES; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_a, a.members[i].coeffs, len*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, len*sizeof(int), cudaMemcpyHostToDevice);

		// Launch kernel to compute and store modded polynomial values
		addMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a, d_b, primes[i-1], len);

		// Copy results from device to host
		cudaMemcpy(result.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
	}
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);

	return result;
}

Poly scalarMultPoly(Poly in, int scalar, int* primes) {
	int len = in.length;
	Poly result = makePolyGivenLength(len);

	// Declare pointers to device arrays
	int *d_in = 0;
	int *d_out = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_in, len*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));

	// Do this for all members
	for (int i = 1; i <= NUMPRIMES; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_in, in.members[i].coeffs, len*sizeof(int), cudaMemcpyHostToDevice);

		// Launch kernel to compute and store modded polynomial values
		scalarMultMods<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_in, scalar, primes[i-1], len);

		// Copy results from device to host
		cudaMemcpy(result.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
	}
	// Free the memory allocated for device arrays
	cudaFree(d_in);
	cudaFree(d_out);

	return result;
}

// Can be accelerated by not copying in between two arithmetic operators
Poly subtractPolys(Poly a, Poly b, int* primes) {
	Poly result = scalarMultPoly(b, -1, primes);
	result = addPolys(a, result, primes);
	return result;
}

Poly multiplyPolys(Poly a, Poly b, int* primes) {
	// Deal with lengths
	int len = a.length + b.length - 1;
	Poly result = makePolyGivenLength(len);
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

	// Do this for all members
	for (int i = 1; i <= NUMPRIMES; i++) {
		// Copy input data from host to device
		cudaMemcpy(d_a, longerPoly.members[i].coeffs, longerPoly.length*sizeof(int), cudaMemcpyHostToDevice);
  		cudaMemset(d_out, 0, len*sizeof(int));

		for (int j = 0; j < shorterPoly.length; j++) {
			cudaMemset(d_temp, 0, len*sizeof(int));
			// Launch kernel to compute and store modded polynomial values
			monomialScalarMultMods<<<(len + TPB - 1) / TPB, TPB>>>(d_temp, d_a, shorterPoly.members[i].coeffs[j], j, primes[i-1], longerPoly.length);
			addMods<<<(len + TPB - 1) / TPB, TPB>>>(d_out, d_temp, d_out, primes[i-1], len);
		}
		// Copy results from device to host
		cudaMemcpy(result.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
	}
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
    cudaFree(d_temp);
	
	return result;
}

int2 LCM(int2 a, int2 b) {
	// Return the product of the two numbers and the larger of the two powers
	int2 c = {a.x*b.x, a.y >= b.y ? a.y : b.y};
	return c;
}

Poly sPoly(Poly a, Poly b, int* primes) {
	int len = a.length >= b.length ? a.length-1 : b.length-1;
	Poly result = makePolyGivenLength(len);

	// Declare pointers to device arrays
	int *d_a = 0;
	int *d_b = 0;
	int *d_out = 0;
	int *d_tempA = 0;
	int *d_tempB = 0;

	//Allocate memory for device arrays
	cudaMalloc(&d_a, (len+1)*sizeof(int));
	cudaMalloc(&d_b, (len+1)*sizeof(int));
	cudaMalloc(&d_out, len*sizeof(int));
	cudaMalloc(&d_tempA, (len+1)*sizeof(int));
	cudaMalloc(&d_tempB, (len+1)*sizeof(int));
	
	// Do this for all members
	for (int i = 1; i <= NUMPRIMES; i++) {
		cudaMemset(d_a, 0, (len+1)*sizeof(int));
		cudaMemset(d_b, 0, (len+1)*sizeof(int));
		cudaMemset(d_out, 0, len*sizeof(int));
		cudaMemset(d_tempA, 0, (len+1)*sizeof(int));
		cudaMemset(d_tempB, 0, (len+1)*sizeof(int));

		// Copy input data from host to device
		cudaMemcpy(d_a, a.members[i].coeffs, a.length*sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_b, b.members[i].coeffs, b.length*sizeof(int), cudaMemcpyHostToDevice);
  		
		// find LCM of highest power monomial in a and b.
		int2 lastA = {a.members[i].coeffs[a.length-1], a.length-1};
		int2 lastB = {b.members[i].coeffs[b.length-1], b.length-1};
		int2 lcm = LCM(lastA, lastB);

		int aScalar = lcm.x / lastA.x;
		int aMonomial = lcm.y - lastA.y;
		int bScalar = lcm.x / lastB.x;
		int bMonomial = lcm.y - lastB.y;
		int currPrime = primes[i-1];
		// Multiply a and b by the requisite scale factor to get the last term to equal the LSM
		monomialScalarMultMods<<<(a.length + TPB - 1)/TPB, TPB>>>(d_tempA, d_a, aScalar, aMonomial, currPrime, len+1);
		monomialScalarMultMods<<<(b.length + TPB - 1) / TPB, TPB>>>(d_tempB, d_b, bScalar, bMonomial, currPrime, len+1);
						
		// Return a - b 
		scalarMultMods<<<(len + TPB - 1) / TPB, TPB >>>(d_tempB, d_tempB, -1, currPrime, len);
		addMods<<<(len + TPB - 1) / TPB, TPB>>>(d_out, d_tempB, d_tempA, currPrime, len);

		// Copy results from device to host
		cudaMemcpy(result.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
	}
  	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_out);
    cudaFree(d_tempA);
	cudaFree(d_tempB);
	
	return result;
}

Poly exponentiate(Poly a, int exp, int* primeArray) {
	Poly result = copyIntoBigger(a, a.length);
	for(int i = 1; i < exp; i++) {
		result = multiplyPolys(a, result, primeArray);
	}
	return result;
}

Poly exponentiateGPU(Poly a, int exp, int* primes) {
	int len = a.length + (exp-1)*(a.length-1);
	Poly result = makePolyGivenLength(len);
	
    // Declare pointers to device arrays
    int *d_a = 0;
    int *d_out = 0;
    int *d_temp = 0;

	//Allocate memory for device arrays
    cudaMalloc(&d_a, len*sizeof(int));
    cudaMalloc(&d_out, len*sizeof(int));
    cudaMalloc(&d_temp, len*sizeof(int));

	// Do this for all members
	for (int i = 1; i <= NUMPRIMES; i++) {
		// Reset memory and set lengths
		cudaMemset(d_out, 0, len*sizeof(int));
		cudaMemset(d_a, 0, len*sizeof(int));
		cudaMemcpy(d_a, a.members[i].coeffs, a.length*sizeof(int), cudaMemcpyHostToDevice);	
		int currentLen = a.length*2 - 1;
		int otherLen = currentLen - a.length + 1;

		for (int numExp = 1; numExp < exp; numExp++) {
			for (int j = 0; j < a.length; j++) {
				// Launch kernel to compute and store modded polynomial values			
				cudaMemset(d_temp, 0, currentLen*sizeof(int));
				monomialScalarMultMods<<<(len + TPB - 1) / TPB, TPB>>>(d_temp, d_a, a.members[i].coeffs[j], j, primes[i-1], otherLen);
				addMods <<<(len + TPB - 1) / TPB, TPB >>>(d_out, d_temp, d_out, primes[i-1], currentLen);
			}
			// Copy d_out into d_a, then reset d_out
			cudaMemcpy(d_a, d_out, len*sizeof(int), cudaMemcpyDeviceToDevice);
			if (numExp != exp-1) {
				cudaMemset(d_out, 0, len*sizeof(int));
			}
			// Update length variables
			otherLen = currentLen;
			currentLen += a.length - 1;
		}
		// Copy results from device to host
		cudaMemcpy(result.members[i].coeffs, d_out, len*sizeof(int), cudaMemcpyDeviceToHost);
	}
	// Free the memory allocated for device arrays
	cudaFree(d_a);
	cudaFree(d_out);
    cudaFree(d_temp);
	
	return result;
}



//  ###########################################################
// I can't believe this works
// ##
// To do:
// 1. Move GPU poly functions to a different .cu file
// 2. Make input poly length dynamic
// 3. Make arithmetic functions return GPU polys

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

// Gets a primesArray on GPU
int* makeGPUPrimes(int* primes) {
	// Create Pointers to arrays on CPU and GPU
	int* GPUresult = 0;

	// Set memory allocation for CPU and GPU results
	cudaMalloc(&GPUresult, NUMPRIMES*sizeof(int));

	// Copy the primes into the GPUresult
	cudaMemcpy(GPUresult, primes, NUMPRIMES*sizeof(int), cudaMemcpyHostToDevice);

	return GPUresult;
}

Poly getGPUPoly(int* d_a) {
	// Grab the length from the device poly
	int len = lenGrabber(d_a);
	int d_len = len*NUMPRIMES + 1;
	
	// Make a new poly based on the length
	Poly result = makePolyGivenLength(len);
	
	// Copy d_a into a CPUresult
	// CONSIDER DOING THIS ABOVE INSTEAD OF JUST LENGRABBER
	int* CPUresult = (int*)calloc(d_len, sizeof(int));
	cudaMemcpy(CPUresult, d_a, d_len*sizeof(int), cudaMemcpyDeviceToHost);

	int pos = 1;
	for (int i = 1; i <= NUMPRIMES; i++) {
		for (int j = 0; j < len; j++) {
			result.members[i].coeffs[j] = CPUresult[pos];
			pos++;
		}
	}

	cudaFree(d_a);
	return result;
}

// WORKS
__global__
void addMods2(int* d_out, int* d_a, int* d_b, int* primes, int len)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i == 0)	{
		d_out[i] = len;
		return;
	}
	if (i > len*NUMPRIMES) return; // Good
	const int r = (i - 1) / len;
	const int mod = primes[r];
	const int x = d_a[i];
	const int y = d_b[i];
	d_out[i] = getRemainder(x + y, mod);
}

// WORKS
void addGPU(int* d_out, int* d_a, int* d_b, int* d_primes) {
	// Grab the length
	int lenA = lenGrabber(d_a);
	int lenB = lenGrabber(d_b);
	int len = lenA >= lenB ? lenA : lenB;
	int d_len = len*NUMPRIMES + 1;
	// Copy smaller GPU Poly into bigger if  needed
	// CAN BE ACCELERATED BY NOT ALWAYS COPYING SMALLER GPU POLY
	// OR CAN BE ACCELERATED BY ADJUSTING ADDMODS FUNCTION (INPUT BOTH LENGTHS)
	int* d_a2 = 0;
	int* d_b2 = 0;
	cudaMalloc(&d_a2, d_len*sizeof(int));
	cudaMalloc(&d_b2, d_len*sizeof(int));
	d_a2 = copyGPUGivenLen(d_a, len);
	d_b2 = copyGPUGivenLen(d_b, len);

	// Do the kernel call
	addMods2<<<(len + TPB - 1)/TPB, TPB>>>(d_out, d_a2, d_b2, d_primes, len);

	cudaFree(d_a2);
	cudaFree(d_b2);
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
// LOTS OF COPYS IN THIS ONE: CAN BE ACCELERATED, SEE BELOW FOR EXAMPLE
int* copyGPUGivenLen(int* d_in, int len) {
	Poly in = getGPUPoly(d_in);
	in = copyIntoBigger(in, len);
	int* result = makeGPUPoly(in);
	return result;
}

/* EXAMPLE OF ACCELERATED GPUCOPY
// Copys a GPU poly into a bigger one
int* copyGPUGivenLen(int* d_in, int len) {
	// Get Length of d_in
	int origLen = lenGrabber(d_in);
	printf("Orig Len = %i\n", origLen);
	int d_origLen = origLen*NUMPRIMES + 1;
	printf("d_origLen = %i\n", d_origLen);
	int d_len = len*NUMPRIMES + 1;

	// Allocate memory for the result
	int* result = 0;
	cudaMalloc(&result, d_len*sizeof(int));
	// Allocate memory for the Host holder
	int* CPUholder = (int*)calloc(d_origLen, sizeof(int));
	// Allocate memory for new Host
	int* CPUresult = (int*)calloc(d_len, sizeof(int));
	// Copy Device data into CPU data
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
	// THIS COPY BREAKS IT
	cudaMemcpy(result, CPUresult, d_len*sizeof(int), cudaMemcpyHostToDevice);

	free(CPUholder);
	free(CPUresult);
	return result;
}
*/
















