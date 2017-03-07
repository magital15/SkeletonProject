#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#define SAFETYLENGTH

// Function that sets primes to reduce how much code is needed
int* setPrimes(int primes[])
{
	int* result = (int*)calloc(NUMPRIMES, sizeof(int));
	for (int i = 0; i < NUMPRIMES; i++)
	{
		result[i] = primes[i];
	}
	return result;
}

// This function makes original polynomials that will be used for arithmetic
Poly makeNewPoly(int coeffArray[], int len, int primes[])
{
	// Create a Poly to store data into
	Poly result;
	
	// Give it the length as an input
	result.length = len;

	// Allocate memory based on the length
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(result.length, sizeof(int));
	}

	// Copy the coeffArray[] into the new Poly in members[0]
	for (int i = 0; i < result.length; i++)
	{
		result.members[0].coeffs[i] = coeffArray[i];
	}

	// Get the modular results to fill the rest of the members
	getMods(result, primes);

	return result;
}

// This function initializes a polynomial that will be added or subtracted
Poly makeAddPoly(Poly a, Poly b)
{
	// Create a Poly to store data into
	Poly result;
	
	// Choose the right length for adding and subtracting
	result.length = a.length <= b.length ? a.length : b.length;
	result.length = a.length <= b.length ? b.length : a.length;	

	// Allocate memory based on the length found
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(result.length, sizeof(int));
	}

	return result;
}

Poly makeScalarPoly(Poly a)
{
	// Create a Poly to store data into
	Poly result;

	// Choose the same length as the input polynomial
	result.length = a.length;

	// Allocate memory based on the length found
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(result.length, sizeof(int));
	}

	return result;
}

Poly makeMultiplyPoly(Poly a, Poly b)
{
	// Create a Poly to store data into
	Poly result;

	// Choose the right length for adding and subtracting
	result.length = a.length + b.length - 1;

	// Allocate memory based on the length found
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(result.length, sizeof(int));
	}

	return result;
}
