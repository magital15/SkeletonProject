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

Poly makePolyGivenLength(int length) {
	// Create a Poly to store data into
	Poly result;

	result.length = length;

	// Allocate memory based on the length found
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(result.length, sizeof(int));
	}

	return result;
}

Poly makeScalarPoly(Poly a)
{
	int length = a.length;
	return makePolyGivenLength(length);
}

// This function initializes a polynomial that will be added or subtracted
// initializes the memory for a poly of length max(length of a, length of b)
Poly makeAddPoly(Poly a, Poly b)
{
	int length = a.length >= b.length ? a.length : b.length;
	return makePolyGivenLength(length);
}

// initializes the memory for a poly of length one less than the product of lengths of a, b
Poly makeMultiplyPoly(Poly a, Poly b)
{
	int length = a.length + b.length - 1;
	return makePolyGivenLength(length);
}

// initializes the memory for a poly with length one less than the max length of a, b
Poly makeSPoly(Poly a, Poly b)
{
	int length = a.length >= b.length ? a.length - 1 : b.length - 1;
	return makePolyGivenLength(length);
}

Poly copyIntoBigger(Poly a, int len)
{
	// Create a Poly to store data into
	Poly result;

	// Give the result the required length
	result.length  = len;

	// Copy the original data into the new poly
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		result.members[i].coeffs = (int*)calloc(result.length, sizeof(int));
		for (int j = 0; j < a.length; j++)
		{
			result.members[i].coeffs[j] = a.members[i].coeffs[j];
		}
	}

	return result;
}
