#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>

// Function that sets the primes
int* setPrimes(int primes[]) {
	int* result = (int*)calloc(NUMPRIMES, sizeof(int));
	for (int i = 0; i < NUMPRIMES; i++) {
		result[i] = primes[i];
	}
	return result;
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
