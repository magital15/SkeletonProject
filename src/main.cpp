#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Still Need to Implement
// a / b = g + (remainder)

// To add another test: Create the bool, add it to the end of the
// testHolder[], then add the testing procedure in tests.cu in a
// format similar to the others
bool testSetPrimes = false;		//	Sets the Prime Array
bool testMakeNewPoly = false;	//	Constructs Original Polynomials
bool testAddPolys = false;		//	a + b = c			--WORKS--
bool testScalarMods = false;	//	a * scalar = c		--WORKS--
bool testSubtractPolys = false;	//	a - b = c			--WORKS--
bool testMultiplyPolys = false;	//	a * b = c			--WORKS--
bool testSPoly = false;			//	sPoly of a, b = c   --WORKS--
bool testExpPoly = false;		//	a ^ exp = c			--WORKS--
bool testExpPolyGPU = false;	//	a ^ exp = c (fast)	--WORKS--
bool testAddOnGPU = false;		//	a + c = d			--WORKS--
bool testMultipleGPU = false;	//	a + c + c... = d	--WORKS--
bool testAddPolySpeeds = false;	//	Speed test			--WORKS--
bool testScalarModsGPU = false;	//	a * scalar = d		--WORKS--
bool testSubtractGPU = false;	//	a - b = d			--WORKS--
bool testMultGPU = false;		//	a * b = d			--TESTING--
bool testHolder[] = {testSetPrimes,testMakeNewPoly,testAddPolys,
testScalarMods,testSubtractPolys,testMultiplyPolys,testSPoly,
testExpPoly,testExpPolyGPU,testAddOnGPU,testMultipleGPU,testAddPolySpeeds,
testScalarModsGPU,testSubtractGPU,testMultGPU};

// Set the prime and device prime array
int Primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3}; 
int* primeArray = setPrimes(Primes);
int* d_primes = makeGPUPrimes(primeArray);

int main() {
// Notation: 1 + 2x + 3x^2 = [1, 2, 3]
// Poly A
	int Polynomial1[] = { 1, 1, 1, 1 };
	int p1len = sizeof(Polynomial1)/sizeof(*Polynomial1);
	Poly a = makeNewPoly(Polynomial1, p1len, primeArray);

// Poly B
	int Polynomial2[] = { 1, 1, 1, 1 };
	int p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b = makeNewPoly(Polynomial2, p2len, primeArray);

// Testing Section
	doTest(testHolder, a, b, primeArray, d_primes);

	return 0;
}
