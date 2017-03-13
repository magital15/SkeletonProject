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
bool testSetPrimes = true;		//	Sets the Prime Array
bool testMakeNewPoly = true;	//	Constructs Original Polynomials
bool testAddPolys = true;		//	a + b = c			--WORKS--
bool testScalarMods = true;	//	a * scalar = c		--WORKS--
bool testSubtractPolys = true;	//	a - b = c			--WORKS--
bool testMultiplyPolys = true;	//	a * b = c			--WORKS--
bool testSPoly = true;			//	sPoly of a, b = c   --WORKS--
bool testExpPoly = true;		//	a ^ exp = c			--WORKS--
bool testExpPolyGPU = true;	//	a ^ exp = c (fast)	--WORKS--
bool testWorkOnGPU = true;		//	[a+b=c] a + c = d	--KINDA WORKS--
bool testHolder[] = {testSetPrimes,testMakeNewPoly,testAddPolys,
testScalarMods,testSubtractPolys,testMultiplyPolys,testSPoly,
testExpPoly,testExpPolyGPU,testWorkOnGPU};

// Set the primes array
int Primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3}; 
int* primeArray = setPrimes(Primes);

int main() {
// Notation: 1 + 2x + 3x^2 = [1, 2, 3]
// Poly A
	int Polynomial1[] = { 221, 21512, 141, 44125 };
	int p1len = sizeof(Polynomial1)/sizeof(*Polynomial1);
	Poly a = makeNewPoly(Polynomial1, p1len, primeArray);

// Poly B
	int Polynomial2[] = { 1, 3, 1, 3, 5 };
	int p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b = makeNewPoly(Polynomial2, p2len, primeArray);

// Testing Section
	doTest(testHolder, a, b, primeArray);

	return 0;
}
