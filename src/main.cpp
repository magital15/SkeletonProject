#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

bool testPrime = false;
bool testPolyDense = false;
bool testGetModsA = false;
bool testGetModsB = false;
bool testAddMods = false;
bool testScalarMods = false;
bool testSubtractMods = false;
bool testMultInverse = false;
bool testMultiplyMods = false;
bool testSPoly = false;
bool testExpPoly = true;
bool testExpPolyGPU = true;


// Setting Primes with new Initializers.cu
int primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3}; 
int* primeArray = setPrimes(primes);

int largePrimes[] = {
	64951, 64969, 64997,  65003,  65011,  65027,  65029,  65033,  65053,  65063,
	65071, 65089, 65099,  65101,  65111,  65119,  65123,  65129,  65141,  65147,
	65167, 65171, 65173,  65179,  65183,  65203,  65213,  65239,  65257,  65267,
	65269, 65287, 65293,  65309,  65323,  65327,  65353,  65357,  65371,  65381,
	65393, 65407, 65413,  65419,  65423,  65437,  65447,  65449,  65479,  65497,
	65519, 65521, 65537 };
int* largePrimeArray = setPrimes(largePrimes);

int main() 
{
// ################################################################ //
// 				Notation:	1 + 2x + 3x^2 = [1, 2, 3]				//
//																	//
//								Poly 1								//
//																	//
	int Polynomial1[] = { 1, 1 };
	int p1len = sizeof(Polynomial1)/sizeof(*Polynomial1);
	Poly a = makeNewPoly(Polynomial1, p1len, primeArray);

//																	//
//								Poly 2								//
//																	//
	int Polynomial2[] = { -301, -302, 303, 11, 12, 15, 71, 234, -1000 };
	int p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b = makeNewPoly(Polynomial2, p2len, primeArray);
//																	//
// ################################################################ //
//																	//
//						Still Need to Implement						//
//																	//

//	a * b = f				--TESTING--
//	a / b = g + (remainder)	--TESTING--
//  sPoly of a, b = h		--TESTING--

//																	//
//##################################################################//
//																	//
//						   Testing Section							//

    if (testPrime == true)
	{
		for (int i = 0; i < NUMPRIMES; i++)
		{
			printf("Prime Array [%i] is: %i\n", i, primeArray[i]);
		}
	}

	if (testPolyDense == true)
	{
		for (int i = 0; i < a.length; i++)
		{
			printf("Dense Poly 1 in a[%i] is: %i\n", i, 
					a.members[0].coeffs[i]);
		}
	}
	
	if (testGetModsA == true)
	{
		printf("Poly a:\n");
		printForReconstruction(a, primeArray);
	}

	if (testGetModsB == true)
	{
		printf("Poly b:\n");
		printForReconstruction(b, primeArray);
	}

	if (testAddMods == true)
	{
		//  a + b = c	--WORKS--
		Poly c = makeAddPoly(a, b);
		addPolys(a, b, c, primeArray);
//		addPolys(a, b, a, primeArray);	// ALSO WORKS

		printf("Poly c:\n");
		printForReconstruction(c, primeArray);
	}

	if (testScalarMods == true)
	{
		//  a * scalar = d	--WORKS--
		Poly d = makeScalarPoly(a);
		scalarMultPoly(a, d, -2, primeArray);
//		scalarMultPoly(a, a, -2, primeArray);	// ALSO WORKS

		printf("Poly d:\n");
		printForReconstruction(d, primeArray);
	}

	if (testSubtractMods == true)
	{
		//  a - b = e	--WORKS--
		Poly e = makeAddPoly(a, b);
		subtractPolys(a, b, e, primeArray);
//		subtractPolys(a, b, a, primeArray);		// PROBABLY WORKS
			
		printf("Poly e:\n");
		printForReconstruction(e, primeArray);	
	}


	if (testMultiplyMods == true)
	{
		//  a * b = f	--TESTING--
		Poly f = makeMultiplyPoly(a, b);
		multiplyPolys(a, b, f, primeArray);

		printf("Poly f:\n");
		printForReconstruction(f, primeArray);
	}

	if (testSPoly == true)
	{
		//  sPoly of a, b = g   --TESTING--
		Poly g = makeSPoly(a, b);
		sPoly(a, b, g, primeArray);
//		sPoly(a, b, a, primeArray);		// Should work the same
	
		printf("Letting A-B be: ");
		for (int i = 0; i < a.length-1; i++)
		{
			printf("%ix^%i + ", a.members[0].coeffs[i], i);
		}
		printf("%ix^%i \n", a.members[0].coeffs[a.length-1], a.length-1);
		printf("Letting B be: ");
		for (int i = 0; i < b.length-1; i++)
		{
			printf("%ix^%i + ", b.members[0].coeffs[i], i);
		}
		printf("%ix^%i \n", b.members[0].coeffs[b.length - 1], b.length - 1);

		printForReconstruction(g, primeArray);
		
	}

	clock_t kernelBegin;
	clock_t kernelEnd;
	double kernelTime;
	int numExp = 4;

	if (testExpPoly == true)
	{
		kernelBegin = clock();
		Poly h = exponentiate(a, numExp, largePrimeArray);
		kernelEnd = clock();

		kernelTime = ((double)(kernelEnd - kernelBegin)) / CLOCKS_PER_SEC;
		printf("Exponention with frequent memcpy's: %f\n", kernelTime*1000.f);

		printForReconstruction(h, largePrimeArray);
	}

	if (testExpPolyGPU == true)
	{
		kernelBegin = clock();
		Poly i = exponentiate_stayOnGPUUntilEnd(a, numExp, largePrimeArray);
		kernelEnd = clock();
		
		kernelTime = ((double)(kernelEnd - kernelBegin)) / CLOCKS_PER_SEC;
		printf("Stay on GPU exponentiation: %f time\n", kernelTime*1000.f);
		
		printForReconstruction(i, largePrimeArray);
	}

//																	//
//##################################################################//

	return 0;
}
