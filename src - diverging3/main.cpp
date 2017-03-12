#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>

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
bool testExpPoly = false;
bool testMemory = false;
bool testMemoryAdd = false;
bool testMovePoly2D = true;


// Setting Primes with new Initializers.cu
int primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3}; 
int* primeArray = setPrimes(primes);

int main() 
{
// ################################################################ //
// 				Notation:	1 + 2x + 3x^2 = [1, 2, 3]				//
//																	//
//								Poly 1								//
//																	//
	int Polynomial1[] = { 1, 1 };
	int p1len = sizeof(Polynomial1)/sizeof(*Polynomial1);
	Poly a = makeOriginalPoly(Polynomial1, p1len, primeArray);

//																	//
//								Poly 2								//
//																	//
	int Polynomial2[] = { 2, 3 };
	int p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b = makeOriginalPoly(Polynomial2, p2len, primeArray);
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
	
	if (testExpPoly == true)
	{
		Poly h = exponentiate(a, 3, primeArray);
		printf("Poly h:\n");
		printForReconstruction(h, primeArray);
	}
/*
	if (testMemoryAdd == true)
	{
		int Polynomial3[] = { 32, 33, 34, 35};
		int p3len = sizeof(Polynomial3)/sizeof(*Polynomial3);
		Poly a = makeOriginalPoly(Polynomial1, p1len, primeArray);
		Poly mem = makeOriginalPoly(Polynomial3, p3len, primeArray);
		int* d_mem = 0;
		moveToDevice(mem, d_mem);
		moveFromDevice(mem, d_mem);
		printf("Poly mem after moving to and from device:\n");
		printForReconstruction(mem, primeArray);
	}
*/
	if (testMovePoly2D == true)
	{
		PolyDevice d_a;
		PolyDevice d_b;
		moveToDevice2D(a, d_a);
		moveFromDevice2D(a, d_a);
		moveToDevice2D(b, d_b);
		Poly u = makeNewPoly();
		PolyDevice d_u;
		moveToDevice2D(u, d_u);
		addPolysTest2(d_u, d_a, d_b, primeArray);
//		moveFromDevice2D(u, d_u, primeArray);
//		moveFromDevice2D(a, d_a, primeArray);
//		printForReconstruction(a, primeArray);
		printForReconstruction(u, primeArray);
//		basicFunction(d_u);
	}


//																	//
//##################################################################//

	return 0;
}
