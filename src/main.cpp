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
bool testSPoly = true;


// Setting Primes with new Initializers.cu
int primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3}; 
int* primeArray = setPrimes(primes);

void printForReconstruction(Poly g) {
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

Poly exponentiate(Poly a, int exp) {
	int len = (a.length-1)*exp;
	
	Poly resEven = copyIntoBigger(a, len);
	Poly resOdd = makePolyGivenLength(len);
	
	// ideally we won't have to move the intermediate results back to the CPU at each step
	for(int i = 1; i <= exp; i++) {
		if (i % 2 == 0) {
			multiplyPolys(resOdd, a, resEven, primeArray); 
		else {
			multiplyPolys(resEven, a, resOdd, primeArray); 
		}
	}
	
	if (exp % 2 == 0) {
		return resEven;
	} else {
		return resOdd;
	}
}

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
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < a.length; i++)
			{
				printf("a.members[%i].c[%i] is: %i\n", j, i, 
						a.members[j].coeffs[i]);
			}
		}	
	}

	if (testGetModsB == true)
	{
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < b.length; i++)
			{
				printf("b.members[%i].c[%i] is: %i\n", j, i, 
						b.members[j].coeffs[i]);
			}
		}	
	}

	if (testAddMods == true)
	{
		//  a + b = c	--WORKS--

		Poly c = makeAddPoly(a, b);
		addPolys(a, b, c, primeArray);
//		addPolys(a, b, a, primeArray);	// ALSO WORKS

		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < c.length; i++)
			{
				printf("c.members[%i].c[%i] is: %i\n", j, i, 
						c.members[j].coeffs[i]);
			}
		}	
	}

	if (testScalarMods == true)
	{
		//  a * scalar = d	--WORKS--

		Poly d = makeScalarPoly(a);
		scalarMultPoly(a, d, -2, primeArray);
//		scalarMultPoly(a, a, -2, primeArray);	// ALSO WORKS

		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < d.length; i++)
			{
				printf("d.members[%i].c[%i] is: %i\n", j, i, 
						d.members[j].coeffs[i]);
			}
		}	
	}

	if (testSubtractMods == true)
	{
		//  a - b = e	--WORKS--

		Poly e = makeAddPoly(a, b);
		subtractPolys(a, b, e, primeArray);
//		subtractPolys(a, b, a, primeArray);		// PROBABLY WORKS
			
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < e.length; i++)
			{
				printf("e.members[%i].c[%i] is: %i\n", j, i, 
						e.members[j].coeffs[i]);
			}
		}	
	}


	if (testMultiplyMods == true)
	{
		//  a * b = f	--TESTING--
		Poly f = makeMultiplyPoly(a, b);
		multiplyPolys(a, b, f, primeArray);

		for (int j = 0; j <= NUMPRIMES; j++)
		{
			printf("p%d f.members[%i] is: ", primeArray[j-1], j);
			for (int i = 0; i < f.length; i++)
			{
				printf("%i,", f.members[j].coeffs[i]);
			}
			printf("\n");
		}	
	}

	if (testSPoly == true)
	{
		//  sPoly of a, b = g   --TESTING--

		Poly g = makeSPoly(a, b);
		sPoly(a, b, g, primeArray);

		// should work the same
		// sPoly(a, b, a, primeArray);
		
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

		printForReconstruction(g);
		
	}

//																	//
//##################################################################//

	return 0;
}
