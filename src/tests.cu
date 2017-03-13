#include "kernel.h"
#include <stdio.h>

void doTest(bool testArray[], Poly a, Poly b, int* primeArray)
{
	// testSetPrimes
	if (testArray[0] == true)
	{
		printf("--Test setPrimes()--\n");
		for (int i = 0; i < NUMPRIMES; i++)
		{
			printf("primeArray[%i] is: %i\n", i, primeArray[i]);
		}
		printf("\n");
	}

	// testMakeNewPoly
	if (testArray[1] == true)
	{
		printf("--Test makeNewPoly()--\n");
		for (int i = 0; i < a.length; i++)
		{
			printf("a.members[0].coeffs[%i] is: %i\n", i, 
					a.members[0].coeffs[i]);
		}
		printForReconstruction(a, primeArray);
		printf("\n");
		for (int i = 0; i < b.length; i++)
		{
			printf("b.members[0].coeffs[%i] is: %i\n", i, 
					b.members[0].coeffs[i]);
		}
		printForReconstruction(b, primeArray);
		printf("\n");
	}
	
	// testAddPolys
	if (testArray[2] == true)
	{
		printf("--Test addPolys()--\n");
		Poly c = addPolys(a, b, primeArray);
		printf("Poly A + B\n");
		printForReconstruction(c, primeArray);
		printf("\n");		
	}
	
	// testScalarMods
	if (testArray[3] == true)
	{
		int scalar = -2;	// Change the scalar input
		printf("--Test scalarMods()--\n");
		Poly c = scalarMultPoly(a, scalar, primeArray);
		printf("Poly %i*A\n", scalar);
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testSubtractPolys 
	if (testArray[4] == true)
	{
		printf("--Test subtractPolys()--\n");
		Poly c = subtractPolys(a, b, primeArray);
		printf("Poly A - B\n");
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testMultiplyPolys
	if (testArray[5] == true)
	{
		printf("--Test multiplyPolys()--\n");
		Poly c = multiplyPolys(a, b, primeArray);
		printf("Poly A * B\n");
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testSPoly
	if (testArray[6] == true)
	{
		printf("--Test sPoly()--\n");
		Poly c = sPoly(a, b, primeArray);
		printf("sPoly of A, B\n");
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testExpPoly
	if (testArray[7] == true)
	{
		int exponent = 3;	// Change the exponent input
		printf("--Test exponentiate()--\n");
		Poly c = exponentiate(a, exponent, primeArray);
		printf("Poly A^%i\n", exponent);
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testExpPolyGPU
	if (testArray[8] == true)
	{
		int exponent = 3;	// Change the exponent input
		printf("--Test exponentiateGPU()--\n");
		Poly c = exponentiateGPU(a, exponent, primeArray);
		printf("Poly A^%i\n", exponent);
		printForReconstruction(c, primeArray);
		printf("\n");
	}
}
