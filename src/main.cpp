#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

bool testPrime = false;
bool testPolyDense = false;
bool testGetMods = false;
bool testAddMods = false;
bool testScalarMods = false;
bool testSubtractMods = false;
bool testMultInverse = true;
bool testMultiplyMods = false;
bool testReconstruct = true;


void primeSetter(int* primeHolder) {
    int primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2};
	for (int i = 0; i < NUMPRIMES; i++)
	{
		primeHolder[i] = primes[i];
	}
}

int gcdExtended(int a, int b, int *x, int *y) {

	// Base Case
	if (a == 0)
	{
		*x = 0, *y = 1;
		return b;
	}

	int x1, y1; // To store results of recursive call
	
	int gcd = gcdExtended(b%a, a, &x1, &y1);

	// Update x and y using results of recursive call
	*x = y1 - (b / a) * x1;
	*y = x1;

	return gcd;
}

void multiplicativeInverse(int* primes, int *k1, int *k2) {
	
	int* prime1 = (int*)calloc(NUMPRIMES - 1, sizeof(int));
	int* prime2 = (int*)calloc(NUMPRIMES - 1, sizeof(int));

	prime1[0] = primes[0];
	prime2[0] = primes[1];
	
	for (int i = 0; i < NUMPRIMES - 1; i++)
	{
		if (i != 0)
		{
			prime1[i] = prime1[i-1] * primes[i];
			prime2[i] = primes[i+1];
		}

		int x, y;
	 	int g = gcdExtended(prime1[i], prime2[i], &x, &y);
		if (g != 1) {
			k1[i] = 0;
			k2[i] = 0;
		}
		else
		{
			int res = (x % prime2[i] + prime2[i]) % prime2[i];
			int otherRes = (y % prime1[i] + prime1[i]) % prime1[i];
			k1[i] = res;
			k2[i] = otherRes;
		}
	}
	free(prime1);
	free(prime2);
}

int reconstruct(int *coeffColumn, int *primes, int *k1, int *k2) {
	int a1 = coeffColumn[0];
	int a2 = coeffColumn[1];
	int p1 = primes[0];
	int p2 = primes[1];
	int a12 = 0;
	int prevAnswer = -1;
	
	int counter = 0;
	while (counter < NUMPRIMES - 1 && a12 != prevAnswer) {
		if (counter != 0) {
			a1 = a12;
			a2 = coeffColumn[counter + 1];
			p1 = p1*p2;
			p2 = primes[counter + 1];
		}

		prevAnswer = a12;
		a12 = (a2*k1[counter]*p1 + a1*k2[counter]*p2) % (p1*p2);

		counter++;
	}

	return a12;
}

int main()
{
// ################################################################ //
// 				Notation:	1 + 2x + 3x^2 = [1, 2, 3]				//
//																	//
//								Poly 1								//
//																	//
	int Polynomial1[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	int p1len = sizeof(Polynomial1)/sizeof(*Polynomial1);
	Poly a;
	a.length = p1len;
//																	//
//								Poly 2								//
//																	//
	int Polynomial2[] = {11, 22, 33, 44, 55, 66, 77, 88, 99, 110};
	int p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b;
	b.length = p2len;
//																	//
//					Poly Memory Allocate & Value Set				//
//																	//
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		a.members[i].coeffs = (int*)calloc(p1len + p2len, sizeof(int));
	}
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		b.members[i].coeffs = (int*)calloc(p2len + p1len, sizeof(int));
	}
	for (int i = 0; i < p1len; i++)
	{
		a.members[0].coeffs[i] = Polynomial1[i];
	}
	for (int i = 0; i < p2len; i++)
	{
		b.members[0].coeffs[i] = Polynomial2[i];
	}
// ################################################################ //
//																	//
//	 						 Prime Affairs						    //
//																	//
	int* primeArray = (int*)calloc(NUMPRIMES, sizeof(int));
	primeSetter(primeArray);

	int* k1modp2 = (int*)calloc(NUMPRIMES - 1, sizeof(int));
	int* k2modp1 = (int*)calloc(NUMPRIMES - 1, sizeof(int));
	multiplicativeInverse(primeArray, k1modp2, k2modp1);
	
//																	//
//																	//
//					Getting Mods of Original Polys  	            //
//																	//
	getMods(a, primeArray);
	getMods(b, primeArray);
//																	//
// ################################################################ //
//																	//
//						  Arithmetic Section						//
//																	//
//  a + b = c	--WORKS--
	
	Poly c;
	if (a.length > b.length)
	{
		c.length = a.length;
	}
	else
	{
		c.length = b.length;
	}
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		c.members[i].coeffs = (int*)calloc(c.length, sizeof(int));
	}
	addPolys(a, b, c, primeArray);
//	addPolys(a, b, a, primeArray);	// ALSO WORKS


//  a * scalar = d	--WORKS--

	Poly d;
	d.length = a.length;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		d.members[i].coeffs = (int*)calloc(a.length, sizeof(int));
	}
	scalarMultPoly(a, d, -2, primeArray);
//	scalarMultPoly(a, a, -2, primeArray);	// ALSO WORKS

//  a - b = e	--TESTING--

	Poly e;
	if (a.length > b.length)
	{
		e.length = a.length;
	}
	else
	{
		e.length = b.length;
	}
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		e.members[i].coeffs = (int*)calloc(e.length, sizeof(int));
	}
	subtractPolys(a, b, e, primeArray);
//	subtractPolys(a, b, a, primeArray);		// PROBABLY WORKS

//  a * b = f	--TESTING--

	Poly f;
	f.length = a.length + b.length - 1;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		f.members[i].coeffs = (int*)calloc(f.length, sizeof(int));
	}
	multiplyPolys(a, b, f, primeArray);
//	multiplyPolys(a, b, a, primeArray);		// DOES NOT WORK
											// BECAUSE OF SIZE


//  a / b = g + (remainder)	--TESTING--



//																	//
//##################################################################//
//																	//																	//
//						Reconstruction Section						//





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
	
	if (testGetMods == true)
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

	if (testAddMods == true)
		{
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
			for (int j = 0; j < NUMPRIMES + 1; j++)
			{
				for (int i = 0; i < e.length; i++)
				{
					printf("e.members[%i].c[%i] is: %i\n", j, i, 
							e.members[j].coeffs[i]);
				}
			}	
		}

	if (testMultInverse == true)
	{
		for (int i = 0; i < NUMPRIMES - 1; i++)
		{
			printf("k1modp2[%i] is: %i\n", i, k1modp2[i]);
			printf("k2modp1[%i] is: %i\n", i, k2modp1[i]);
		}
	}

	if (testMultiplyMods == true)
	{
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < f.length; i++)
			{
				printf("f.members[%i].c[%i] is: %i\n", j, i, 
						f.members[j].coeffs[i]);
			}
		}	
	}

	if (testReconstruct == true)
	{
		int* coeffColumn = (int*)calloc(NUMPRIMES, sizeof(int));
		int realAnswer = 6000;
		for (int i = 0; i < NUMPRIMES; i++)
		{
			coeffColumn[i] = (realAnswer%primeArray[i] + primeArray[i]) % primeArray[i];
			printf("coeffColumn[%i] is: %i\n", i, coeffColumn[i]);
		}		
		int answer = reconstruct(coeffColumn, primeArray, k1modp2, k2modp1);
		printf("Answer for reconstruct is: %i\n", answer);
	}
//																	//
//##################################################################//

	return 0;
}
