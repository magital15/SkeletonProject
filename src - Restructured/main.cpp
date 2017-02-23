#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>

bool testPrime = false;
bool testPolyDense = false;
bool testGetMods = false;
bool testAddMods = false;
bool testScalarMods = false;
bool testSubtractMods = true;


void primeSetter(int* primeHolder) {
    int primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2};
	for (int i = 0; i < NUMPRIMES; i++)
	{
		primeHolder[i] = primes[i];
	}
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
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		a.members[i].coeffs = (int*)calloc(p1len, sizeof(int));
		a.members[i].length = p1len;
	}
	for (int i = 0; i < p1len; i++)
	{
		a.members[0].coeffs[i] = Polynomial1[i];
	}
//																	//
//								Poly 2								//
//																	//
	int Polynomial2[] = {11, 22, 33, 44, 55, 66, 77, 88, 99, 110};
	int p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		b.members[i].coeffs = (int*)calloc(p2len, sizeof(int));
		b.members[i].length = p2len;
	}
	for (int i = 0; i < p2len; i++)
	{
		b.members[0].coeffs[i] = Polynomial2[i];
	}
//																	//
//																	//
// ################################################################ //
//																	//
//	 					 Prime Holding Array					    //
//																	//
	int* primeArray = (int*)calloc(NUMPRIMES, sizeof(int));
	primeSetter(primeArray);
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
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		c.members[i].coeffs = (int*)calloc(a.members[0].length, 
										   sizeof(int));
		c.members[i].length = a.members[0].length;
	}
	addPolys(a, b, c, primeArray);
//	addPolys(a, b, a, primeArray);	// ALSO WORKS


//  a * scalar = d	--WORKS--

	Poly d;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		d.members[i].coeffs = (int*)calloc(a.members[0].length, 
										   sizeof(int));
		d.members[i].length = a.members[0].length;
	}
	scalarMultPoly(a, d, -2, primeArray);
//	scalarMultPoly(a, a, -2, primeArray);	// ALSO WORKS

//  a - b = e	--TESTING--

	Poly e;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		e.members[i].coeffs = (int*)calloc(a.members[0].length, 
										   sizeof(int));
		e.members[i].length = a.members[0].length;
	}
	subtractPolys(a, b, e, primeArray);

//  a * b = f	--TESTING--




//  a / b = g + (remainder)	--TESTING--



//																	//
//##################################################################//
//																	//																	//
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
		for (int i = 0; i < a.members[0].length; i++)
		{
			printf("Dense Poly 1 in a[%i] is: %i\n", i, 
					a.members[0].coeffs[i]);
		}
	}
	
	if (testGetMods == true)
	{
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < a.members[0].length; i++)
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
				for (int i = 0; i < c.members[0].length; i++)
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
				for (int i = 0; i < d.members[0].length; i++)
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
				for (int i = 0; i < e.members[0].length; i++)
				{
					printf("e.members[%i].c[%i] is: %i\n", j, i, 
							e.members[j].coeffs[i]);
				}
			}	
		}

//																	//
//##################################################################//

	return 0;
}
