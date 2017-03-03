#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

bool testPrime = false;
bool testPolyDense = false;
bool testGetMods = false;
bool testAddMods = false;
bool testScalarMods = false;
bool testSubtractMods = false;
bool testMultInverse = false;
bool testMultiplyMods = false;
bool testReconstruct = true;


void primeSetter(double* primeHolder) {
	//double primes[] = {31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2};
    //double primes[] = {29989, 29983, 29959, 29947, 29927, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2};
	//double primes[] = {610331, 610327, 610301, 610289, 610279};
	double primes[] = {739, 733, 727, 719, 709, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2};
	//double primes[] = {11, 7, 5, 3, 2};
	for (int i = 0; i < NUMPRIMES; i++)
	{
		primeHolder[i] = primes[i];
	}
}

double gcdExtended(double a, double b, int *x, int *y) {

	// Base Case
	if (a == 0)
	{
		*x = 0, *y = 1;
		return b;
	}

	int x1, y1; // To store results of recursive call
	
	double gcd = gcdExtended(double(int(fmod(b,a))), a, &x1, &y1);

	// Update x and y using results of recursive call
	*x = y1 - int((b / a)) * x1;
	*y = x1;
	//printf("*x is: %i\n*y is: %i\n", y1 - int((b / a)) * x1, x1);
	return gcd;
}

void multiplicativeInverse(double* primes, double *k1, double *k2) {
	
	double* prime1 = (double*)calloc(NUMPRIMES - 1, sizeof(double));
	double* prime2 = (double*)calloc(NUMPRIMES - 1, sizeof(double));

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
	 	double g = gcdExtended(prime1[i], prime2[i], &x, &y);
		if (g != 1) {
			k1[i] = 0;
			k2[i] = 0;
		}
		else
		{
			double res = double(int(fmod(x,prime2[i])));
			double otherRes = double(int(fmod(y, prime1[i])));
			k1[i] = res;
			k2[i] = otherRes;
		}
		//printf("FOR I IS: %i\nx is: %f\ny is: %f\ng is: %f\n",i,x,y,g);
	}
	free(prime1);
	free(prime2);
}

double reconstruct(double *coeffColumn, double *primes, double *k1, double *k2) {
	double a1 = coeffColumn[0];
	double a2 = coeffColumn[1];
	double p1 = primes[0];
	double p2 = primes[1];
	double a12 = 0;
	double prevAnswer = -1;
	
	int counter = 0;
	while (counter < NUMPRIMES - 1 && a12 != prevAnswer) {
		if (counter != 0) {
			a1 = a12;
			a2 = coeffColumn[counter + 1];
			p1 = p1*p2;
			p2 = primes[counter + 1];
		}

		prevAnswer = a12;
		a12 = fmod((a2*k1[counter]*p1 + a1*k2[counter]*p2), (p1*p2));
		/*
		if (p1*p2 < 2147483647) {
			printf("Failure at counter = %i\np1*p2 is %f\n",counter, p1*p2);
		}
		*/
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
	double Polynomial1[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	double p1len = sizeof(Polynomial1)/sizeof(*Polynomial1);
	Poly a;
	a.length = p1len;
//																	//
//								Poly 2								//
//																	//
	double Polynomial2[] = {11, 22, 33, 44, 55, 66, 77, 88, 99, 110};
	double p2len = sizeof(Polynomial2)/sizeof(*Polynomial2);
	Poly b;
	b.length = p2len;
//																	//
//					Poly Memory Allocate & Value Set				//
//																	//
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		a.members[i].coeffs = (double*)calloc(p1len + p2len, sizeof(double));
	}
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		b.members[i].coeffs = (double*)calloc(p2len + p1len, sizeof(double));
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
	double* primeArray = (double*)calloc(NUMPRIMES, sizeof(double));
	primeSetter(primeArray);

	double* k1modp2 = (double*)calloc(NUMPRIMES - 1, sizeof(double));
	double* k2modp1 = (double*)calloc(NUMPRIMES - 1, sizeof(double));
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
		c.members[i].coeffs = (double*)calloc(c.length, sizeof(double));
	}
	addPolys(a, b, c, primeArray);
//	addPolys(a, b, a, primeArray);	// ALSO WORKS


//  a * scalar = d	--WORKS--

	Poly d;
	d.length = a.length;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		d.members[i].coeffs = (double*)calloc(a.length, sizeof(double));
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
		e.members[i].coeffs = (double*)calloc(e.length, sizeof(double));
	}
	subtractPolys(a, b, e, primeArray);
//	subtractPolys(a, b, a, primeArray);		// PROBABLY WORKS

//  a * b = f	--TESTING--

	Poly f;
	f.length = a.length + b.length - 1;
	for (int i = 0; i < NUMPRIMES + 1; i++)
	{
		f.members[i].coeffs = (double*)calloc(f.length, sizeof(double));
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
			printf("Prime Array [%i] is: %f\n", i, primeArray[i]);
		}
	}

	if (testPolyDense == true)
	{
		for (int i = 0; i < a.length; i++)
		{
			printf("Dense Poly 1 in a[%i] is: %f\n", i, 
					a.members[0].coeffs[i]);
		}
	}
	
	if (testGetMods == true)
	{
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < a.length; i++)
			{
				printf("a.members[%i].c[%i] is: %f\n", j, i, 
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
					printf("c.members[%i].c[%i] is: %f\n", j, i, 
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
					printf("d.members[%i].c[%i] is: %f\n", j, i, 
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
					printf("e.members[%i].c[%i] is: %f\n", j, i, 
							e.members[j].coeffs[i]);
				}
			}	
		}

	if (testMultInverse == true)
	{
		for (int i = 0; i < NUMPRIMES - 1; i++)
		{
			printf("k1modp2[%i] is: %f\n", i, k1modp2[i]);
			printf("k2modp1[%i] is: %f\n", i, k2modp1[i]);
		}
	}

	if (testMultiplyMods == true)
	{
		for (int j = 0; j < NUMPRIMES + 1; j++)
		{
			for (int i = 0; i < f.length; i++)
			{
				printf("f.members[%i].c[%i] is: %f\n", j, i, 
						f.members[j].coeffs[i]);
			}
		}	
	}

	if (testReconstruct == true)
	{
		int timesWrong = 0;
		for (int j = 1000; j < 15000; j++) {
			double* coeffColumn = (double*)calloc(NUMPRIMES, sizeof(double));
			double realAnswer = j;
			for (int i = 0; i < NUMPRIMES; i++)
			{
				coeffColumn[i] = double(int(fmod(realAnswer, primeArray[i])));
				//printf("coeffColumn[%i] is: %f\n", i, coeffColumn[i]);
			}		
			double answer = reconstruct(coeffColumn, primeArray, k1modp2, k2modp1);
			
			if (answer == realAnswer) {
				if (timesWrong > 0) {
					printf("Times Wrong = %i\n", timesWrong);
					timesWrong = 0;
				}
				printf("This answer was found: %f\n", answer);
			}
			else {
				timesWrong++;
			}
		}
	}
//																	//
//##################################################################//

	return 0;
}
