#ifndef KERNEL_H
#define KERNEL_H
#define NUMPRIMES 5 // Can only be as long as the array in Prime Setter

typedef struct PolyDense_t {
	int base;		// base
	int* coeffs;	// coefficients
}PolyDense;

typedef struct Poly_t {
	int length;		// length
	PolyDense members[NUMPRIMES + 1];
}Poly;

// Kernel wrapper functions
void getMods(Poly p, int* primes);
void addPolys(Poly a, Poly b, Poly c, int* primes);
void scalarMultPoly(Poly a, Poly c, int scalar, int* primes);
void subtractPolys(Poly a, Poly b, Poly c, int* primes);
void multiplyPolys(Poly a, Poly b, Poly c, int* primes);
void sPoly(Poly a, Poly b, Poly c, int* primes);

// Initializer wrapper functions
int* setPrimes(int primes[]);
Poly makeNewPoly(int coeffArray[], int len, int primes[]);
Poly makeAddPoly(Poly a, Poly b);
Poly makeScalarPoly(Poly a);
Poly makeMultiplyPoly(Poly a, Poly b);

#endif
