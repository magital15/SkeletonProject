#ifndef KERNEL_H
#define KERNEL_H
#define NUMPRIMES 5 // Can only be as long as the array in Prime Setter

typedef struct PolyDense_t {
	int length;		// length
	int base;		// base
	int* coeffs;	// coefficients
}PolyDense;

typedef struct Poly_t {
	PolyDense members[NUMPRIMES + 1];
}Poly;

// Kernel wrapper for functions
void getMods(Poly p, int* primes);
void addPolys(Poly a, Poly b, Poly c, int* primes);
void scalarMultPoly(Poly a, Poly c, int scalar, int* primes);
void subtractPolys(Poly a, Poly b, Poly c, int* primeArray);

#endif
