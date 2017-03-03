#ifndef KERNEL_H
#define KERNEL_H
#define NUMPRIMES 16 // Can only be as long as the array in Prime Setter

typedef struct PolyDense_t {
	double base;		// base
	double* coeffs;	// coefficients
}PolyDense;

typedef struct Poly_t {
	double length;		// length
	PolyDense members[NUMPRIMES + 1];
}Poly;

// Kernel wrapper for functions
void getMods(Poly p, double* primes);
void addPolys(Poly a, Poly b, Poly c, double* primes);
void scalarMultPoly(Poly a, Poly c, double scalar, double* primes);
void subtractPolys(Poly a, Poly b, Poly c, double* primes);
void multiplyPolys(Poly a, Poly b, Poly c, double* primes);

#endif
