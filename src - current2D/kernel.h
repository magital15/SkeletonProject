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

typedef struct PolyDenseDevice_t {
	int* coeffs;
}PolyDenseDevice;

typedef struct PolyDevice_t {
	int length;	//IMPORTANTE
	PolyDenseDevice members[NUMPRIMES];
}PolyDevice;

// Kernel wrapper functions
void getMods(Poly p, int* primes);
void addPolys(Poly a, Poly b, Poly c, int* primes);
void scalarMultPoly(Poly a, Poly c, int scalar, int* primes);
void subtractPolys(Poly a, Poly b, Poly c, int* primes);
void multiplyPolys(Poly a, Poly b, Poly c, int* primes);
void sPoly(Poly a, Poly b, Poly c, int* primes);
Poly exponentiate(Poly a, int exp, int* primes);
void moveToDevice(Poly in, int* d_in);
void moveFromDevice(Poly out, int d_out);


	// Testing
void moveToDevice2D(Poly in, PolyDevice d_in);
void moveFromDevice2D(Poly out, PolyDevice d_out);
void addPolysTest(Poly a, Poly b, int* primes);
void addPolysTest2(PolyDevice c, PolyDevice a, PolyDevice b, int* primes);
PolyDevice copyIntoBiggerDevice(PolyDevice a, int len);
void basicFunction(PolyDevice in);
void printForReconstruction2(Poly in, int* primeArray);

// Initializer wrapper functions
int* setPrimes(int primes[]);
Poly makeOriginalPoly(int coeffArray[], int len, int primes[]);
Poly makePolyGivenLength(int length);
Poly makeAddPoly(Poly a, Poly b);
Poly makeScalarPoly(Poly a);
Poly makeMultiplyPoly(Poly a, Poly b);
Poly makeSPoly(Poly a, Poly b);
Poly copyIntoBigger(Poly a, int len);
Poly makeNewPoly();
void printForReconstruction(Poly g, int* primeArray);




#endif
