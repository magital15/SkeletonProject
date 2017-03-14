#ifndef KERNEL_H
#define KERNEL_H
#define NUMPRIMES 5 // Can only be as long as the array in Primes on Main.cpp

typedef struct PolyDense_t {
	int* coeffs;
}PolyDense;

typedef struct Poly_t {
	int length;
	PolyDense members[NUMPRIMES + 1];
}Poly;

// Test function
void doTest(bool array[], Poly a, Poly b, int* primeArray, int* d_primes);

// Kernel wrapper functions
void getMods(Poly p, int* primes);
Poly addPolys(Poly a, Poly b, int* primes);
Poly scalarMultPoly(Poly a, int scalar, int* primes);
Poly subtractPolys(Poly a, Poly b, int* primes);
Poly multiplyPolys(Poly a, Poly b, int* primes);
Poly sPoly(Poly a, Poly b, int* primes);
Poly exponentiate(Poly a, int exp, int* primes);
Poly exponentiateGPU(Poly a, int exp, int* primes);

// Initializer wrapper functions
int* setPrimes(int primes[]);
Poly makeNewPoly(int coeffArray[], int len, int primes[]);
Poly makePolyGivenLength(int length);
Poly copyIntoBigger(Poly a, int len);
void printForReconstruction(Poly g, int* primeArray);

// Initializer wrapper functions (GPU)
int* makeGPUPrimes(int* primes);
int* makeGPUPoly(Poly a);
Poly getGPUPoly(int* d_in);
int* copyGPUGivenLen(int* d_in, int len);
int lenGrabber(int* d_in);

// GPU wrapper functions
int* addGPU(int* d_a, int* d_b, int* d_primes);
int* scalarModsGPU(int* d_in, int scalar, int* d_primes);
int* subtractGPU(int* d_a, int* d_b, int* d_primes);
int* multGPU(int* d_a, int* d_b, int* d_primes);
int* appendGPU(int* d_a, int* d_b);


#endif
