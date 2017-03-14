#include "kernel.h"
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

void doTest(bool testArray[], Poly a, Poly b, int* primeArray, int* d_primes)
{
	// testSetPrimes
	if (testArray[0] == true) {
		printf("--Test setPrimes()--\n");
		for (int i = 0; i < NUMPRIMES; i++) {
			printf("primeArray[%i] is: %i\n", i, primeArray[i]);
		}
		printf("\n");
	}

	// testMakeNewPoly
	if (testArray[1] == true) {
		printf("--Test makeNewPoly()--\n");
		for (int i = 0; i < a.length; i++) {
			printf("a.members[0].coeffs[%i] is: %i\n", i, 
					a.members[0].coeffs[i]);
		}
		printForReconstruction(a, primeArray);
		printf("\n");
		for (int i = 0; i < b.length; i++) {
			printf("b.members[0].coeffs[%i] is: %i\n", i, 
					b.members[0].coeffs[i]);
		}
		printForReconstruction(b, primeArray);
		printf("\n");
	}
	
	// testAddPolys
	if (testArray[2] == true) {
		printf("--Test addPolys()--\n");
		Poly c = addPolys(a, b, primeArray);
		printf("Poly A + B\n");
		printForReconstruction(c, primeArray);
		printf("\n");		
	}
	
	// testScalarMods
	if (testArray[3] == true) {
		int scalar = -2;	// Change the scalar input
		printf("--Test scalarMods()--\n");
		Poly c = scalarMultPoly(a, scalar, primeArray);
		printf("Poly %i*A\n", scalar);
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testSubtractPolys 
	if (testArray[4] == true) {
		printf("--Test subtractPolys()--\n");
		Poly c = subtractPolys(a, b, primeArray);
		printf("Poly A - B\n");
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testMultiplyPolys
	if (testArray[5] == true) {
		printf("--Test multiplyPolys()--\n");
		Poly c = multiplyPolys(a, b, primeArray);
		printf("Poly A * B\n");
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testSPoly
	if (testArray[6] == true) {
		printf("--Test sPoly()--\n");
		Poly c = sPoly(a, b, primeArray);
		printf("sPoly of A, B\n");
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testExpPoly
	if (testArray[7] == true) {
		int exponent = 3;	// Change the exponent input
		printf("--Test exponentiate()--\n");
		Poly c = exponentiate(a, exponent, primeArray);
		printf("Poly A^%i\n", exponent);
		printForReconstruction(c, primeArray);
		printf("\n");
	}
	
	// testExpPolyGPU
	if (testArray[8] == true) {
		int exponent = 3;	// Change the exponent input
		printf("--Test exponentiateGPU()--\n");
		Poly c = exponentiateGPU(a, exponent, primeArray);
		printf("Poly A^%i\n", exponent);
		printForReconstruction(c, primeArray);
		printf("\n");
	}

	// testAddOnGPU
	if (testArray[9] == true) {
		printf("--Test addOnGPU()--\n");
		int* d_a = makeGPUPoly(a);
		int* d_b = makeGPUPoly(b);	
		int* d_c = addGPU(d_a, d_b, d_primes);
		Poly c = getGPUPoly(d_c);
		printForReconstruction(c, primeArray);
		printf("\n");
	}

	// testMultipleGPU
	if (testArray[10] == true) {
		int timesAdd = 32;
		printf("--Test multipleGPU()--\n");
		int* d_a = makeGPUPoly(a);
		int* d_b = makeGPUPoly(b);	
		int* d_c = 0;
		for (int i = 0; i < timesAdd; i++) {
			d_c = addGPU(d_a, d_c, d_primes);
		}
		Poly c = getGPUPoly(d_c);
		printForReconstruction(c, primeArray);
		printf("\n");
	}

	// testAddPolySpeeds
	if (testArray[11] == true) {
		// Setup number of times iterations
		int times[] = {5, 10, 20, 50, 75, 100, 250, 500, 750, 1000};
		int lenTimes = sizeof(times)/sizeof(*times);

		// Initialize
		Poly c = addPolys(a, b, primeArray);
		int* d_a = makeGPUPoly(a);
		int* d_b = makeGPUPoly(b);
		int* d_c = 0;

		float milliseconds1[lenTimes];
		float milliseconds2[lenTimes];
		
		// Test 1
		cudaEvent_t start1, stop1, start2, stop2;
		cudaEventCreate(&start1);
		cudaEventCreate(&stop1);
		cudaEventCreate(&start2);
		cudaEventCreate(&stop2);

		for (int i = 0; i < lenTimes; i++) {
			cudaEventRecord(start1);
			for (int j = 0; j < times[i]; j++) {
				c = addPolys(a, c, primeArray);
			}
			cudaEventRecord(stop1);
			cudaEventSynchronize(stop1);
			milliseconds1[i] = 0;
			cudaEventElapsedTime(&milliseconds1[i], start1, stop1);
		}
		
		for (int i = 0; i < lenTimes; i++) {
			cudaEventRecord(start2);
			for (int j = 0; j < times[i]; j++) {
				d_c = addGPU(d_a, d_c, d_primes);
			}
			cudaEventRecord(stop2);
			cudaEventSynchronize(stop2);
			milliseconds2[i] = 0;
			cudaEventElapsedTime(&milliseconds2[i], start2, stop2);
		}

		// Print Results
		for (int i = 0; i < lenTimes; i++) {
			printf("Times Added = %i\n", times[i]);
			printf("CPU transfer (ms) = %f\n", milliseconds1[i]);
			printf("GPU only (ms) = %f\n", milliseconds2[i]);
		}
	}

	// testScalarModsGPU
	if (testArray[12] == true) {
		printf("--Test scalarModsGPU()--\n");
		int* d_a = makeGPUPoly(a);	
		int* d_c = scalarModsGPU(d_a, -1, d_primes);
		Poly c = getGPUPoly(d_c);
		printForReconstruction(c, primeArray);
		printf("\n");
	}

	// testSubtractGPU
	if (testArray[13] == true) {
		printf("--Test subtractGPU()--\n");
		int* d_a = makeGPUPoly(a);
		int* d_b = makeGPUPoly(b);	
		int* d_c = subtractGPU(d_a, d_b, d_primes);
		Poly c = getGPUPoly(d_c);
		printForReconstruction(c, primeArray);
		printf("\n");
	}

	// testSubtractGPU
	if (testArray[14] == true) {
		printf("--Test multGPU()--\n");
		int* d_a = makeGPUPoly(a);
		int* d_b = makeGPUPoly(b);	
		int* d_c = multGPU(d_a, d_b, d_primes);
		Poly c = getGPUPoly(d_c);
		printForReconstruction(c, primeArray);
		printf("\n");
	}
}
