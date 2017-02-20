#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

// TURN THESE CONSTANTS ON AND OFF TO WORK ON DIFFERENT GOALS
const bool testing = false;
const bool multipleCoeff = false;
const bool convertPolToNewBaseAndBack = true;

int2 modInverse(int a, int m);
int gcdExtended(int a, int b, int *x, int *y);
int reconstruct(int *shortenedRes, int *mods);

// Function to find modulo inverse of a
// PRECONDITION: a and m are coprime
// code authors: GeeksForGeeks
int2 modInverse(int a, int m)
{
	int x, y;
 	int g = gcdExtended(a, m, &x, &y);
	if (g != 1) {
		printf("Inverse doesn't exist for %d, %d", a, m);
		return int2{ 0, 0 };
	}
	else
	{
		// m is added to handle negative x
		int res = (x%m + m) % m;
		printf("%d * %d = 1 (mod %d)\n", a, res, m);
 		printf("So %d is the multiplicative inverse of %d (mod %d)\n", res, a, m);

		int otherRes = (y%a + a) % a;
		printf("%d * %d = 1 (mod %d)\n", m, otherRes, a);
		printf("So %d is the multiplicative inverse of %d (mod %d)\n", otherRes, m, a);

		return int2{ res, otherRes };
	}
}

// C function for extended Euclidean Algorithm
// code authors: GeeksForGeeks
int gcdExtended(int a, int b, int *x, int *y)
{
	// Base Case
	if (a == 0)
	{
		*x = 0, *y = 1;
		return b;
	}

	int x1, y1; // To store results of recursive call
	int gcd = gcdExtended(b%a, a, &x1, &y1);

	// Update x and y using results of recursive
	// call
	*x = y1 - (b / a) * x1;
	*y = x1;

	return gcd;
}

// PRECONDITION: shortenedRes and mods both contain at least 2 elements
// calculation stops once result stops changing
// TBD negative results
int reconstruct(int *shortenedRes, int *mods, int arrSize) {
	int a1 = shortenedRes[0];
	int a2 = shortenedRes[1];
	int p1 = mods[0];
	int p2 = mods[1];
	int a12 = 0;
	int prevAnswer = -1;
	
	int counter = 1;
	while (counter < numMods && a12 != prevAnswer) {
		if (counter != 1) {
			a1 = a12;
			a2 = shortenedRes[counter];
			p1 = p1*p2;
			p2 = mods[counter];
		}

		int2 multiplicativeInverses = modInverse(p1, p2);
		int k1modp2 = multiplicativeInverses.x;
		int k2modp1 = multiplicativeInverses.y;

		prevAnswer = a12;
		a12 = (a2*k1modp2*p1 + a1*k2modp1*p2) % (p1*p2);
		printf("clue %d gives the number %d \n \n", counter, a12);
		counter++;
	}

	return a12;
}

typedef struct polDense_t {
	int base;
	int length;
	int arr[0];
}*polDense;

typedef struct poly_t {
	polDense members[6];
}*poly;


polDense init(int mod, int *input, int input_length) {
	polDense p = (polDense)malloc(sizeof(int)*(input_length + 2));
	//polDense p = init(input_length);
	(*p).base = mod;
	(*p).length = input_length;
	memcpy(input, (*p).arr, input_length*sizeof(int));
	return p;
}
polDense copy(polDense p) {
	return init((*p).base, (*p).arr, (*p).length);
}

// returns previous mod base
int modAllCoeff(polDense p, int newMod) {
	int temp = (*p).base;
	for (int i = 0; i < (*p).length; i++) {
		(*p).arr[i] %= newMod;
	}
	(*p).base = newMod;
	return temp;
}

polDense copyMod(polDense p, int newMod) {
	polDense res = init((*p).base, (*p).arr, (*p).length);

	//printf("I copied ");
	//printArray((*p).arr, (*p).length);

	modAllCoeff(p, newMod);

	//printf("I modded by %d", newMod);
	//printArray((*res).arr, (*p).length);

	return res;
}

poly init(int *input, int input_length) {
	poly p = (poly)malloc(sizeof(int)*(6 * (input_length + 2)));

	(*p).members[0] = init(0, input, input_length);
	(*p).members[1] = copyMod((*p).members[0], prime0);
	(*p).members[2] = copyMod((*p).members[0], prime1);
	(*p).members[3] = copyMod((*p).members[0], prime2);
	(*p).members[4] = copyMod((*p).members[0], prime3);
	(*p).members[5] = copyMod((*p).members[0], prime4);

	return p;
}

void printArray(int *arr, int size) {
	std::string build = std::to_string(arr[0]);
	for (int i = 0; i < size; i++) {
		build += ", " + std::to_string(arr[i]);
	}
	std::cout << "print array: " << build << "\n";
}

void printArray(polDense p) {
	printArray((*p).arr, (*p).length);
}

// Driver Program
int main()
{
	if (multipleCoeff) {
		// stores result of polynomial division
		int outputCoeff[numCoeff];
		int shortenedRes[numCoeff][numMods]; // get these from cumodp
		int mods[numMods]; // choose these ourselves

		// MAKE THIS PARALLEL
		for (int i = 0; i < numCoeff; i++) {
			outputCoeff[i] = reconstruct(shortenedRes[i], mods, numMods);
		}
	}
	else if (testing) {
		int mods[] = { 3, 5, 7, 11, 13, 17, 19 };
		int clues[numMods];
		for (int i = 0; i < numMods; i++) {
			clues[i] = -4672 % mods[i];
		}
		// SEND THIS TO THE GPU - AT SOME POINT VIA MERGE SORT
		int result = reconstruct(clues, mods, numMods);
	}
	else if (convertPolToNewBaseAndBack) {
		int input[numCoeff];
		for (int i = 0; i < numCoeff; i++) {
			input[i] = 10 * i;
		}
		printArray(input, numCoeff);

		poly p = init(input, numCoeff);
		
		/*for (int i = 0; i < (*(*p).members[0]).length; i++)
		{
			if (i == 0) printf("-Original polDense in poly \n");
			printf("Coefficient %i is: %i\n", i, (*(*p).members[0]).arr[i]);
		}*/
		printArray((*p).members[0]);
		printArray((*p).members[1]);
		printArray((*p).members[2]);
		printArray((*p).members[3]);
		printArray((*p).members[4]);
		printArray((*p).members[5]);

		
		/*
		polDense pol(0,input,numCoeff);
		//printArray(pol.arr, numCoeff);
		
		polDense polModA(0, dummy, numCoeff);
		polModA.copyDataFrom(pol);
		polModA.modAllCoeff(17);
		printArray(polModA.arr, numCoeff);

		polDense polModB(0, dummy, numCoeff);
		polModB.copyDataFrom(pol);
		polModB.modAllCoeff(19);
		printArray(polModB.arr, numCoeff);

		// TEST RECONSTRUCTED ARRAY - WOAH CRASHED
		int temp[numCoeff];
		for (int i = 0; i < numCoeff; i++) {
			int nums[] = { polModA.arr[i], polModB.arr[i] };
			int mods[] = { 17, 19 };
			temp[i] = reconstruct(nums, mods, 2);
		}
		printArray(temp, numCoeff);
		*/
	}


	return 0;
}

