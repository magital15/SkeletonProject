#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

const bool testing = true;
const bool multipleCoeff = false;
//extern istream cin;

int modInverse(int a, int m);
int gcdExtended(int a, int b, int *x, int *y);
int findCoeff(int *shortenedRes, int *mods);
//int findCoeff(int a12, int modProduct, int counter, int *shortenedRes, int *mods, int arrSize);

// Function to find modulo inverse of a
// PRECONDITION: a and m are coprime
int modInverse(int a, int m)
{
	int x, y;
	int g = gcdExtended(a, m, &x, &y);
	if (g != 1)
		printf("Inverse doesn't exist");
	else
	{
		// m is added to handle negative x
		int res = (x%m + m) % m;
		printf("%d * %d = 1 (mod %d)\n", a, res, m);
		printf("So %d is the multiplicative inverse of %d (mod %d)\n", res, a, m);
		return res;
	}
	return 0;
}

// C function for extended Euclidean Algorithm
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
int findCoeff(int *shortenedRes, int *mods, int arrSize) {
	int a1 = shortenedRes[0];
	int a2 = shortenedRes[1];
	int p1 = mods[0];
	int p2 = mods[1];

	// should return int2 rather than individual value
	int k1modp2 = modInverse(p1, p2);
	int k2modp1 = modInverse(p2, p1);

	int a12 = (a2*k1modp2*p1 + a1*k2modp1*p2) % (p1*p2);
	
	int counter = 2;
	while (counter < numMods) {
		int a1 = a12;
		int a2 = shortenedRes[counter];
		int p1 = p1*p2;
		int p2 = mods[counter];
		k1modp2 = modInverse(p1, p2);
		k2modp1 = modInverse(p2, p1);
		a12 = (a2*k1modp2*p1 + a1*k2modp1*p2) % (p1*p2);
	}

	return a12;
}

/*
// recursive method that takes two array pointers and pointer at next calculated index
int findCoeff(int a12, int modProduct, int counter, int *shortenedRes, int *mods, int arrSize) {
	if (counter >= arrSize) {
		return a12;
	}
	int a1 = a12;
	int a2 = shortenedRes[counter];
	int p1 = modProduct;
	int p2 = mods[counter];

	int new_a12 = (a2*modInverse(p1, p2)*p1 + a1*modInverse(p2, p1)*p2) % (p1*p2);

	return findCoeff(new_a12, p1*p2, ++counter, shortenedRes, mods, arrSize);
}*/


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
			outputCoeff[i] = findCoeff(shortenedRes[i], mods, numMods);
		}
	}
	else if (testing) {
		int mods[] = { 5, 7, 13, 17, 19 };
		int clues[] = { 0, 4, 12, 8, 6 };
		int result = findCoeff(mods, clues, numMods);
		printf("Riddle answer %d", result);
	}
	//char wait;
	//cin>>wait; // keeps console window open during debug mode
	return 0;
}