#ifndef STRUCTHELPERFUNCTIONS_H
#define STRUCTHELPERFUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <sstream>
#include <iostream>
using namespace std;

#define prime0  7
#define prime1  11
#define prime2  13
#define prime3  17
#define prime4  19

typedef struct polDense_t {
	int base;
	int length;
	int arr[0];
}*polDense;

typedef struct poly_t {
	polDense members[6];
}*poly;

polDense init(int mod, int *input, int input_length);
polDense copy(polDense p);
int modAllCoeff(polDense p, int newMod);
polDense copyMod(polDense p, int newMod);
poly init(int *input, int input_length);
void printArray(int *arr, int size);
void printArray(polDense p);

polDense init(int mod, int *input, int input_length) {
	polDense p = (polDense)malloc(sizeof(int)*(input_length + 2));
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

#endif