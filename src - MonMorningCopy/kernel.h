#ifndef KERNEL_H
#define KERNEL_H
#define numCoeff 5
#define numMods 7

#define prime0  7
#define prime1  11
#define prime2  13
#define prime3  17
#define prime4  19

#include <sstream>
#include <iostream>
using namespace std;

void innerProduct(float * a, float * b, float * output);

//polDense init(int mod, int *input, int input_length);
//polDense copy(polDense p);
//int modAllCoeff(polDense p, int newMod);
//polDense copyMod(polDense p, int newMod);
//poly init(int *input, int input_length);



#endif