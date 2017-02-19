#include "kernel.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#define TPB 64
__global__ void innerProductKernel(float *d_out, float *a, float *b)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	//if (i >= N) return;

}

void innerProduct(float * a, float * b, float * output)
{

}