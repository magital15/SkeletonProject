NVCC = /usr/local/cuda/bin/nvcc
NVCC_FLAGS = -g -G -Xcompiler -Wall

all: main.exe

main.exe: main.o kernel.o initializers.o tests.o kernelGPU.o
	$(NVCC) $^ -o $@

main.o: main.cpp kernel.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

kernel.o: kernel.cu kernel.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

initializers.o: initializers.cu kernel.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ 

tests.o: tests.cu kernel.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ 

kernelGPU.o: kernelGPU.cu kernel.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@
