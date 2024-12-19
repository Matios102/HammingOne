# Compiler and flags
CPP = g++
C = gcc
NVCC = nvcc
CFLAGS = -Wall -O2
CPPFLAGS = -Wall -O2 -std=c++11
CUDAFLAGS = -arch=sm_60 -Xcompiler -Wall -O2

# Source files
INPUT_GENERATOR_SRC = input_generator.c
CPP_SOLUTION_SRC = cpp_solution.cpp
CUDA_SOLUTION_SRC = cuda_solution.cu

# Executable files
INPUT_GENERATOR_EXE = input_generator
CPP_SOLUTION_EXE = cpp_solution
CUDA_SOLUTION_EXE = cuda_solution

# Default target
all: run_input run_cpu run_gpu

# Target for generating input
$(INPUT_GENERATOR_EXE): $(INPUT_GENERATOR_SRC)
	$(C) $(CFLAGS) -o $(INPUT_GENERATOR_EXE) $(INPUT_GENERATOR_SRC)

input: $(INPUT_GENERATOR_EXE)

run_input: $(INPUT_GENERATOR_EXE)
	./$(INPUT_GENERATOR_EXE) 1000 100000 10

# Target for CPU implementation
$(CPP_SOLUTION_EXE): $(CPP_SOLUTION_SRC)
	$(CPP) $(CPPFLAGS) -o $(CPP_SOLUTION_EXE) $(CPP_SOLUTION_SRC)

run_cpu: $(CPP_SOLUTION_EXE)
	./$(CPP_SOLUTION_EXE)

# Target for GPU implementation
$(CUDA_SOLUTION_EXE): $(CUDA_SOLUTION_SRC)
	$(NVCC) $(CUDAFLAGS) -o $(CUDA_SOLUTION_EXE) $(CUDA_SOLUTION_SRC)

run_gpu: $(CUDA_SOLUTION_EXE)
	./$(CUDA_SOLUTION_EXE)

# Clean up generated files
clean:
	rm -f $(INPUT_GENERATOR_EXE) $(CPP_SOLUTION_EXE) $(CUDA_SOLUTION_EXE) input_sequences.txt cpu_output_pairs.txt gpu_output_pairs.txt

.PHONY: all run_input run_cpu run_gpu clean
