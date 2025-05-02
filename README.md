# 15418 Final Project

This repository contains solvers for tridiagonal matrix systems. The test cases and solvers are organized as described below.

## Test Case Format (MPI & OpenMP)

Each test case file is structured as follows:

- **First line:** Size of the matrix (n x n). 
- **Second line:** Main diagonal entries
- **Third line:** Sub-diagonal entries
- **Fourth line:** Super-diagonal entries
- **Fifth line:** RHS (Right Hand Side) vector

If it is a multi-system solver, the first line will have m and n. "m" for the number of systems, and n for the size of each matrix. Then we will have the second-fifth lines for every m that exists.

## Running the Solvers

### Sequential Solver

**Compilation**: `g++ -std=c++11 -o sequentialSolver sequentialSolver.cpp`  
**Execution**: `./sequentialSolver < inputs/{input_name}.txt`

### MPI Recursive Doubling Solver

**Compilation**: `mpic++ -o mpiRecursiveDoubling mpiRecursiveDoubling.cpp`  
**Execution**: `mpirun -np {num_threads} ./mpiRecursiveDoubling inputs/{input_name}.txt`

### MPI Brugnano Solver

**Compilation**: `mpic++ -o mpiBrugnano mpiBrugnano.cpp`  
**Execution**: `mpirun -np {num_threads} ./mpiBrugnano inputs/{input_name}.txt`

### MPI Differential Solver

**Compilation**: `mpicxx -std=c++11 -o diffSolver diffSolver.cpp`
**Execution**: `mpirun -np {num_threads} ./diffSolver inputs/{input_name}.txt`

### OpenMP Recursive Doubling Solver

**Compilation**: `g++ -O3 -fopenmp openMPRecursiveDoubling.cpp -o openMPRecursiveDoubling`
**Execution**: `./openMPRecursiveDoubling -f inputs/{input_name} -n {num_threads}`

### OpenMP Brugnano Solver

**Compilation**: `g++ -O3 -fopenmp openMPBrug.cpp -o openMPBrug`
**Execution**: `./openMPBrug -f inputs/{input_name} -n {num_threads}`

### CUDA Parallel Parallel Cyclic Reduction (PCR)

**Compilation**: `nvcc -O3 -arch=sm_60 -o cudaPCR cudaPCR.cu`
**Execution**: `./cudaPCR inputs/{input_name}.txt`

### CUDA Parallel Cyclic Reduction (CR)

**Compilation**: `nvcc -O3 -arch=sm_60 -o cudaCR cudaCR.cu`
**Execution**: `./cudaCR inputs/{input_name}.txt`

### CUDA Parallel PCR/CR Hybrid

**Compilation**: `nvcc -O3 -arch=sm_60 -o cudaHybridPCR cudaHybridPCR.cu`
**Execution**: `./cudaHybridPCR inputs/{input_name}.txt`