# 15418 Final Project

This repository contains solvers for tridiagonal matrix systems. The test cases and solvers are organized as described below.

## Test Case Format

Each test case file is structured as follows:

- **First line:** Size of the matrix (n x n)
- **Second line:** Main diagonal entries
- **Third line:** Sub-diagonal entries
- **Fourth line:** Super-diagonal entries
- **Fifth line:** RHS (Right Hand Side) vector

## Running the Solvers

### Sequential Solver

**Compilation**: `g++ -std=c++11 -o sequentialSolver sequentialSolver.cpp`  
**Execution**: `./sequentialSolver < inputs/{input_name}.txt`

### MPI Solver (Recursive Doubling)

**Compilation**: `mpic++ -o mpiRecursiveDoubling mpiRecursiveDoubling.cpp`  
**Execution**: `mpirun -np {num_threads} ./mpiRecursiveDoubling inputs/{input_name}.txt`

### MPI Brugnano Solver

**Compilation**: `mpic++ -o mpiBrugnano mpiBrugnano.cpp`  
**Execution**: `mpirun -np {num_threads} ./mpiBrugnano inputs/{input_name}.txt`
