# 15418-final-project

Format of test cases:
# first number: size of matrix (nxn)
# second row: main diagonal
# third row: super diagonal
# fourth row: subdiagonal
# fifth row: RHS vector


Running sequentialSolver:
Compile: g++ -std=c++11 -o mpiSolver mpiSolver.cpp
Run: ./sequentialSolver < inputs/{input_name}.txt

Running mpiSolver:
Compile: mpic++ -o mpiSolver mpiSolver.cpp
Run: mpirun -np {num_threads} ./mpiSolver inputs/{input_name}.txt