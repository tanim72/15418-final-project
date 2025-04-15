#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <chrono>

using namespace std;


void ThomasAlgorithm_P(int pid, int nproc, int N, double *a, double *b, double *c, double *q, double *x) {
    
    const double EPSILON = 1e-15;
    int i;
    double *l = new double[N];
    double *d = new double[N];
    double *y = new double[N];
    
    for(i = 0; i < N; i++){
        l[i] = d[i] = y[i] = 0.0;
    }
    
    double R[2][2], U[2][2];
    R[0][0] = 1.0; R[0][1] = 0.0;
    R[1][0] = 0.0; R[1][1] = 1.0;
    
    
    int localRows = N / nproc;
    int offset = pid * localRows;
    
    // Stage 1: Compute local coefficients.
    
    if (pid == 0) {
        // Special case for first row on process 0
        double tmpVal = a[offset] * R[0][0];
        R[1][0] = R[0][0];  
        R[1][1] = R[0][1]; 
        R[0][1] = a[offset] * R[0][1];
        R[0][0] = tmpVal;
        
        
        // Process remaining rows
        for(i = 1; i < localRows; i++){
            int idx = i + offset;
            double a_val = a[idx];
            double mult = b[idx - 1] * c[idx - 1];
            
            
            double tmp0 = a_val * R[0][0] - mult * R[1][0];
            double tmp1 = a_val * R[0][1] - mult * R[1][1];
            
          
            
            R[1][0] = R[0][0];
            R[1][1] = R[0][1];
            R[0][0] = tmp0;
            R[0][1] = tmp1;

            // Normalization
            double scale = std::max({fabs(R[0][0]), fabs(R[0][1]), fabs(R[1][0]), fabs(R[1][1])});
            if (scale > 0.0) {
                R[0][0] /= scale;
                R[0][1] /= scale;
                R[1][0] /= scale;
                R[1][1] /= scale;
            }

        }
    } else {
        // For processes other than 0
        for(i = 0; i < localRows; i++){
            int idx = i + offset;
            if (idx > 0) {  
                double a_val = a[idx];
                double mult = b[idx - 1] * c[idx - 1];
                
                
                double tmp0 = a_val * R[0][0] - mult * R[1][0];
                double tmp1 = a_val * R[0][1] - mult * R[1][1];
                
               
                
                R[1][0] = R[0][0]; 
                R[1][1] = R[0][1];  
                R[0][0] = tmp0;
                R[0][1] = tmp1;

                // Normalization
                double scale = std::max({fabs(R[0][0]), fabs(R[0][1]), fabs(R[1][0]), fabs(R[1][1])});
                if (scale > 0.0) {
                    R[0][0] /= scale;
                    R[0][1] /= scale;
                    R[1][0] /= scale;
                    R[1][1] /= scale;
                }

            } 
        }
    }
    
    
    // Apply scailing 
    double max_val = std::max(std::max(std::abs(R[0][0]), std::abs(R[0][1])), 
                             std::max(std::abs(R[1][0]), std::abs(R[1][1])));
    
    if (max_val > 1e10 || max_val < 1e-10) {
        double scale = (max_val > 1.0) ? 1.0/max_val : 1.0;
        R[0][0] *= scale;
        R[0][1] *= scale;
        R[1][0] *= scale;
        R[1][1] *= scale;
    }
    
    // Stage 2: Recursive-Doubling Communication.
  
    int stages = (int)log2(nproc);
    for(i = 0; i < stages; i++){
        
        int partner = pid + (1<<i);
        if (partner < nproc) {
            MPI_Send(&R[0][0], 4, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
        }
        
        partner = pid - (1<<i);
        if (partner >= 0) {
            MPI_Recv(&U[0][0], 4, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            
            double tmp0 = R[0][0]*U[0][0] + R[0][1]*U[1][0];
            R[0][1] = R[0][0]*U[0][1] + R[0][1]*U[1][1];
            R[0][0] = tmp0;
            tmp0 = R[1][0]*U[0][0] + R[1][1]*U[1][0];
            R[1][1] = R[1][0]*U[0][1] + R[1][1]*U[1][1];
            R[1][0] = tmp0;
            
           
            
            
            max_val = std::max(std::max(std::abs(R[0][0]), std::abs(R[0][1])), 
                              std::max(std::abs(R[1][0]), std::abs(R[1][1])));
            
            if (max_val > 1e10 || max_val < 1e-10) {
                double scale = (max_val > 1.0) ? 1.0/max_val : 1.0;
                R[0][0] *= scale;
                R[0][1] *= scale;
                R[1][0] *= scale;
                R[1][1] *= scale;
            }
        }
    }
    
    
    // Stage 3: Compute local boundary coefficient d
    
    if (nproc > 1) {
        double denom = R[1][0] + R[1][1];
        if (std::abs(denom) < EPSILON) {
            denom = (denom >= 0) ? EPSILON : -EPSILON;
        }
        
        d[offset + localRows - 1] = (R[0][0] + R[0][1]) / denom;
        
        if (pid == 0) {
            MPI_Send(&d[offset + localRows - 1], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&d[offset - 1], 1, MPI_DOUBLE, pid - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            if (pid != nproc - 1) {
                MPI_Send(&d[offset + localRows - 1], 1, MPI_DOUBLE, pid + 1, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        // Single process case
        double denom = R[1][0] + R[1][1];
        if (std::abs(denom) < EPSILON) {
            denom = (denom >= 0) ? EPSILON : -EPSILON;
        }
        
        d[offset + localRows - 1] = (R[0][0] + R[0][1]) / denom;
        
    }
    
    // Stage 4: Compute local arrays d and l.
    
    if (pid == 0) {
        l[0] = 0.0;
        d[0] = a[0];
        
        for (i = 1; i < localRows - 1; i++){
            int idx = offset + i;
            if (std::abs(d[idx - 1]) < EPSILON) {
                l[idx] = b[idx - 1] / (d[idx - 1] + ((d[idx - 1] >= 0) ? EPSILON : -EPSILON));
            } else {
                l[idx] = b[idx - 1] / d[idx - 1];
            }
            
            d[idx] = a[idx] - l[idx] * c[idx - 1];
            
            
            if (std::isnan(l[idx]) || std::isnan(d[idx])) {
            }
        }
        
        int idx = offset + localRows - 1;
        if (std::abs(d[idx - 1]) < EPSILON) {
            l[idx] = b[idx - 1] / (d[idx - 1] + ((d[idx - 1] >= 0) ? EPSILON : -EPSILON));
        } else {
            l[idx] = b[idx - 1] / d[idx - 1];
        }
        
        
       
    } else {
        for (i = 0; i < localRows - 1; i++){
            int idx = offset + i;
            if (std::abs(d[idx - 1]) < EPSILON) {
                l[idx] = b[idx - 1] / (d[idx - 1] + ((d[idx - 1] >= 0) ? EPSILON : -EPSILON));
            } else {
                l[idx] = b[idx - 1] / d[idx - 1];
            }
            
            d[idx] = a[idx] - l[idx] * c[idx - 1];
            
            
          
        }
        
        int idx = offset + localRows - 1;
        if (std::abs(d[idx - 1]) < EPSILON) {
            l[idx] = b[idx - 1] / (d[idx - 1] + ((d[idx - 1] >= 0) ? EPSILON : -EPSILON));
        } else {
            l[idx] = b[idx - 1] / d[idx - 1];
        }
        
        
        
    }
    
    if (pid > 0)
        d[offset - 1] = 0.0; // adjust boundary condition
    
    
    // Stage 5: Distribute the full d and l arrays.
    
    double *tmpArr = new double[N];
    for (i = 0; i < N; i++)
        tmpArr[i] = d[i];
    MPI_Allreduce(tmpArr, d, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (i = 0; i < N; i++)
        tmpArr[i] = l[i];
    MPI_Allreduce(tmpArr, l, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] tmpArr;
    
    
    // Stage 6: Forward and Backward Substitution.
    
    if (pid == 0) {
        // Forward substitution
        y[0] = q[0];
        
        for (i = 1; i < N; i++){
            y[i] = q[i] - l[i] * y[i - 1];
            
            if (i % 100 == 0 || std::isnan(y[i])) {
            }
            
            
        }
        
        // Backward substitution
        if (std::abs(d[N - 1]) < EPSILON) {
            x[N - 1] = y[N - 1] / (d[N - 1] + ((d[N - 1] >= 0) ? EPSILON : -EPSILON));
        } else {
            x[N - 1] = y[N - 1] / d[N - 1];
        }
        
        
        for (i = N - 2; i >= 0; i--){
            if (std::abs(d[i]) < EPSILON) {
                x[i] = (y[i] - c[i] * x[i + 1]) / (d[i] + ((d[i] >= 0) ? EPSILON : -EPSILON));
            } else {
                x[i] = (y[i] - c[i] * x[i + 1]) / d[i];
            }
            
        }
        
    }
    
    delete[] l;
    delete[] d;
    delete[] y;
    
}


int main(int argc, char *argv[]) {
    int pid, nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    if (pid == 0) {
        printf("Starting MPI Parallel Thomas Algorithm Solver with %d processes\n", nproc);
    }
    
    int N; // system size read from file
    double *a, *b, *c, *q, *x;
    
    // Make sure an input file was specified.
    if(argc < 2) {
        if(pid == 0)
            cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        MPI_Finalize();
        return 1;
    }
    
    // Process 0 reads from the input file.
    if(pid == 0) {
        printf("Process %d: Reading input file: %s\n", pid, argv[1]);
        ifstream fin(argv[1]);
        if(!fin) {
            cerr << "Error: Cannot open file " << argv[1] << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        fin >> N;
        printf("Process %d: System size N = %d\n", pid, N);
        
        // Check if N is compatible with number of processes
        if (N % nproc != 0) {
            printf("ERR: N (%d) is not divisible by number of processes (%d)\n", N, nproc);
        }
        
        // Check if number of processes is a power of 2
        int log2_p = (int)log2(nproc);
        if ((1 << log2_p) != nproc) {
            printf("WARNING: Number of processes (%d) is not a power of 2\n", nproc);
        }
        
        a = new double[N];
        b = new double[N];
        c = new double[N];
        q = new double[N];
        
        for (int i = 0; i < N; i++) {
            fin >> a[i];
           
        }
        
        for (int i = 0; i < N - 1; i++) {
            fin >> b[i];
            
        }
        b[N - 1] = 0.0;  // dump
        
        for (int i = 0; i < N - 1; i++) {
            fin >> c[i];
           
        }
        c[N - 1] = 0.0;  // dump
        
        for (int i = 0; i < N; i++) {
            fin >> q[i];
            
        }
        fin.close();
        
       
    }

    
    
    // Broadcast the system size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("Process %d: Received N = %d\n", pid, N);
    
    // All processes allocate arrays
    if(pid != 0) {
        a = new double[N];
        b = new double[N];
        c = new double[N];
        q = new double[N];
    }
    
    // Broadcast the matrix and vector data
    MPI_Bcast(a, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(c, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(q, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    printf("Process %d: Received all data via broadcast\n", pid);
    
    const auto compute_start = std::chrono::steady_clock::now();

    x = new double[N];
    for (int i = 0; i < N; i++) {
        x[i] = 0.0;
    }
    
    ThomasAlgorithm_P(pid, nproc, N, a, b, c, q, x);

    double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(
        std::chrono::steady_clock::now() - compute_start).count();
    if (pid == 0)
        std::cout << "Computation time (sec): " << std::fixed << std::setprecision(10) << compute_time << "\n";
    
    // comment out if you don't want to print it
    if (pid == 0) {
        cout << "Solution vector x:" << endl;
        for (int i = 0; i < N; i++) {
            cout << x[i] << " ";
        }
        cout << endl;
    }
    
    // Clean up memory.
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] q;
    delete[] x;
    
    MPI_Finalize();
    return 0;
}
