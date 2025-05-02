

// Differential equation solver using MPI
#include <mpi.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <cmath>
using std::vector;

// Serial Thomas solver for a tridiagonal system:
void thomasSolve(const vector<double>& a,
                 const vector<double>& b,
                 const vector<double>& c,
                 const vector<double>& r,
                 vector<double>& x) {
    int n = b.size();
    if (n == 0) return;
    vector<double> w(n), g(n);
    w[0] = c[0] / b[0];
    g[0] = r[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * w[i-1];
        w[i] = c[i] / denom;
        g[i] = (r[i] - a[i] * g[i-1]) / denom;
    }
    x.resize(n);
    x[n-1] = g[n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = g[i] - w[i] * x[i+1];
    }
}

int main(int argc, char** argv) {
    const auto init_start = std::chrono::steady_clock::now();
    MPI_Init(&argc, &argv);
    int pid;
    int P;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (argc != 2 && pid == 0) {
        std::cerr << "Usage: mpirun -np <P> " << argv[0] << " input.txt\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int N = 0;
    vector<double> A, B, C, R;

    if (pid == 0) {
        std::ifstream fin(argv[1]);
        if (!fin) throw std::runtime_error("Cannot open input file");

        fin >> N;
        B.resize(N);
        R.resize(N);

        // Main diagonal
        for (int i = 0; i < N; ++i) {
            fin >> B[i];
        }
        // Sub-diagonal 
        A.assign(N, 0.0);
        for (int i = 1; i < N; ++i) {
            fin >> A[i];
        }
        // Super-diagonal 
        C.assign(N, 0.0);
        for (int i = 0; i < N-1; ++i) {
            fin >> C[i];
        }

        // RHS (N entries)
        for (int i = 0; i < N; ++i) {
            fin >> R[i];
        }

        if (N % P != 0)
            throw std::runtime_error("N must be divisible by number of processes");
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    const auto compute_start = std::chrono::steady_clock::now();

    if (pid != 0) {
        A.resize(N);
        B.resize(N);
        C.resize(N);
        R.resize(N);
    }

    int M = N / P;                 
    vector<double> a(M), b(M), c(M), r(M);

    // Scatter global to local
    MPI_Scatter(A.data(), M, MPI_DOUBLE, a.data(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B.data(), M, MPI_DOUBLE, b.data(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(C.data(), M, MPI_DOUBLE, c.data(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(R.data(), M, MPI_DOUBLE, r.data(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    vector<double> omega(M), gamma(M), xR(M), xUH(M), xLH(M);

    omega[0] = c[0] / b[0];
    gamma[0] = r[0] / b[0];
    for (int i = 1; i < M; ++i) {
        double denom = b[i] - a[i] * omega[i-1];
        omega[i] = c[i] / denom;
        gamma[i] = (r[i] - a[i] * gamma[i-1]) / denom;
    }

    xR[M-1] = gamma[M-1];
    xLH[M-1] = -omega[M-1];
    xUH[M-1] = a[M-1] / b[M-1];
    for (int i = M-2; i >= 0; --i) {
        xR[i] = gamma[i] - omega[i] * xR[i+1];
        xLH[i] = -omega[i] * xLH[i+1];
        double denom = b[i] - c[i] * xUH[i+1];
        xUH[i] = -a[i] / denom;
    }
    xUH[0] = -xUH[0];
    for (int i = 1; i < M; ++i) {
        xUH[i] = -xUH[i] * xUH[i-1];
    }
    double uhc = 0.0;
    double lhc = 0.0;

    if (P > 1) {
        // Build & solve reduced system via log2P Sendrecv passes 
        int log2P = 0; 
        while ((1<<log2P) < P) ++log2P;
        int total = 8 * (1<<log2P);
        vector<double> out(total, 0.0);

        // pack 8 values for each process
        out[0] = -1.0;      
        out[1] = xUH[0];
        out[2] = xLH[0];    
        out[3] = -xR[0];
        out[4] = xUH[M-1];  
        out[5] = xLH[M-1];
        out[6] = -1.0;      
        out[7] = -xR[M-1];

        for (int step = 0; step < log2P; ++step) {
            int chunk = 8 * (1<<step);
            int left  = (pid - (1<<step) + P) % P;
            int right = (pid + (1<<step)) % P;
            MPI_Sendrecv(
                out.data(), chunk, MPI_DOUBLE, left,  0,
                out.data()+chunk, chunk, MPI_DOUBLE, right, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE
            );
        }

        // extract 2P-2 sized reduced system
        int Rdim = 2*P - 2;
        vector<double> ra(Rdim), rb(Rdim), rc(Rdim), rr(Rdim);
        int base0 = 8*(P-pid) + 4;
        for (int i = 0; i < Rdim; ++i) {
            int idx = (base0 + 4*i) % (8*P);
            ra[i] = out[idx];
            rb[i] = out[idx+1];
            rc[i] = out[idx+2];
            rr[i] = out[idx+3];
        }

        vector<double> coeffs;
        thomasSolve(ra, rb, rc, rr, coeffs);

        if (pid > 0) uhc = coeffs[2*pid - 2];
        if (pid < P-1) lhc = coeffs[2*pid - 1];
    }

    // Gather the solution
    vector<double> xloc(M);
    for (int i = 0; i < M; ++i) {
        xloc[i] = xR[i] + uhc * xUH[i] + lhc * xLH[i];
    }

    vector<double> X;
    if (pid == 0) X.resize(N);
    MPI_Gather(xloc.data(), M, MPI_DOUBLE, X.data(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
    double total_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();

    if (pid == 0) {
        std::cout << "Computation time (sec): " << std::fixed << std::setprecision(10) << compute_time << "\n";
        std::cout << "Total time (sec): " << std::fixed << std::setprecision(10) << total_time << "\n";
        // comment out below for large N
        std::cout << "Solution x:\n";
        for (int i = 0; i < N; ++i){
            std::cout << X[i] << (i+1==N? "\n":" ");
        }
    }

    MPI_Finalize();
    return 0;
}