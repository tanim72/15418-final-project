#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <mpi.h>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <chrono>

using namespace std;

/**
 * Local Modified Thomas Algorithm 
 * 
 * @param m Size of the local system (number of equations per process).
 * @param a Sub-diagonal of size m
 * @param b Main diagonal of size m
 * @param c Super-diagonal of size m
 * @param d Right-hand side (RHS) vector of size m
 */
void modifiedThomasAlgorithm(int m, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
    // Normalize first row
    d[0] = d[0] / b[0];
    c[0] = c[0] / b[0];
    a[0] = a[0] / b[0];
    
    // Normalize second row
    d[1] = d[1] / b[1];
    c[1] = c[1] / b[1];
    a[1] = a[1] / b[1];
    
    // Forward elimination 
    for (int i = 2; i < m; i++) {
        double r = 1.0 / (b[i] - a[i] * c[i-1]);
        d[i] = r * (d[i] - a[i] * d[i-1]);
        c[i] = r * c[i];
        a[i] = -r * (a[i] * a[i-1]);
    }
    
    // Backward substitution 
    for (int i = m-3; i >= 1; i--) {
        d[i] = d[i] - c[i] * d[i+1];
        c[i] = -c[i] * c[i+1];
        a[i] = a[i] - c[i] * a[i+1];
    }
    
   
    // unsure if this is the desired result as paper is unclear about the final step
    double r = 1.0 / (1.0 - a[1]*c[0]); // could also be (b[1] - a[1]*c[0]) to avoid NaN

    d[0] = r*(d[0]-a[0]*d[1]);
    c[0] = -r*c[0]*c[1];
    a[0] = r*(a[0]);
}

// Standard Thomas algorithm for solving the reduced system

void standardThomasSolver(int size, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
    vector<double> gamma(size, 0.0);
    
    // Forward elimination
    gamma[0] = c[0] / b[0];
    d[0] = d[0] / b[0];
    
    for (int i = 1; i < size; i++) {
        double denom = b[i] - a[i] * gamma[i-1];
        if (i < size-1) {
            gamma[i] = c[i] / denom;
        }
        d[i] = (d[i] - a[i] * d[i-1]) / denom;
    }
    
    // Back substitution
    for (int i = size-2; i >= 0; i--) {
        d[i] = d[i] - gamma[i] * d[i+1];
    }
}

// Update the remaining unknowns after solving the reduced system

void updateSolution(int m, const vector<double>& a, const vector<double>& c, vector<double>& d, double d_first, double d_last) {
    for (int i = 1; i < m-1; i++) {
        d[i] = d[i] - a[i] * d_first - c[i] * d_last;
    }
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int pid, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 2) {
        if (pid == 0) {
            cout << "Usage: mpirun -np <number_of_processes> " << argv[0] << " <input_file>" << endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    int N = 0;                    
    int m = 0;                   
    vector<double> global_a;      // Global sub-diagonal
    vector<double> global_b;      // Global main diagonal
    vector<double> global_c;      // Global super-diagonal
    vector<double> global_d;      // Global RHS vector
    vector<double> global_x;      // Global solution vector
    
   
    if (pid == 0) {
        string input_file = argv[1];
        ifstream fin(input_file);

        if (!fin) {
            cerr << "Error: Cannot open input file " << input_file << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        // Read N
        fin >> N;

        global_a.resize(N, 0.0);
        global_b.resize(N, 0.0);
        global_c.resize(N, 0.0);
        global_d.resize(N, 0.0);
        global_x.resize(N, 0.0);

        // Read main diagonal
        for (int i = 0; i < N; i++) {
            fin >> global_b[i];
        }

        // Read sub-diagonal
        global_a[0] = 0.0;
        for (int i = 1; i < N; i++) {
            fin >> global_a[i];
        }

        // Read super-diagonal
        for (int i = 0; i < N - 1; i++) {
            fin >> global_c[i];
        }
        global_c[N - 1] = 0.0;

        // Read RHS vector 
        for (int i = 0; i < N; i++) {
            fin >> global_d[i];
        }

        fin.close();
    }

    
    // Broadcast the total size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    const auto compute_start = std::chrono::steady_clock::now();
    
    m = N / size;
    int remainder = N % size;
    
    if (pid < remainder) {
        m++;
    }
    
    // Calculate start index for each process
    int start_idx = pid * (N / size) + min(pid, remainder);
    
    vector<double> local_a(m, 0.0);
    vector<double> local_b(m, 0.0);
    vector<double> local_c(m, 0.0);
    vector<double> local_d(m, 0.0);
    
    // Distribute data to all processes
    int* recvcounts = new int[size];
    int* displs = new int[size];
    
    for (int i = 0; i < size; i++) {
        recvcounts[i] = N / size;
        if (i < remainder) {
            recvcounts[i]++;
        }
        displs[i] = i * (N / size) + min(i, remainder);
    }
    
    MPI_Scatterv(global_b.data(), recvcounts, displs, MPI_DOUBLE, local_b.data(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(global_c.data(), recvcounts, displs, MPI_DOUBLE, local_c.data(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(global_a.data(), recvcounts, displs, MPI_DOUBLE, local_a.data(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(global_d.data(), recvcounts, displs, MPI_DOUBLE, local_d.data(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    // Step 1: Apply the modified Thomas algorithm to local system
    modifiedThomasAlgorithm(m, local_a, local_b, local_c, local_d);
    
    // Step 2: Gather the first and last rows from each process to build the reduced system
    int reduced_size = 2 * size;
    vector<double> reduced_a(reduced_size, 0.0);
    vector<double> reduced_b(reduced_size, 0.0);
    vector<double> reduced_c(reduced_size, 0.0);
    vector<double> reduced_d(reduced_size, 0.0);
    
    // Local coefficients to be gathered
    vector<double> local_coefs(6, 0.0);  // a_1, c_1, d_1, a_m, c_m, d_m
    local_coefs[0] = local_a[0];
    local_coefs[1] = local_c[0];
    local_coefs[2] = local_d[0];
    local_coefs[3] = local_a[m-1];
    local_coefs[4] = local_c[m-1];
    local_coefs[5] = local_d[m-1];
    
    // Gather coefficients to root process
    vector<double> all_coefs;
    if (pid == 0) {
        all_coefs.resize(6 * size, 0.0);
    }
    
    MPI_Gather(local_coefs.data(), 6, MPI_DOUBLE, all_coefs.data(), 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Step 3: Root process constructs and solves the reduced system
    vector<double> reduced_solution;
    if (pid == 0) {
        // Build the reduced system
        for (int i = 0; i < size; i++) {
            int idx1 = 2 * i;     // First equation from process i
            int idx2 = 2 * i + 1; // Last equation from process i
            
            reduced_a[idx1] = all_coefs[6*i + 0];  // a_1
            reduced_b[idx1] = 1.0;               
            reduced_c[idx1] = all_coefs[6*i + 1];  // c_1
            reduced_d[idx1] = all_coefs[6*i + 2];  // d_1
            
            reduced_a[idx2] = all_coefs[6*i + 3];  // a_m
            reduced_b[idx2] = 1.0;                
            if (i < size - 1) {
                reduced_c[idx2] = all_coefs[6*i + 4];  // c_m
            } else {
                reduced_c[idx2] = 0.0; 
            }
            reduced_d[idx2] = all_coefs[6*i + 5];  // d_m
        }
        
        for (int i = 1; i < size; i++) {
            int prev_last = 2 * i - 1;  // Last equation from previous process
            int curr_first = 2 * i;     // First equation from current process
            
            reduced_c[prev_last] = -reduced_a[curr_first];
            reduced_a[curr_first] = -reduced_c[prev_last];
        }
        
        standardThomasSolver(reduced_size, reduced_a, reduced_b, reduced_c, reduced_d);
        
        reduced_solution = reduced_d;
    }
    
    // Step 4: Distribute the solution of the reduced system back to each process
    vector<double> local_reduced_solution(2, 0.0);  // d_1, d_m
    MPI_Scatter(reduced_solution.data(), 2, MPI_DOUBLE, local_reduced_solution.data(), 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Step 5: Update the remaining unknowns
    double d_first = local_reduced_solution[0];
    double d_last = local_reduced_solution[1];
    
    local_d[0] = d_first;
    local_d[m-1] = d_last;
    
    // Update  remaining solutions (d_2 to d_(m-1))
    updateSolution(m, local_a, local_c, local_d, d_first, d_last);
    
    // Gather all solutions back 
    MPI_Gatherv(local_d.data(), m, MPI_DOUBLE, global_x.data(), recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(
                    std::chrono::steady_clock::now() - compute_start).count();
    if (pid == 0)
        std::cout << "Computation time (sec): " << std::fixed << std::setprecision(10) << compute_time << "\n";

    // comment this out if you don't want to print it
    if (pid == 0) {
        cout << "Solution x: ";
        for (int i = 0; i < N; i++) {
            cout << global_x[i] << " ";
        }
        cout << endl;
    }
    
    // Clean up
    delete[] recvcounts;
    delete[] displs;
    
    MPI_Finalize();
    return 0;
}
