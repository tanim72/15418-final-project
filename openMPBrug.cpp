

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <unistd.h>   
#include <getopt.h>  
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <omp.h>
#include <immintrin.h> 

using namespace std;

//Local Modified Thomas Algorithm 
void modifiedThomasAlgorithm(int m, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
    // Normalize first row
    d[0] = d[0] / b[0];
    c[0] = c[0] / b[0];
    a[0] = a[0] / b[0];
    
    // Normalize second row
    if (m > 1) {
        d[1] = d[1] / b[1];
        c[1] = c[1] / b[1];
        a[1] = a[1] / b[1];
    }
    
    // Forward elimination
    for (int i = 2; i < m; i++) {
        // Prefetch data for next iteration (going backward)
        if (i + 2 < m) {
            _mm_prefetch((const char*)&a[i+2], _MM_HINT_T0);
            _mm_prefetch((const char*)&b[i+2], _MM_HINT_T0);
            _mm_prefetch((const char*)&c[i+2], _MM_HINT_T0);
            _mm_prefetch((const char*)&d[i+2], _MM_HINT_T0);
        }
       
        if (i - 2 >= 0) {
            _mm_prefetch((const char*)&a[i-2], _MM_HINT_T0);
            _mm_prefetch((const char*)&c[i-2], _MM_HINT_T0);
            _mm_prefetch((const char*)&d[i-2], _MM_HINT_T0);
        }
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
    
    if (m >= 2) {
        double r = 1.0 / (1.0 - a[1]*c[0]);
        d[0] = r * (d[0] - a[0] * d[1]);
        c[0] = -r * c[0] * c[1];
        a[0] = r * a[0];
    }
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
        // Prefetch data for next iteration (going backward)
        if (i - 2 >= 0) {
            _mm_prefetch((const char*)&gamma[i-2], _MM_HINT_T0);
            _mm_prefetch((const char*)&d[i-2], _MM_HINT_T0);
        }
        
        double gamma_i = gamma[i];
        double d_next = d[i+1];
        d[i] = d[i] - gamma_i * d_next;
    }
}

int main(int argc, char* argv[]) {
    auto init = chrono::steady_clock::now();
    string input_file;
    int num_threads = 0;

    int opt;
    while ((opt = getopt(argc, argv, "f:n:")) != -1) {
        switch (opt) {
            case 'f':
                input_file = optarg;
                break;
            case 'n':
                num_threads = atoi(optarg);
                break;
            default:
                cerr << "Usage: " << argv[0] << " -f <input_file> -n <num_threads>\n";
                return 1;
        }
    }

    if (input_file.empty() || num_threads <= 0) {
        cerr << "Usage: " << argv[0] << " -f <input_file> -n <num_threads>\n";
        return 1;
    }

    omp_set_num_threads(num_threads);

    ifstream fin(input_file);
    if (!fin) {
        cerr << "Error: Cannot open input file " << input_file << "\n";
        return 1;
    }

    int N;
    fin >> N;

    vector<double> global_a(N), global_b(N), global_c(N), global_d(N), global_x(N);

    // Read diagonals and RHS
    for (int i = 0; i < N; i++) fin >> global_b[i];
    global_a[0] = 0.0;
    for (int i = 1; i < N; i++) fin >> global_a[i];
    for (int i = 0; i < N-1; i++) fin >> global_c[i];
    global_c[N-1] = 0.0;
    for (int i = 0; i < N; i++) fin >> global_d[i];
    fin.close();

    cout << "Running with " << num_threads << " threads" << endl;
    auto t0 = chrono::steady_clock::now();

    // Partition the work
    vector<int> chunk_sizes(num_threads), start_indices(num_threads);
    int base = N / num_threads, rem = N % num_threads;
    int cur = 0;
    for (int i = 0; i < num_threads; i++) {
        chunk_sizes[i]   = base + (i < rem ? 1 : 0);
        start_indices[i] = cur;
        cur += chunk_sizes[i];
    }

    //  buffer to collect boundary coefficients
    vector<vector<double>> all_coefs(num_threads, vector<double>(6));

    #pragma omp parallel
    {
        int tid       = omp_get_thread_num();
        int start_idx = start_indices[tid];
        int m         = chunk_sizes[tid];

        vector<double> local_a(m), local_b(m), local_c(m), local_d(m);
        for (int i = 0; i < m; i++) {
            // Prefetch the next data elements
            if (i + 4 < m) {
                int gi_prefetch = start_idx + i + 4;
                _mm_prefetch((const char*)&global_a[gi_prefetch], _MM_HINT_T0);
                _mm_prefetch((const char*)&global_b[gi_prefetch], _MM_HINT_T0);
                _mm_prefetch((const char*)&global_c[gi_prefetch], _MM_HINT_T0);
                _mm_prefetch((const char*)&global_d[gi_prefetch], _MM_HINT_T0);
            }
            
            int gi = start_idx + i;
            local_a[i] = global_a[gi];
            local_b[i] = global_b[gi];
            local_c[i] = global_c[gi];
            local_d[i] = global_d[gi];
        }

        if (m == 1) {
            local_d[0] /= local_b[0];
            local_a[0] = local_c[0] = 0.0;
        } else {
            modifiedThomasAlgorithm(m, local_a, local_b, local_c, local_d);
        }

        // reduced system coefficients
        all_coefs[tid][0] = local_a[0];
        all_coefs[tid][1] = local_c[0];
        all_coefs[tid][2] = local_d[0];
        all_coefs[tid][3] = local_a[m-1];
        all_coefs[tid][4] = local_c[m-1];
        all_coefs[tid][5] = local_d[m-1];

        #pragma omp barrier
        #pragma omp single
        {
            int R = 2 * num_threads;
            vector<double> ra(R), rb(R,1.0), rc(R), rd(R);

            // building reduced system
            for (int i = 0; i < num_threads; i++) {
                int e1 = 2*i, e2 = 2*i+1;
                
                // Prefetch next coefficients
                if (i + 1 < num_threads) {
                    _mm_prefetch((const char*)&all_coefs[i+1][0], _MM_HINT_T0);
                }
                
                ra[e1] = all_coefs[i][0]; 
                rc[e1] = all_coefs[i][1]; 
                rd[e1] = all_coefs[i][2];
                ra[e2] = all_coefs[i][3]; 
                rc[e2] = (i < num_threads-1 ? all_coefs[i][4] : 0.0);
                rd[e2] = all_coefs[i][5];
            }
            // link the boundaries
            for (int i=1; i<num_threads; ++i) {
                int prev = 2*i-1, nxt=2*i;
                rc[prev]  = -ra[nxt];
                ra[nxt]   = -rc[prev];
            }
            standardThomasSolver(R, ra, rb, rc, rd);
            // put back values
            for (int i=0; i<num_threads; ++i) {
                int s = start_indices[i], e = s+chunk_sizes[i]-1;
                global_x[s] = rd[2*i];
                global_x[e] = rd[2*i+1];
            }
        }
        #pragma omp barrier

        int s = start_idx;
        double d0 = global_x[s], dN = global_x[s+m-1];
        for (int i = 1; i < m-1; i++) {
            // Prefetch next interior values
            if (i + 2 < m - 1) {
                _mm_prefetch((const char*)&local_d[i+2], _MM_HINT_T0);
                _mm_prefetch((const char*)&local_a[i+2], _MM_HINT_T0);
                _mm_prefetch((const char*)&local_c[i+2], _MM_HINT_T0);
            }
            
            double local_d_i = local_d[i];
            double local_a_i = local_a[i];
            double local_c_i = local_c[i];
            
            global_x[s+i] = local_d_i - local_a_i * d0 - local_c_i * dN;
        }
    }

    double secs = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t0).count();
    double totalTime = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - init).count();
    cout << "Computation time (sec): " << fixed << setprecision(10) << secs << "\n";
    cout << "Total time (sec): " << fixed << setprecision(10) << totalTime << "\n";

    // comment out below for large N
    cout << "Solution x: ";
    for (int i = 0; i < N; i++) cout << global_x[i] << " ";
    cout << endl;

    return 0;
}