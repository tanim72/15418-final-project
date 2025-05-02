

// openMP recursive doubling implementation of Thomas algorithm
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <chrono>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include <iomanip>

void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog << " -f <input_file> -n <num_threads>\n";
  exit(EXIT_FAILURE);
}

void ThomasAlgorithm_OMP(int N,
                         double* a, double* b, double* c,
                         double* q, double* x,
                         int num_threads)
{
    const double EPS = 1e-15;
    int localRows = N / num_threads;

    std::vector<std::array<double,4>> Rstore(num_threads);
    std::vector<double> d(N), l(N);

    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();

       
        int base = N / num_threads;
        int rem  = N % num_threads;
        int localRows = (tid < rem ? base+1 : base);
        int offset = (tid < rem ? tid*(base+1) : rem*(base+1) + (tid-rem)*base);

        
        int prefetch_distance = 8; // can be tuned

        // stage 1: local R accumulation 
        double R00=1.0, R01=0.0, R10=0.0, R11=1.0;
        for(int i = 0; i < localRows; ++i) {
            int idx = offset + i;

            double aval = a[idx];
            double mult = (idx>0 ? b[idx-1]*c[idx-1] : 0.0);
            double tmp0 = aval*R00 - mult*R10;
            double tmp1 = aval*R01 - mult*R11;
            R10=R00; R11=R01;
            R00=tmp0; R01=tmp1;

            double sc = std::max({fabs(R00), fabs(R01), fabs(R10), fabs(R11)});
            if(sc>0){ R00/=sc; R01/=sc; R10/=sc; R11/=sc; }
        }
        Rstore[tid] = {R00,R01,R10,R11};

        // stage 2: recursive doubling 
        int stages = (int)std::log2(num_threads);
        for(int s=0; s<stages; ++s) {
            int dist = 1<<s;
            auto myR = Rstore[tid];
            int partner = (tid>=dist ? tid-dist : tid+dist);
            auto pR = Rstore[partner];

            double n00 = myR[0]*pR[0] + myR[1]*pR[2];
            double n01 = myR[0]*pR[1] + myR[1]*pR[3];
            double n10 = myR[2]*pR[0] + myR[3]*pR[2];
            double n11 = myR[2]*pR[1] + myR[3]*pR[3];

            double sc = std::max({fabs(n00), fabs(n01), fabs(n10), fabs(n11)});
            if(sc>0){ n00/=sc; n01/=sc; n10/=sc; n11/=sc; }

            #pragma omp barrier
            Rstore[tid] = {n00,n01,n10,n11};
            #pragma omp barrier
        }

        // stage 3: boundary d at end of segment 
        {
            auto &R = Rstore[tid];
            int boundary = offset + localRows - 1;
            double denom = R[2] + R[3];
            if(fabs(denom)<EPS) denom = (denom>=0?EPS:-EPS);
            d[boundary] = (R[0] + R[1]) / denom;
        }
        #pragma omp barrier

        // stage 4: build full d[] and l[] in parallel 
        if(tid==0) {
            l[0]=0; d[0]=a[0];
            for(int i=1;i<localRows;++i){
                double dp=d[i-1];
                double dm = fabs(dp)<EPS ? dp+(dp>=0?EPS:-EPS) : dp;
                l[i] = b[i-1]/dm;
                d[i] = a[i] - l[i]*c[i-1];
            }
        } else {
            for(int i=0;i<localRows;++i){
                int idx = offset+i;

                if (idx + prefetch_distance < N) {
                    __builtin_prefetch(&b[idx + prefetch_distance - 1], 0, 1);
                    __builtin_prefetch(&c[idx + prefetch_distance - 1], 0, 1);
                    __builtin_prefetch(&a[idx + prefetch_distance], 0, 1);
                }

                double dp = d[idx-1];
                double dm = fabs(dp)<EPS ? dp+(dp>=0?EPS:-EPS) : dp;
                l[idx] = b[idx-1]/dm;
                d[idx] = a[idx] - l[idx]*c[idx-1];
            }
        }
        #pragma omp barrier

        // stage 6: forward/backward (single thread) 
        #pragma omp single
        {
            std::vector<double> y(N);
            y[0] = q[0];
            for(int i=1;i<N;++i) {
                y[i] = q[i] - l[i]*y[i-1];
            }

            x[N-1] = y[N-1] /
                (fabs(d[N-1])<EPS ? d[N-1]+(d[N-1]>=0?EPS:-EPS) : d[N-1]);
            for(int i=N-2;i>=0;--i){
                double dm = fabs(d[i])<EPS ? d[i]+(d[i]>=0?EPS:-EPS) : d[i];
                x[i] = (y[i] - c[i]*x[i+1]) / dm;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    auto tinit = std::chrono::steady_clock::now();
    std::string input_file;
    int num_threads = 0;
    int opt;
    while((opt = getopt(argc, argv, "f:n:")) != -1) {
      switch(opt) {
        case 'f': input_file = optarg; break;
        case 'n': num_threads = atoi(optarg); break;
        default: print_usage(argv[0]);
      }
    }
    if(input_file.empty() || num_threads<1)
      print_usage(argv[0]);

    std::cout << "Number of threads: " << num_threads << "\n";

    std::ifstream fin(input_file);
    if(!fin) {
      std::cerr<<"Cannot open "<<input_file<<"\n";
      return 1;
    }
    int N; fin>>N;
    std::vector<double> a(N), b(N), c(N), q(N), x(N,0.0);
    for(int i=0;i<N;++i) fin>>a[i];
    for(int i=0;i<N-1;++i) fin>>b[i]; b[N-1]=0;
    for(int i=0;i<N-1;++i) fin>>c[i]; c[N-1]=0;
    for(int i=0;i<N;++i) fin>>q[i];
    fin.close();

    auto t0 = std::chrono::steady_clock::now();
    ThomasAlgorithm_OMP(N, a.data(), b.data(), c.data(), q.data(), x.data(), num_threads);
    double t = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    double total_time = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tinit).count();

    std::cout << "Computation time (sec): "
              << std::fixed << std::setprecision(10) << t << "\n";
    
    std::cout << "Total time (sec): "
              << std::fixed << std::setprecision(10) << total_time << "\n";
    std::cout << "Solution x:\n";
    for(double xi: x) std::cout<< xi << " ";
    std::cout<<"\n";
    return 0;
}


