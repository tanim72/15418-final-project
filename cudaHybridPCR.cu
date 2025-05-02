// CUDA implementation of the cyclic reduction algorithm
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cuda_runtime.h>
#include "CycleTimer.h"
#include <iostream>

#define CUDA_CHECK(call)                                                 \
  do {                                                                   \
    cudaError_t err = call;                                              \
    if (err != cudaSuccess) {                                            \
      fprintf(stderr,"%s:%d %s\n",                                       \
              __FILE__,__LINE__,cudaGetErrorString(err));                \
      exit(1);                                                           \
    }                                                                    \
  } while (0)


#define CR_LEVELS 1                 // CR_LEVELS will keep increasing by a factor of 2
constexpr int TPB = 256;            // threads per block

__device__ __host__ inline int log2ceil(int x) {
    int l = 0;
    while ((1 << l) < x) ++l;
    return l;
}

template<int TILE>
__global__ void hybrid_solve_multi_tiled(int N,
                                         const float* __restrict__ a_in,
                                         const float* __restrict__ b_in,
                                         const float* __restrict__ c_in,
                                         const float* __restrict__ d_in,
                                         float*       __restrict__ x_out)
{
    extern __shared__ float s[];
    float* sa = s + 0*TILE;
    float* sb = s + 1*TILE;
    float* sc = s + 2*TILE;
    float* sd = s + 3*TILE;
    float* sx = s + 4*TILE;

    int tid   = threadIdx.x;           
    int sys   = blockIdx.x;            
    int row0  = blockIdx.y * TILE;     
    int row   = row0 + tid;            
    int base  = sys * N;               

    // load into shared memory
    if (row < N) {
        sa[tid] = a_in[base + row];
        sb[tid] = b_in[base + row];
        sc[tid] = c_in[base + row];
        sd[tid] = d_in[base + row];
    }
    __syncthreads();

    // CR forward
    for (int lvl = 0; lvl < CR_LEVELS; ++lvl) {
        int stride = 1 << lvl;
        if (row < N && (row & (2*stride - 1)) == 0 && row + stride < N && tid + stride < TILE)
        {
            float k = sa[tid + stride] / sb[tid];
            sb[tid] -= sc[tid] * k;
            sd[tid] -= sd[tid + stride] * k;
            sc[tid] = 0.0f;
            sa[tid + stride] = 0.0f;
        }
        __syncthreads();
    }

    // PCR forward
    for (int stride = 1<<CR_LEVELS; stride < N; stride <<= 1) {
        float old_sa = (row < N ? sa[tid] : 0.0f);
        float old_sb = (row < N ? sb[tid] : 1.0f);
        float old_sc = (row < N ? sc[tid] : 0.0f);
        float old_sd = (row < N ? sd[tid] : 0.0f);

        int iL = row - stride, iR = row + stride;

        float sbL=1, saL=0, scL=0, sdL=0;
        if (row < N && iL >= 0) {
            if ((iL >= row0) && (iL <  row0 + TILE)) {
                sbL = sb[tid - stride];
                saL = sa[tid - stride];
                scL = sc[tid - stride];
                sdL = sd[tid - stride];
            } else {
                sbL = b_in[base + iL];
                saL = a_in[base + iL];
                scL = c_in[base + iL];
                sdL = d_in[base + iL];
            }
        }

        float sbR=1, saR=0, scR=0, sdR=0;
        if (row < N && iR < N) {
            if ((iR >= row0) && (iR <  row0 + TILE)) {
                sbR = sb[tid + stride];
                saR = sa[tid + stride];
                scR = sc[tid + stride];
                sdR = sd[tid + stride];
            } else {
                sbR = b_in[base + iR];
                saR = a_in[base + iR];
                scR = c_in[base + iR];
                sdR = d_in[base + iR];
            }
        }

        float alpha = (row < N && iL >= 0) ? -old_sa / sbL : 0.0f;
        float gamma = (row < N && iR <  N) ? -old_sc / sbR : 0.0f;

        float nsa = alpha * saL;
        float nsc = gamma * scR;
        float nsb = old_sb + alpha * scL + gamma * saR;
        float nsd = old_sd + alpha * sdL + gamma * sdR;

        __syncthreads();
        if (row < N) {
            sa[tid] = nsa;
            sb[tid] = nsb;
            sc[tid] = nsc;
            sd[tid] = nsd;
        }
        __syncthreads();
    }

    // direct solve
    if (row < N) sx[tid] = sd[tid] / sb[tid];
    __syncthreads();

    // CR backward
    for (int lvl = CR_LEVELS-1; lvl >= 0; --lvl) {
        int stride = 1 << lvl;
        if (row < N
            && (row & (2*stride - 1)) == stride
            && tid >= stride
            && tid + stride < TILE)
        {
            float xim = sx[tid - stride];
            float xip = ((row + stride < N) ? sx[tid + stride] : 0.0f);
            sx[tid] = (sd[tid] - sa[tid]*xim - sc[tid]*xip) / sb[tid];
        }
        __syncthreads();
    }

    if (row < N) x_out[base + row] = sx[tid];
}


int main(int argc,char**argv)
{
    double t0_all = CycleTimer::currentSeconds();

    if (argc!=2) {
      fprintf(stderr,"Usage: %s <in.txt>\n",argv[0]);
      return 1;
    }
    std::ifstream fin(argv[1]);
    if (!fin) { perror("open"); return 1; }

    int M,N;  
    fin>>M>>N;
    if ((N & (N-1))!=0)
        fprintf(stderr,"Warning: N=%d not power-of-two, but it will still work.\n",N);

    std::vector<float> h_a(M*N), h_b(M*N),
                       h_c(M*N), h_d(M*N),
                       h_x(M*N,0.f);

    for (int system=0; system<M; ++system) {
        for (int i=0; i<N; ++i) fin>>h_b[system*N+i];
        h_a[system*N+0] = 0.f;
        for (int i=1; i<N; ++i) fin>>h_a[system*N+i];
        for (int i=0; i<N-1; ++i) fin>>h_c[system*N+i];
        h_c[system*N+N-1] = 0.f;
        for (int i=0; i<N; ++i) fin>>h_d[system*N+i];
    }

    float *d_a,*d_b,*d_c,*d_d,*d_x;
    size_t bytes = M*N*sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_a,bytes));
    CUDA_CHECK(cudaMalloc(&d_b,bytes));
    CUDA_CHECK(cudaMalloc(&d_c,bytes));
    CUDA_CHECK(cudaMalloc(&d_d,bytes));
    CUDA_CHECK(cudaMalloc(&d_x,bytes));

    CUDA_CHECK(cudaMemcpy(d_a,h_a.data(),bytes,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b,h_b.data(),bytes,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_c,h_c.data(),bytes,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_d,h_d.data(),bytes,cudaMemcpyHostToDevice));

    double t0 = CycleTimer::currentSeconds();
    dim3 block(TPB);
    dim3 grid(M,(N + TPB - 1)/TPB);
    size_t shmem = 5 * TPB * sizeof(float);

    
    hybrid_solve_multi_tiled<TPB><<<grid,block,shmem>>>
      (N, d_a,d_b,d_c,d_d, d_x);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    

    CUDA_CHECK(cudaMemcpy(h_x.data(),d_x,bytes,cudaMemcpyDeviceToHost));

    double t1 = CycleTimer::currentSeconds();
    std::cout<<"Computation time: "<<(t1-t0)<<" s\n";
    std::cout<<"Total time:       "<<(t1-t0_all)<<" s\n";


     // print result (comment out for large N)
    for (int system=0; system<M; ++system) {
        int base = system * N;
        for (int i=0; i<N; ++i) {
        if (i == 0) std::cout << "x" << system << " = ";
        std::cout << h_x[base+i] << " ";
        }
        std::cout << "\n";
    }



    // free device memory
    cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);
    cudaFree(d_d); cudaFree(d_x);
    return 0;
}
