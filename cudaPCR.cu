
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include "CycleTimer.h"
#include <iostream>

constexpr int TPB = 256; // threads/block

#define CUDA_CHECK(call) do {                                      \
  cudaError_t err = call;                                          \
  if (err != cudaSuccess) {                                        \
    fprintf(stderr, "%s:%d error: %s\n",                           \
      __FILE__, __LINE__, cudaGetErrorString(err));                \
    exit(1);                                                       \
  }                                                                \
} while(0)

__global__ void pcr_stride_kernel(
    int N, int M, int stride,
    const float* __restrict__ a_in,
    const float* __restrict__ b_in,
    const float* __restrict__ c_in,
    const float* __restrict__ d_in,
    float*       __restrict__ a_out,
    float*       __restrict__ b_out,
    float*       __restrict__ c_out,
    float*       __restrict__ d_out)
{
    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= M * N) return;
    int system = gid / N;
    int i   = gid % N;
    int base = system * N;

    float ai = a_in[base + i];
    float bi = b_in[base + i];
    float ci = c_in[base + i];
    float di = d_in[base + i];

    float alpha = 0.f, gamma = 0.f;
    if (i >= stride) {
        float bL = b_in[base + i - stride];
        alpha = - ai / bL;
    }
    if (i + stride < N) {
        float bR = b_in[base + i + stride];
        gamma = - ci / bR;
    }

    float anew = (i >= stride) ? alpha * a_in[base + i - stride] : 0.f;
    float cnew = (i + stride < N) ? gamma * c_in[base + i + stride] : 0.f;
    float bnew = bi + (i >= stride ? alpha * c_in[base + i - stride] : 0.f) + (i + stride < N ? gamma * a_in[base + i + stride] : 0.f);
    float dnew = di + (i >= stride ? alpha * d_in[base + i - stride] : 0.f) + (i + stride < N ? gamma * d_in[base + i + stride] : 0.f);

    a_out[base + i] = anew;
    b_out[base + i] = bnew;
    c_out[base + i] = cnew;
    d_out[base + i] = dnew;
}

__global__ void final_solve_kernel(
    int N, int M,
    const float* __restrict__ b_in,
    const float* __restrict__ d_in,
    float* __restrict__ x_out)
{
    int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= M * N) return;
    x_out[gid] = d_in[gid] / b_in[gid];
}

int main(int argc, char** argv) {
    double t0_all = CycleTimer::currentSeconds();
    if (argc != 2) {
      fprintf(stderr, "Usage: %s <testcase.txt>\n", argv[0]);
      return 1;
    }
    std::ifstream fin(argv[1]);
    if (!fin) { perror("open"); return 1; }

    int M, N;
    fin >> M >> N;
    if ((N & (N - 1)) != 0) {
      fprintf(stderr, "Warning: N=%d is not power-of-two. PCR still works but likely suboptimal.\n", N);
    }

    size_t total = size_t(M) * N;
    size_t bytes = total * sizeof(float);

    std::vector<float> h_a(total), h_b(total),
                       h_c(total), h_d(total),
                       h_x(total, 0.f);

    for (int system = 0; system < M; ++system) {
      int base = system * N;
      for (int i = 0; i < N; ++i) fin >> h_b[base + i]; // main diagonal
      h_a[base + 0] = 0.f;
      for (int i = 1; i < N; ++i) fin >> h_a[base + i]; // sub-diagonal
      for (int i = 0; i < N-1; ++i) fin >> h_c[base + i]; // super-diagonal
      h_c[base + N-1] = 0.f;
      for (int i = 0; i < N; ++i) fin >> h_d[base + i]; // right-hand side
    }

    float *d_a1, *d_b1, *d_c1, *d_d1,
          *d_a2, *d_b2, *d_c2, *d_d2,
          *d_x;
    CUDA_CHECK(cudaMalloc(&d_a1, bytes));
    CUDA_CHECK(cudaMalloc(&d_b1, bytes));
    CUDA_CHECK(cudaMalloc(&d_c1, bytes));
    CUDA_CHECK(cudaMalloc(&d_d1, bytes));
    CUDA_CHECK(cudaMalloc(&d_a2, bytes));
    CUDA_CHECK(cudaMalloc(&d_b2, bytes));
    CUDA_CHECK(cudaMalloc(&d_c2, bytes));
    CUDA_CHECK(cudaMalloc(&d_d2, bytes));
    CUDA_CHECK(cudaMalloc(&d_x , bytes));

    CUDA_CHECK(cudaMemcpy(d_a1, h_a.data(), bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b1, h_b.data(), bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_c1, h_c.data(), bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_d1, h_d.data(), bytes, cudaMemcpyHostToDevice));

    double t0 = CycleTimer::currentSeconds();
    int numBlocks = (total + TPB - 1) / TPB;
    dim3 grid(numBlocks), block(TPB);

    // PCR forward‐reduction - we use one kernel per stride
    float *a_in=d_a1, *b_in=d_b1, *c_in=d_c1, *d_in=d_d1;
    float *a_out=d_a2,*b_out=d_b2,*c_out=d_c2,*d_out=d_d2;
    
    for (int stride = 1; stride < N; stride <<= 1) {
      pcr_stride_kernel
        <<<grid, block>>>(N, M, stride,
                         a_in, b_in, c_in, d_in,
                         a_out,b_out,c_out,d_out);
      CUDA_CHECK(cudaGetLastError());
      CUDA_CHECK(cudaDeviceSynchronize());
      std::swap(a_in, a_out);
      std::swap(b_in, b_out);
      std::swap(c_in, c_out);
      std::swap(d_in, d_out);
    }

    // final solve
    final_solve_kernel
      <<<grid, block>>>(N, M, b_in, d_in, d_x);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    

    // copy back
    CUDA_CHECK(cudaMemcpy(h_x.data(), d_x, bytes, cudaMemcpyDeviceToHost));
    double t1 = CycleTimer::currentSeconds();

    std::cout<<"PCR time:        "<<(t1-t0)   <<" s\n";
    std::cout<<"Total program:   "<<(CycleTimer::currentSeconds()-t0_all)<<" s\n";

    // print result (comment out for large N)
    for (int system=0; system<M; ++system) {
      int base = system * N;
      for (int i=0; i<N; ++i) {
      if (i == 0) std::cout << "x" << system << " = ";
      std::cout << h_x[base+i] << " ";
      }
      std::cout << "\n";
    }


    // cleanup
    cudaFree(d_a1); cudaFree(d_b1); cudaFree(d_c1); cudaFree(d_d1);
    cudaFree(d_a2); cudaFree(d_b2); cudaFree(d_c2); cudaFree(d_d2);
    cudaFree(d_x);

    return 0;
}
