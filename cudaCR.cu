
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include "CycleTimer.h"
#include <iostream>

constexpr int TPB = 256;  // threads per block

#define CUDA_CHECK(call) do {                                 \
  cudaError_t err = call;                                     \
  if (err != cudaSuccess) {                                   \
    fprintf(stderr,                                       \
      "%s:%d error: %s\n",                             \
      __FILE__, __LINE__, cudaGetErrorString(err));          \
    exit(1);                                                  \
  }                                                           \
} while(0)

__device__
inline bool is_pivot(int i, int stride) {
  return (i % (2*stride)) == 0;
}

__global__
void CR_main_kernel(
    int N, int M, int stride,
    const float* __restrict__ a_in,
    const float* __restrict__ b_in,
    const float* __restrict__ c_in,
    const float* __restrict__ d_in,
    float* __restrict__ a_out,
    float* __restrict__ b_out,
    float* __restrict__ c_out,
    float* __restrict__ d_out)
{
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = M * N;
  if (gid >= total) return;
  int system = gid / N;
  int i = gid % N;
  int base = system * N;

  float ai = a_in[gid];
  float bi = b_in[gid];
  float ci = c_in[gid];
  float di = d_in[gid];

  if (is_pivot(i, stride)) {
    // get the coupling coefficients
    float k1 = (i >= stride) ? ai / b_in[base + i - stride] : 0.0f;
    float k2 = (i + stride < N) ? ci / b_in[base + i + stride]: 0.0f;
    float anew = (i >= stride) ? -a_in[base + i - stride] * k1 : 0.0f;
    float cnew = (i + stride < N) ? -c_in[base + i + stride] * k2 : 0.0f;
    float bnew = bi - (i >= stride ? c_in[base + i - stride] * k1 : 0.0f) - (i + stride < N ? a_in[base + i + stride] * k2 : 0.0f);
    float dnew = di - (i >= stride ? d_in[base + i - stride] * k1 : 0.0f) - (i + stride < N ? d_in[base + i + stride] * k2 : 0.0f);

    a_out[gid] = anew;
    b_out[gid] = bnew;
    c_out[gid] = cnew;
    d_out[gid] = dnew;
  } else {
    // if it's a non-pivot, just copy
    a_out[gid] = 0.0f;
    b_out[gid] = bi;
    c_out[gid] = 0.0f;
    d_out[gid] = di;
  }
}

__global__
void final_solve_kernel(
    int N, int M,
    const float* __restrict__ b_in,
    const float* __restrict__ d_in,
    float* __restrict__ x_out)
{
  // solve the last 1x1 systemtem
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = M * N;
  if (gid >= total) return;
  x_out[gid] = d_in[gid] / b_in[gid];
}

int main(int argc,char**argv) {
  double init = CycleTimer::currentSeconds();
  if (argc != 2) {
    fprintf(stderr,"Usage: %s <inputs/{input_name}.txt>\n",argv[0]);
    return 1;
  }
  std::ifstream fin(argv[1]);
  if (!fin) { perror("open"); return 1; }

  int M, N;
  fin >> M >> N;
  if ((N & (N-1)) != 0) {
    fprintf(stderr,"Error: N=%d not power-of-two\n", N);
    return 1;
  }

  size_t total = size_t(M)*N;
  size_t bytes = total * sizeof(float);

  std::vector<float> h_a(total), h_b(total), h_c(total), h_d(total), h_x(total,0.0f);

  // read input
  for (int system=0; system<M; ++system) {
    int base = system * N;
    for (int i=0; i<N; ++i) {
      fin >> h_b[base+i];
    }
    h_a[base+0] = 0.0f;
    for (int i=1; i<N; ++i){
      fin >> h_a[base+i];
    }
    for (int i=0; i<N-1; ++i){
      fin >> h_c[base+i];
    }
    h_c[base+N-1] = 0.0f;
    for (int i=0; i<N; ++i){
      fin >> h_d[base+i];
    }
  }

  // allocate device (double-buffer for CR)
  float *a1, *b1, *c1, *d1, *a2, *b2, *c2, *d2, *dx;

  CUDA_CHECK(cudaMalloc(&a1,bytes));
  CUDA_CHECK(cudaMalloc(&b1,bytes));
  CUDA_CHECK(cudaMalloc(&c1,bytes));
  CUDA_CHECK(cudaMalloc(&d1,bytes));
  CUDA_CHECK(cudaMalloc(&a2,bytes));
  CUDA_CHECK(cudaMalloc(&b2,bytes));
  CUDA_CHECK(cudaMalloc(&c2,bytes));
  CUDA_CHECK(cudaMalloc(&d2,bytes));
  CUDA_CHECK(cudaMalloc(&dx,bytes));

  // copy inputs to a1,b1,c1,d1
  CUDA_CHECK(cudaMemcpy(a1,h_a.data(),bytes,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(b1,h_b.data(),bytes,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(c1,h_c.data(),bytes,cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d1,h_d.data(),bytes,cudaMemcpyHostToDevice));

  double t0 = CycleTimer::currentSeconds();
  int numBlocks = (total + TPB-1)/TPB;
  dim3 grid(numBlocks), block(TPB);

  float *ain=a1;
  float *bin=b1;
  float *cin=c1;
  float *din=d1;
  float *aout=a2;
  float *bout=b2;
  float *cout=c2;
  float *dout=d2;

  int levels = 0;
  while ((1<<levels) < N) ++levels;


  // one kernel launch per CR level
  for (int lvl=0; lvl<levels; ++lvl) {
    int stride = 1<<lvl;
    CR_main_kernel
      <<<grid,block>>>(N,M,stride,
                      ain,bin,cin,din,
                      aout,bout,cout,dout);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    // swap buffers
    std::swap(ain,aout);
    std::swap(bin,bout);
    std::swap(cin,cout);
    std::swap(din,dout);
  }

  // final 1×1 solve
  final_solve_kernel
    <<<grid,block>>>(N,M, bin, din, dx);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  // copy result into h_x
  CUDA_CHECK(cudaMemcpy(h_x.data(),dx,bytes,cudaMemcpyDeviceToHost));
  double t1 = CycleTimer::currentSeconds();

  std::cout<<"CR time: "<<(t1-t0)<<" s\n";
  std::cout<<"Total program: "<<(CycleTimer::currentSeconds()-init)<<" s\n";

  // print result (comment out for large N)
  for (int system=0; system<M; ++system) {
    int base = system * N;
    for (int i=0; i<N; ++i) {
      if (i == 0){
        std::cout << "x" << system << " = ";
      }
      std::cout << h_x[base+i] << " ";
    }
    std::cout << "\n";
  }

  // free everything 
  cudaFree(a1);
  cudaFree(b1);
  cudaFree(c1);
  cudaFree(d1);
  cudaFree(a2);
  cudaFree(b2);
  cudaFree(c2);
  cudaFree(d2);
  cudaFree(dx);

  return 0;
}
