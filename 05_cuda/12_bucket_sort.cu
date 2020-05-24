#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void initialize(int *bucket){
  int i = threadIdx.x;
  bucket[i] = 0;
}

__global__ void bucket_add(int *key, int *bucket) {
  int i = threadIdx.x;
  int content = key[i];
  atomicAdd(&bucket[content],1);
}

__global__ void bucket_return(int *key, int num, int offset){
  int i = threadIdx.x;
  key[i+offset] = num;
}

int main() {
  int n = 50;
  int range = 5;
  int *key;
  cudaMallocManaged(&key,n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int *bucket;
  cudaMallocManaged(&bucket,range*sizeof(int));
  initialize<<<1,range>>>(bucket);
  //for (int i=0; i<range; i++) {
  //  bucket[i] = 0;
  //}
  cudaDeviceSynchronize();
  bucket_add<<<1,n>>>(key,bucket);
  cudaDeviceSynchronize();
  //for (int i=0; i<n; i++) {
  //  bucket[key[i]]++;
  //}
  int offset = 0;
  for (int i=0;i<range;i++){
    bucket_return<<<1,bucket[i]>>>(key,i,offset);
    offset += bucket[i];
  }
  cudaDeviceSynchronize();
  //for (int i=0, j=0; i<range; i++) {
  //  for (; bucket[i]>0; bucket[i]--) {
  //    key[j++] = i;
  //  }
  //}

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
