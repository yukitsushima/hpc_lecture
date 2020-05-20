#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    //fx[i] = fy[i] = 0;
  }
  __m256 mvec  = _mm256_load_ps(m);
  for(int i=0; i<N; i++) {
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 yi = _mm256_set1_ps(y[i]);
    __m256 xj = _mm256_load_ps(x);
    __m256 yj = _mm256_load_ps(y);
    __m256 rx = _mm256_sub_ps(xi,xj);
    __m256 ry = _mm256_sub_ps(yi,yj);
    __m256 r2 = _mm256_add_ps(_mm256_mul_ps(rx,rx),_mm256_mul_ps(ry,ry));
    //_CMP_NEQ_OQ is also OK?
    __m256 mask = _mm256_cmp_ps(r2,_mm256_setzero_ps(),_CMP_GT_OQ);
    //rinv = 1 / r
    __m256 rinv = _mm256_rsqrt_ps(r2);
    //dfx = m * rx / r^3
    __m256 dfx = _mm256_mul_ps(mvec,rx);
    dfx = _mm256_mul_ps(dfx,rinv);
    dfx = _mm256_mul_ps(dfx,rinv);
    dfx = _mm256_mul_ps(dfx,rinv);
    dfx = _mm256_blendv_ps(_mm256_setzero_ps(),dfx,mask);
    __m256 dfy = _mm256_mul_ps(mvec,ry);
    dfy = _mm256_mul_ps(dfy,rinv);
    dfy = _mm256_mul_ps(dfy,rinv);
    dfy = _mm256_mul_ps(dfy,rinv);
    dfy = _mm256_blendv_ps(_mm256_setzero_ps(),dfy,mask);
    //reduction
    __m256 red = _mm256_hadd_ps(dfx,dfy);
    red = _mm256_hadd_ps(red,red);
    float red_unpack[N];
    _mm256_store_ps(red_unpack,red);
    fx[i] -= red_unpack[0] + red_unpack[4];
    fy[i] -= red_unpack[1] + red_unpack[5];
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
