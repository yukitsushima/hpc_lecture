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
<<<<<<< HEAD
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
=======
  __m256 zero = _mm256_setzero_ps();
  for(int i=0; i<N; i+=8) {
    __m256 xi = _mm256_load_ps(x+i);
    __m256 yi = _mm256_load_ps(y+i);
    __m256 fxi = zero;
    __m256 fyi = zero;
    for(int j=0; j<N; j++) {
      __m256 dx = _mm256_set1_ps(x[j]);
      __m256 dy = _mm256_set1_ps(y[j]);
      __m256 mj = _mm256_set1_ps(m[j]);
      __m256 r2 = zero;
      dx = _mm256_sub_ps(xi, dx);
      dy = _mm256_sub_ps(yi, dy);
      r2 = _mm256_fmadd_ps(dx, dx, r2);
      r2 = _mm256_fmadd_ps(dy, dy, r2);
      __m256 mask = _mm256_cmp_ps(r2, zero, _CMP_GT_OQ);
      __m256 invR = _mm256_rsqrt_ps(r2);
      invR = _mm256_blendv_ps(zero, invR, mask);
      mj = _mm256_mul_ps(mj, invR);
      invR = _mm256_mul_ps(invR, invR);
      mj = _mm256_mul_ps(mj, invR);
      fxi = _mm256_fmadd_ps(dx, mj, fxi);
      fyi = _mm256_fmadd_ps(dy, mj, fyi);
    }
    _mm256_store_ps(fx+i, fxi);
    _mm256_store_ps(fy+i, fyi);
>>>>>>> upstream/master
  }
  for(int i=0; i<N; i++)
    printf("%d %g %g\n",i,fx[i],fy[i]);
}
