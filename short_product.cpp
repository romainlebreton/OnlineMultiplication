/*
Copyright (c) 2014 Romain Lebreton and Éric Schost.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; see file COPYING. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <NTL/lzz_pX.h>

#define KARX 16

NTL_CLIENT

/*--------------------------------------*/
/* Algorithms for the low short product */
/* c = a*b mod x^s                      */
/*--------------------------------------*/

/*---------------------------------------------------------------*/
/* plain short product, using the same kind of preconditionning  */
/* as in plain product                                           */
/*---------------------------------------------------------------*/
void plainmul_short(zz_p *c, const zz_p *a, const zz_p *b, long s){

  long p = zz_p::modulus();
  double pinv = zz_p::ModulusInverse();

  // first column
  long t1 = rep(a[0]);
  mulmod_precon_t bpinv = PrepMulModPrecon(t1, p, pinv); 
  for (long j = 0; j < s; j++)
    c[j].LoopHole() = MulModPrecon(rep(b[j]), t1, p, bpinv);

  // all others, done with add-muls.
  for (long i = 1; i < s; i++){
    long t1 = rep(a[i]);
    mulmod_precon_t bpinv = PrepMulModPrecon(t1, p, pinv);
    for (long j = i, l = 0; j < s; j++, l++){
      long t2 = MulModPrecon(rep(b[l]), t1, p, bpinv);
      c[j].LoopHole() = AddMod(t2, rep(c[j]), p);
    }
  }
}

/*---------------------------------------------------------------*/
/* karatsuba short product, using Hanrot-Zimmermann's "variant"  */
/*---------------------------------------------------------------*/
void karmul_short(zz_p *c, const zz_p *a, const zz_p *b, long s, zz_p *wk){

  long k, l, h;
  if (s < KARX){
    plainmul_short(c, a, b, s);
    return;
  }

  h = (s+1) >> 1;
  l = s >> 1;
  k = l << 1;

  // decomposes a
  zz_p *G0 = wk;
  zz_p *G1 = wk + h*2;
  zz_p *G2 = wk + h*4;
  for(long i = 0, j = 0; i < s-1; i+=2, j++){
    G0[j] = a[i];
    G1[j] = a[i]+a[i+1];
    G2[j] = a[i+1];
  }

  // decomposes b
  zz_p *GG0 = wk + h;
  zz_p *GG1 = wk + h*3;
  zz_p *GG2 = wk + h*5;
  for(long i = 0, j = 0; i < s-1; i+=2, j++){
    GG0[j] = b[i];
    GG1[j] = b[i]+b[i+1];
    GG2[j] = b[i+1];
  }

  // fix when s is odd
  if (k < s){
    G0[l] = a[k];
    GG0[l] = b[k];
  }

  zz_p *P0 = wk + h*6;
  zz_p *P1 = wk;
  zz_p *P2 = wk + h;
  wk += 7*h;

  // 3 recursive calls
  karmul_short(P0, G0, GG0, h, wk);

  if (k < s)
    h--;
  karmul_short(P1, G1, GG1, h, wk);
  karmul_short(P2, G2, GG2, h, wk);

  // sets the results
  c[0] = P0[0];
  c[1] = P1[0] - P0[0] - P2[0];
  for(long i = 2, j = 1; i < s-1; i+=2, j++){
    c[i]   = P0[j] + P2[j-1];
    c[i+1] = P1[j] - P0[j] - P2[j];
  }

  // fix when s is odd
  if (k < s)
    c[k] = P0[l] + P2[l-1];
}

/*---------------------------------------------------------------*/
/* first does Karatsuba, then plain mul in low degree            */
/*---------------------------------------------------------------*/
void mul_short(zz_pX& c, const zz_pX& a, const zz_pX& b, long s){
  c.rep.SetLength(s);
  if (s+1 <= NTL_zz_pX_MUL_CROSSOVER) {
    vec_zz_p wk;
    wk.SetLength(7*2*s); // Not Optimal
    karmul_short(c.rep.elts(), a.rep.elts(), b.rep.elts(), s, wk.elts());
    return;
  }
  MulTrunc(c, a, b, s + 1);
}

void mul_short(zz_pX& c, const zz_pX& a, const zz_pX& b) {
  long s = max (deg(a), deg(b));
  if (deg(a) <= NTL_zz_pX_MUL_CROSSOVER) {
    mul_short(c, a, b, s + 1);
    return;
  }
  MulTrunc(c, a, b, s + 1);
}


/*----------------------------------------------------*/
/* Algorithms for the high short product and addition */
/* c = c + a*b div x^(s-1)                            */
/* This short product is required in relax_middle     */
/* It is computed as rev( short_mul (rev(a) * rev(b)))*/
/*----------------------------------------------------*/

/*-------------------------------------------------*/
/* Naive algorithm                                 */
/* We use preconditionning as in plain product :   */
/* NTL precomputes the inverse of b[i] to speed up */
/* the modular multiplication by b[i] mod p        */
/*-------------------------------------------------*/
void
PlainHighShortAddMul (zz_p *xp, const zz_p *ap, const zz_p *bp, long s){
  if (s == 0) return;
  long p = zz_p::modulus();
  double pinv = zz_p::ModulusInverse();
  for (long i = 0; i < s; i++) {
    long t1 = rep(bp[i]);
    // Precomputation of b[i]^-1
    mulmod_precon_t bpinv = PrepMulModPrecon(t1, p, pinv); 
    zz_p *xp1 = xp+i;
    for (long j = s - 1 - i; j < s; j++) {
      long t2;
      t2 = MulModPrecon(rep(ap[j]), t1, p, bpinv);
      xp1[j].LoopHole() = AddMod(t2, rep(xp1[j]), p);
    } }
}

/*------------------------------*/
/* karatsuba high short product */
/* by reversing karmul_short    */
/*------------------------------*/
void karmul_highshort(zz_p *c, const zz_p *a, const zz_p *b, long s){  
  long k, l, h;
  if (s < KARX){
    PlainHighShortAddMul (c, a, b, s);
    return;
  }
  vec_zz_p wkv;
  wkv.SetLength(7*2*s); // Not optimal
  zz_p *wk = wkv.elts();
  h = (s+1) >> 1;
  l = s >> 1;
  k = l << 1;

  // decomposes a
  zz_p *G0 = wk;
  zz_p *G1 = wk + h*2;
  zz_p *G2 = wk + h*4;
  for(long j = 0, i = s-1; i > 0; j++, i-=2){
    G0[j] = a[i];
    G1[j] = a[i]+a[i-1];
    G2[j] = a[i-1];
  }

  // decomposes b
  zz_p *GG0 = wk + h;
  zz_p *GG1 = wk + h*3;
  zz_p *GG2 = wk + h*5;
  for(long j = 0, i = s-1; i > 0; j++, i-=2){
    GG0[j] = b[i];
    GG1[j] = b[i]+b[i-1];
    GG2[j] = b[i-1];
  }

  // fix when s is odd
  if (k < s){
    G0[l]  = a[0];
    GG0[l] = b[0];
  }

  zz_p *P0 = wk + h*6;
  zz_p *P1 = wk;
  zz_p *P2 = wk + h;
  wk += 7*h;

  // 3 recursive calls
  karmul_short(P0, G0, GG0, h, wk);

  if (k < s)
    h--;
  karmul_short(P1, G1, GG1, h, wk);
  karmul_short(P2, G2, GG2, h, wk);

  // sets the results
  long top = 2*s-2;
  c[top] = P0[0];
  c[top-1] = P1[0] - P0[0] - P2[0];
  for(long j = 1, i = top - 2; i > s-1; j++, i-=2){
    c[i]   = P0[j] + P2[j-1];
    c[i-1] = P1[j] - P0[j] - P2[j];
  }

  // fix when s is odd
  if (k < s)
    c[top-k] = P0[l] + P2[l-1];
}

/*------------------------------------------------------------*/
/* FFT 'false' short product :                                */
/* Let 2^k be the next power of two >= max(deg(a),deg(b)) + 1 */
/* Let q = (a * b) div x^2^k, deg(q) < 2^k                    */
/* Let r = (a * b) mod x^2^k, deg(r) < 2^k                    */
/* Output : c = (a * b) mod x^2^k    if p == q                */
/*          c = (a * b) div x^2^k    if p == r                */
/* We use (a * b) mod (x^2^k - 1) = q + r                     */
/*                                                            */
/* Remark : The output is the high short product when         */
/* max(deg(a),deg(b)) + 2 is a power of two and p == r        */
/*------------------------------------------------------------*/
void FFTMul_short(zz_pX& x, const zz_pX& a, const zz_pX& b, const zz_pX& p) {

  long k, d;
  
  if (IsZero(a) || IsZero(b)) {
    clear(x);
    return;
  }
  
  d = max (deg(a), deg(b));
  k = NextPowerOfTwo(d+1);
  
  fftRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);
  
  TofftRep(R1, a, k);
  TofftRep(R2, b, k);
  mul(R1, R1, R2);
  FromfftRep(x, R1, 0, 1 << k);
  x -= p;
}


void mul_highshort(zz_pX& c, const zz_pX& a, const zz_pX& b, long s) {
  if (s < NTL_zz_pX_MUL_CROSSOVER) {
    karmul_highshort (c.rep.elts(), a.rep.elts(), b.rep.elts(), s+1);
    return; 
  }
  zz_pX ap, bp, rp, tp;
  ap.rep.SetLength(s+2);
  for (long i = 0; i < s + 1; i ++)
    ap.rep[i+1] = a.rep[i];
  bp.rep.SetLength(s+2);
  for (long i = 0; i < s + 1; i ++)
    bp.rep[i+1] = b.rep[i];
  rp.rep.SetLength(s+2);
  for (long i = 0; i < s; i ++)
    rp.rep[i+2] = c.rep[i];
  FFTMul_short(tp, ap, bp, rp);    
  for (long i = s, j = 0; i <= 2*s; i++, j++)
    c.rep[i] = tp.rep[j];
  return;
}
