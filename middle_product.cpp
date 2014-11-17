/*
Copyright (c) 2014 Romain Lebreton, Grégoire Lecerf and Éric Schost.

Part of this code comes from
G. Lecerf and É. Schost. Experimental Tellegen support for NTL, 2003.
see http://lecerf.perso.math.cnrs.fr/software/tellegen/index.html

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
#include <assert.h>


#define KARX (32)

NTL_CLIENT

/*----------------------------------------*/
/* reverse in degree < d                  */
/*----------------------------------------*/
zz_pX rev(const zz_pX& a, long d){
  zz_pX rA;
  if (deg(a) >= d)
    Error("degree too large to reverse");
  rA.rep.SetLength(d);
  for (long i = 0; i < d; i++)
    rA.rep[i] = coeff(a, d-1-i);
  rA.normalize();
  return rA;
}

/*------------------------------------------*/
/* naive transposed product of (a,b)        */
/* returns x=(a*b div x^(N-1)) mod x^N      */
/*------------------------------------------*/
void tPlainMul(zz_p *xp, const zz_p *ap, const zz_p *bp, long N){
   if (N == 0) return;

   for (long i = 0; i < N; i++)
     clear(xp[i]);

   long p = zz_p::modulus();
   double pinv = zz_p::ModulusInverse();

   for (long i = 0; i < N; i++) {
      long t1 = rep(bp[i]);
      mulmod_precon_t bpinv = PrepMulModPrecon(t1, p, pinv); // ((double) t1)*pinv;
      for (long j = 0; j < N; j++) {
         long t2 = MulModPrecon(rep(ap[i+j]), t1, p, bpinv);
         xp[j].LoopHole() = AddMod(t2, rep(xp[j]), p);
      }
   }
}


/*------------------------------------------*/
/* Karatsuba's transposed product of (a,b)  */
/* returns x=(a*b div x^(N-1)) mod x^N      */
/*------------------------------------------*/
void KarSub1(zz_p *T, const zz_p *b, long sb){
   long i;

   for (i = 0; i < sb; i++)
      sub(T[i], T[i], b[i]);
}

void KarFoldSub1(zz_p *T, const zz_p *b, long sb, long hsa){
   long m = sb - hsa;
   long i;

   for (i = 0; i < m; i++)
      sub(T[i], b[hsa+i], b[i]);

   for (i = m; i < hsa; i++)
     NTL::negate(T[i], b[i]);
}

void KarAdd2(zz_p *T, const zz_p *a, const zz_p *b, long s){
   long i;
   for (i = 0; i < s; i++)
      add(T[i], a[i], b[i]);
}

void KarAdd1(zz_p *T, const zz_p *b, long sb){
   long i;

   for (i = 0; i < sb; i++)
      add(T[i], T[i], b[i]);
}


// sa=sb
void TKarMul(zz_p *b, const zz_p *c, long sc, const zz_p *a, long s, zz_p *stk){
  if (s < KARX){
    tPlainMul(b, c, a, s);
    return;
  }

  long hs = (s + 1) >> 1;
  long hs2 = hs << 1;
  long hs3 = hs2 + hs;

  zz_p *T = stk;
  stk += hs3;
  zz_p *T1 = T+hs;
  zz_p *T2 = T1+hs;

  KarFoldSub1(T, a, s, hs);
  TKarMul(b, c+hs, min(hs2-1,sc-hs), T, hs, stk);

  long i = sc - hs;
  KarAdd2(T, c + hs, c, i);
  for(; i < min(sc, s+s-1-hs); ++i) {T[i] = c[i]; cout << "1" <<endl;}
  for(; i < s+s-1-hs; ++i) {clear(T[i]); cout << "2" <<endl;}
  TKarMul(b+hs, T1, s+s-1-hs2, a+hs, s-hs, stk);
  // cout << b[hs+1] << endl;
  KarSub1(b+hs, b, s-hs);

  // cout << b[hs+1] << endl;

  TKarMul(T2, T, hs2-1, a, hs, stk);
  KarAdd1(b, T2, hs);
}



/*----------------------------------------*/
/* middle product via FFT                 */
/* returns x=(a*b div x^(N-1)) mod x^N    */
/*----------------------------------------*/
void middle_FFT(zz_pX& x, const zz_pX& a, const zz_pX& b, long N){
  long k = NextPowerOfTwo(2*N-1);
  fftRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);
  TofftRep(R1, a, k);
  TofftRep(R2, b, k);
  mul(R1, R1, R2);
  FromfftRep(x, R1, N-1, 2*N-2);
}


/*--------------------------------------*/
/* middle product of (a,b)              */
/*   a has length <= 2*N-1              */
/*   b has length <= N                  */
/*   x has length <= N                  */
/* returns x=(a*b div x^(N-1)) mod x^N  */
/*--------------------------------------*/
zz_pX middle_product(const zz_pX& c, const zz_pX& a, long N){
  
  zz_pX b;

  if (c == 0 || a == 0){
    clear(b);
    return b;
  }

  if (deg(c) >= 2*N-1 || deg(a) >= N)
    Error("degree mismatch for middle product");

  if (N > NTL_zz_pX_MUL_CROSSOVER){
    middle_FFT(b, a, c, N);
    return b;
  }

  zz_pX ra = rev(a, N);
  const zz_p *ap, *cp;
  zz_p *bp;
  
  ap = ra.rep.elts();
  cp = c.rep.elts();
  
  b.rep.SetLength(N);
  bp = b.rep.elts();
  
  if (N < KARX)
    tPlainMul(bp, cp, ap, N);
  else {
    long n, hn, sp, depth;
    
    n = N;
    sp = 0;
    depth = 0;
    do {
      hn = (n+1) >> 1;
      sp += (hn << 2) - 1;
      n = hn;
      depth++;
    } while (n >= KARX);
    vec_zz_p stk;
    stk.SetLength(sp);
    TKarMul(bp, cp, 2*N-1, ap, N, stk.elts());
  }
  b.normalize();
  return b;
}
