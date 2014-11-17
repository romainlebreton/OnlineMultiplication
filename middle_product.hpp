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

#ifndef _middle_product__HPP_
#define _middle_product__HPP_

#include <NTL/lzz_pX.h>

NTL_CLIENT

/*--------------*/
/* Declarations */
/*--------------*/
zz_pX middle_product(const zz_pX& a, const zz_pX& b, long N);
void TKarMul(zz_pX& b, const zz_pX& c, const zz_pX& a, long N);

/*---------------------------------------------------*/
/* naive transposed addition multiplication of (a,b) */
/* returns x = x + (a*b div x^(N-1)) mod x^N         */
/*---------------------------------------------------*/
inline void
tPlainAddMul (zz_p *xp, const zz_p *ap, const zz_p *bp, long N){
  if (N == 0) return;
  long p = zz_p::modulus();
  double pinv = zz_p::ModulusInverse();
  for (long i = 0; i < N; i++) {
    long t1 = rep(bp[N - 1 - i]);                          // reverse polynomial
    mulmod_precon_t bpinv = PrepMulModPrecon(t1, p, pinv); // ((double) t1)*pinv;
    for (long j = 0; j < N; j++) {
      long t2 = MulModPrecon(rep(ap[i + j]), t1, p, bpinv);
      xp[j].LoopHole() = AddMod(t2, rep(xp[j]), p);
    } }
};

#endif
