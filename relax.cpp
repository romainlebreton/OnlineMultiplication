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
#include "middle_product.hpp"
#include "short_product.hpp"

#define KARX (16)
NTL_CLIENT

/*---------------------------------------------------------------*/
/* Plain addition and multiplication                             */
/* x = x + a * b                                                 */
/* We use preconditionning as in plain product :                 */
/* NTL precomputes the inverse of b[i] to speed up               */
/* the modular multiplication by b[i] mod p                      */
/*---------------------------------------------------------------*/
inline void
PlainAddMul (zz_p *xp, const zz_p *ap, const zz_p *bp, long length){
    if (length == 0) return;
    long p = zz_p::modulus();
    double pinv = zz_p::ModulusInverse();
    for (long i = 0; i < length; i++) {
        long t1 = rep(bp[i]);
        // Precomputation of b[i]^-1
        mulmod_precon_t bpinv = PrepMulModPrecon(t1, p, pinv);
        zz_p *xp1 = xp+i;
        for (long j = 0; j < length; j++) {
            long t2;
            t2 = MulModPrecon(rep(ap[j]), t1, p, bpinv);
            xp1[j].LoopHole() = AddMod(t2, rep(xp1[j]), p);
        } }
}

/*-----------------------------------------------------------------------*/
/*       Semi-relaxed algorithms for Addition-Multiplication of (a,b)    */
/*-----------------------------------------------------------------------*/
/*  The following algorithms perform only the i-th step of the add-mul   */
/*  After steps 0 to n-1, we have c = (c + a * b) mod x^n                */
/*  The algorithms are relaxed only w.r.t. b                             */
/*  The polynomial a is considered known in advance                      */
/*-----------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*  Algorithm of                                                            */
/*  [van der Hoeven, Relaxed multiplication using the middle product, 2003]	*/
/*--------------------------------------------------------------------------*/
void semi_relax_middle(zz_pX &c, zz_pX& a, zz_pX& b, long i){
    long i1 = i, ell = 0;
    // ell is the largest power of two dividing i
    while ((i1 & 1) == 1){
        ell++;
        i1 >>= 1;
    }
    ell = 1<<ell;
    // ell is the size of the middle product
    // if ell is small, we use naive middle product
    if (ell <= KARX) {
        const zz_p *ap, *bp;
        zz_p *cp;
        ap = a.rep.elts();
        bp = b.rep.elts() + i - ell + 1;
        cp = c.rep.elts() + i;
        tPlainAddMul(cp, ap, bp, ell);
        return;
    }
    // otherwise we use general middle product (based on FFT or Karatsuba)
    zz_pX bB, tmp;
    bB.rep.SetLength(ell);
    for (long j=0; j<ell; j++)
        bB.rep[j].LoopHole() = rep(b.rep[j+i-ell+1]);
    tmp = middle_product(trunc(a, 2*ell-1), bB, ell);
    long l = deg(tmp)+1;
    for (long j = 0; j < l; j++)
        c.rep[j+i] += tmp.rep[j];
}

/*------------------------------------------------------------------*/
/* Algorithm of                                                     */
/* [Fischer, Stockmeyer, Fast on-line integer multiplication, 1974] */
/*------------------------------------------------------------------*/
void semi_relax(zz_pX &c, zz_pX& a, zz_pX& b, long i){
    long i1 = i, sA = 0, sB = i, ell = 1;
    long go = 0;
    zz_pX bA, bB;
    long l;
    a.rep.SetLength(max(a.rep.length(), i+1));
    b.rep.SetLength(max(b.rep.length(), i+1));
    zz_p *crep = c.rep.elts()+i;
    // Succession of val_2(i+1) plain products
    do{
        if (ell <= KARX) {    // Optimization if the size of the product is too small
            const zz_p *arep = a.rep.elts()+sA, *brep = b.rep.elts()+sB;
            PlainAddMul(crep, arep, brep, ell);
        } else {
            bA.rep.SetLength(ell);
            bB.rep.SetLength(ell);
            for (long j=0; j<ell; j++)
                bA.rep[j] = a.rep[j+sA];
            for (long j=0; j<ell; j++)
                bB.rep[j] = b.rep[j+sB];
            mul(bB,bA,bB);
            l = deg(bB)+1;
            for (long j = 0; j < l; j++)
                c.rep[j+i] += bB.rep[j];
        }
        sA += ell;
        sB -= ell;
        ell <<= 1;
        go = i1 & 1;
        i1 >>= 1;
    } while (go == 1);
}

/*-----------------------------------------------------------------------*/
/*            Relaxed algorithms for Addition-Multiplication of (a,b)    */
/*-----------------------------------------------------------------------*/
/*  The following algorithms perform only the i-th step of the add-mul   */
/*  After steps 0 to n-1, we have c = (c + a * b) mod x^n                */
/*  The algorithms are relaxed w.r.t. a and b                            */
/*-----------------------------------------------------------------------*/

/*---------------------------------------------------------------*/
/*  Plain multiplication, which is online                        */
/*---------------------------------------------------------------*/
void
lazy_mul(zz_pX& c, const zz_pX& a, const zz_pX& b, long i){
    for (long k = 0; k <= i; k++)
        c.rep[i] += a.rep[k] * b.rep[i-k];
}

/*------------------------------------------------------------------*/
/* Algorithm of                                                     */
/* [Fischer, Stockmeyer, Fast on-line integer multiplication, 1974] */
/*------------------------------------------------------------------*/
void relax_mul(zz_pX& c, zz_pX& a, zz_pX& b, long i){
    if(i == 0){
        c.rep[0] += a.rep[0]*b.rep[0];
        return;
    }
    a.rep.SetLength(max(a.rep.length(), i+1));
    b.rep.SetLength(max(b.rep.length(), i+1));
    zz_p *crep = c.rep.elts()+i;
    // Optimization : We do the 1x1 products separately
    crep[0] += a.rep[0]*b.rep[i] + a.rep[i]*b.rep[0];
    zz_pX bA, bB, tmp;
    long k, l;
    long length;
    long u_length, k_length;
    length = 2;
    k = i + 2;
    // Succession of val_2(i+2) steps of plain products
    while((k&1) == 0){
        k = k >> 1;
        u_length = length-1;
        k_length = i-u_length;
        const zz_p *arep = a.rep.elts()+u_length, *brep = b.rep.elts()+k_length;
        // First of the possibly two products
        if (length <= KARX)    // Optimization if the size of the product is too small
            PlainAddMul(crep, arep, brep, length);
        else {
            bA.rep.SetLength(length);
            bB.rep.SetLength(length);
            zz_p *bArep = bA.rep.elts(), *bBrep = bB.rep.elts();
            for (long j=0; j<length; j++){
                bArep[j].LoopHole() = rep(arep[j]);
                bBrep[j].LoopHole() = rep(brep[j]);
            }
            mul(tmp,bA,bB);
            l = tmp.rep.length();
            const zz_p *tmprep = tmp.rep.elts();
            for (long j = 0; j < l; j++)
                crep[j] += tmprep[j];
        }
        // if k != 2, we have to do a second product
        if (k == 2)
            break;
        // Second product
        arep = a.rep.elts()+k_length;
        brep = b.rep.elts()+u_length;
        if (length <= KARX)    // Optimization if the size of the product is too small
            PlainAddMul(crep, arep, brep, length);
        else{
            zz_p *bArep = bA.rep.elts(), *bBrep = bB.rep.elts();
            for (long j=0; j<length; j++){
                bArep[j].LoopHole() = rep(arep[j]);
                bBrep[j].LoopHole() = rep(brep[j]);
            }
            mul(tmp,bA,bB);
            l = tmp.rep.length();
            const zz_p *tmprep = tmp.rep.elts();
            for (long j = 0; j < l; j++)
                crep[j] += tmprep[j];
        }
        length <<= 1;
    }
}

/*-------------------------------------------------------------------------*/
/* Algorithm from                                                          */
/* [Lebreton, Schost, A simple and fast online power series multiplication */
/*                    and its analysis]                                    */
/*-------------------------------------------------------------------------*/
void relax_middle(zz_pX &c, zz_pX& a, zz_pX& b, long i){
    if (i < 2) {
        lazy_mul(c,a,b,i);
        return;
    }
    long i1 = i+2, ell = 0;
    // ell is the largest power of two dividing i + 2
    while ((i1 & 1) == 0){
        ell++;
        i1 >>= 1;
    }
    ell = 1<<ell;
    // When i+2 is a power of two, we only perform a short product
    if (i1 == 1) {
        mul_highshort (c, a, b, i);
        return;
    }
    // Otherwise, we have two middle products to compute
    if (ell <= KARX) {    // Optimization if the size of the product is too small
        const zz_p *ap, *bp;
        zz_p *cp;
        ap = a.rep.elts();
        bp = b.rep.elts() + i - ell + 1;
        cp = c.rep.elts() + i;
        tPlainAddMul(cp, ap, bp, ell);
        ap = a.rep.elts() + i - ell + 1;
        bp = b.rep.elts();
        tPlainAddMul(cp, bp, ap, ell);
    } else {
        zz_pX bB, tmp;
        bB.rep.SetLength(ell);
        for (long j=0; j<ell; j++)
            bB.rep[j] = b.rep[j+i-ell+1];
        tmp  = middle_product(trunc(a, 2*ell-1), bB, ell);
        for (long j=0; j<ell; j++)
            bB.rep[j] = a.rep[j+i-ell+1];
        tmp += middle_product(trunc(b, 2*ell-1), bB, ell);
        long l = deg(tmp)+1;
        for (long j = 0; j < l; j++)
            c.rep[j+i] += tmp.rep[j];
    }
}
