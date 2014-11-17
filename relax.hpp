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

#ifndef _relax_HPP_
#define _relax_HPP_

#include <NTL/lzz_pX.h>
#include "middle_product.hpp"
#include "short_product.hpp"

#define KARX (16)
NTL_CLIENT

void PlainAddMul (zz_p *xp, const zz_p *ap, const zz_p *bp, long length);

void semi_relax_middle(zz_pX &c, zz_pX& a, zz_pX& b, long i);
void semi_relax(zz_pX &c, zz_pX& a, zz_pX& b, long i);

void lazy_mul(zz_pX& c, const zz_pX& a, const zz_pX& b, long i);
void relax_mul(zz_pX& c, zz_pX& a, zz_pX& b, long i);
void relax_middle(zz_pX &c, zz_pX& a, zz_pX& b, long i);
  
#endif
