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

#ifndef __SHORT_PRODUCT_HPP
#define __SHORT_PRODUCT_HPP

#include <NTL/lzz_pX.h>

NTL_CLIENT

void mul_short(zz_pX& c, const zz_pX& a, const zz_pX& b, long s);
void mul_short(zz_pX& c, const zz_pX& a, const zz_pX& b);
void mul_highshort (zz_pX& c, const zz_pX& a, const zz_pX& b, long step);
#endif
