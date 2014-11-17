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

#include "relax.hpp"


void test(const zz_pX& a,const zz_pX& b, long n) {
    if (IsZero(trunc(a-b,n)))
        cout << "\t\t true";
    else
        cout << "\t\t false";
}

void
check (long length) {
    zz_p::init(268435459); // First prime after 2^28
    zz_pX a, b, c;
    cout << length;

    /*-------------------------*/
    /* Relaxed Multiplications */
    /*-------------------------*/
    random(a, length);
    random(b, length);

    /* Naive multiplication */
    c.rep.SetLength(2*length-1);
    for(long i = 0; i < length; i++)
        lazy_mul(c, a, b, i);
    test(c,a*b,length);

    /* [Fischer, Stockmeyer, 1974] algorithm */
    c.kill();
    c.rep.SetLength(2*length-1);
    for(long i = 0; i < length; i++)
        relax_mul(c, a, b, i);
    test(c,a*b,length);

    /* [Lebreton, Schost, 2014] algorithm */
    c.kill();
    c.rep.SetLength(2*length + 1);
    for(long i = 0; i < length; i++)
        relax_middle(c, a, b, i);
    test(c,a*b,length);


    /*------------------------------*/
    /* Semi-relaxed multiplications */
    /*------------------------------*/
    random(a, 2*length);
    random(b, length);

    /* [Fischer, Stockmeyer, 1974] algorithm */
    c.kill();
    c.rep.SetLength(3*length-2);
    for(long i = 0; i < length; i++)
        semi_relax(c, a, b, i);
    test(c,a*b,length);

    /* [van der Hoeven, 2003] algorithm */
    c.kill();
    c.rep.SetLength(2*length);
    for(long i = 0; i < length; i++)
        semi_relax_middle(c, a, b, i);
    test(c,a*b,length);

    cout << endl;
}


int main(){
    std::cout.precision(2);

    cout << "\n\nCheck correctness of relaxed algorithms\n";
    cout << "Length\t\t Lazy mul\t relax_mul\t relax_middle\t semi_relax\t semi_relax_middle\n";
    cout << "------\t\t --------\t ---------\t ------------\t ----------\t -----------------\n";
    for (long l = 1; l < 10; l++)
        check(l);
}
