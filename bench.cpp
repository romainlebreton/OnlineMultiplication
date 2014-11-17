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

void 
benchmarks (long length) {
    zz_p::init(268435459);
    zz_pX a, b, c;
    double u;
    long nc = 10, count;
    cout << length << "\t";

    /*------------------------*/
    /* Offline Multiplication */
    /*------------------------*/

    random(a, length);
    random(b, length);

    /* NTL multiplication */
    c.rep.SetLength(2*length-1);
    u = GetTime();
    for (count = 0; GetTime() < u+0.1; )
        for (long j = 0; j < nc; j++, count++)
            mul(c, a, b);
    u = GetTime() - u;
    double reftime = u/count;
    cout << (u/count) << "\t";

    /*-------------------------*/
    /* Relaxed Multiplications */
    /*-------------------------*/

    /* Naive multiplication */
    c.kill();
    c.rep.SetLength(length);
    u = GetTime();
    for (count = 0; GetTime() < u+0.1; )
        for (long j = 0; j < nc; j++, count ++)
            for(long i = 0; i < length; i++)
                lazy_mul(c, a, b, i);
    u = GetTime() - u;
    cout << (u / (count * reftime)) << "\t";

    /* [Fischer, Stockmeyer, 1974] algorithm */
    c.kill();
    c.rep.SetLength(2*length-1);
    u = GetTime();
    for (count = 0; GetTime() < u+0.1; )
        for (long j = 0; j < nc; j++, count ++)
            for(long i = 0; i < length; i++)
                relax_mul(c, a, b, i);
    u = GetTime() - u;
    cout << (u / (count * reftime)) << "\t";

    /* [Lebreton, Schost, 2014] algorithm */
    c.kill();
    c.rep.SetLength(2*length + 1);
    u = GetTime();
    for (count = 0; GetTime() < u+0.1; )
        for (long j = 0; j < nc; j++, count ++)
            for(long i = 0; i < length; i++)
                relax_middle(c, a, b, i);
    u = GetTime() - u;
    cout << (u / (count * reftime)) << "\t";

    /*------------------------------*/
    /* Semi-relaxed multiplications */
    /*------------------------------*/

    random(a, 2*length);
    random(b, length);

    /* [Fischer, Stockmeyer, 1974] algorithm */
    c.kill();
    c.rep.SetLength(3*length-2);
    u = GetTime();
    for (count = 0; GetTime() < u+0.1; )
        for (long j = 0; j < nc; j++, count++)
            for(long i = 0; i < length; i++)
                semi_relax(c, a, b, i);
    u = GetTime() - u;
    cout << (u / (count * reftime)) << "\t";

    /* [van der Hoeven, 2003] algorithm */
    c.kill();
    c.rep.SetLength(2*length);
    u = GetTime();
    for (count = 0; GetTime() < u+0.1; )
        for (long j = 0; j < nc; j++, count++)
            for(long i = 0; i < length; i++)
                semi_relax_middle(c, a, b, i);
    u = GetTime() - u;
    cout << (u / (count * reftime)) << "\t";

    cout << endl;
}

int main(){
    std::cout.precision(2);

    cout << "Length\tNTL mult (sec)\tlazy_mul/NTL\trelax_mul/NTL\trelax_middle/NTL\tsemi_relax/NTL\tsemi_relax_middle/NTL\n";
    cout << "------\t--------------\t------------\t-------------\t----------------\t--------------\t---------------------\n";
//    cout << "Length    NTL mult (sec)    Ratio Lazy mul    Ratio relax_mul    Ratio relax_middle   ";
//    cout << " Ratio semi_relax    Ratio semi_relax_middle \n";
    for (long pl= 8; pl < 12; pl++)
        for (int i = 0; i < 8; i++) {
            int j = (8+i) << (pl-3);
            benchmarks (j);
        }
}
