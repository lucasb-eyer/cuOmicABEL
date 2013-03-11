/*
 * Copyright (c) 2010-2012, Diego Fabregat-Traver and Paolo Bientinesi.
 * All rights reserved.
 *
 * This file is part of OmicABEL.
 * 
 * OmicABEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * OmicABEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with OmicABEL. If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * Translated from http://www.netlib.org/fmm/fmin.f by:
 *   Diego Fabregat-Traver (fabregat@aices.rwth-aachen.de)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "optimization.h"


static double dsign( double a, double b )
{
	return b >= 0.0 ? fabs(a) : -fabs(a);
}

double minimize( double ax, double bx, double tol, double (*f)(double, void *), void *data )
{
	
      double a, b, c, d, e,
			 eps, xm, p, q, r,
			 tol1, tol2, u, v, w;
      double fu, fv, fw, fx, x;

// c is the squared inverse of the golden ratio
      c = 0.5 * ( 3.0 - sqrt(5.0) );

// eps is approximately the square root of the relative machine
//  precision.
      eps = 1.0;
mach:   eps = eps / 2.0;
      tol1 = 1.0 + eps;
      if ( tol1 > 1.0) 
		  goto mach;
      eps = sqrt( eps );

//  initialization
      a = ax;
      b = bx;
      v = a + c*(b - a);
      w = v;
      x = v;
      e = 0.0;
      fx = f( x, data );
      fv = fx;
      fw = fx;

//  main loop starts here
loop: 
	  /*printf("Iteration - sigma: %f\n", *sigma);*/
	  xm = 0.5 * (a + b);
      tol1 = eps * fabs(x) + tol / 3.0;
      tol2 = 2.0 * tol1;

//  check stopping criterion
      if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) 
		  goto end;

// is golden-section necessary
      if (fabs(e) <= tol1) 
		  goto golden;

//  fit parabola
      r = (x - w)*(fx - fv);
      q = (x - v)*(fx - fw);
      p = (x - v)*q - (x - w)*r;
      q = 2.0*(q - r);
      if (q > 0.0) p = -p;
      q = fabs(q);
      r = e;
      e = d;

//  is parabola acceptable
      if (fabs(p) >= fabs(0.5*q*r)) 
		  goto golden;
      if (p <= q*(a - x)) 
		  goto golden;
      if (p >= q*(b - x)) 
		  goto golden;

//  a parabolic interpolation step
      d = p/q;
      u = x + d;

//  f must not be evaluated too close to ax or bx
      if ((u - a) < tol2) d = dsign(tol1, xm - x);
      if ((b - u) < tol2) d = dsign(tol1, xm - x);
      goto fnotclose;

//  a golden-section step
golden:   
	  if (x >= xm) e = a - x;
      if (x < xm) e = b - x;
      d = c*e;

//  f must not be evaluated too close to x
fnotclose:   if (fabs(d) >= tol1) u = x + d;
      if (fabs(d) < tol1) u = x + dsign(tol1, d);
      fu = f( u, data );

//  update  a, b, v, w, and x
      if (fu > fx) 
		  goto sixty;
      if (u >= x) 
		  a = x;
      if (u < x)
		  b = x;
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
      goto loop;
sixty:   if (u < x) a = u;
      if (u >= x) b = u;
      if (fu <= fw) 
		  goto sevty;
      if (w == x) 
		  goto sevty;
      if (fu <= fv) 
		  goto eigty;
      if (v == x) 
		  goto eigty;
      if (v == w) 
		  goto eigty;
      goto loop;
sevty:   v = w;
      fv = fw;
      w = u;
      fw = fu;
      goto loop;
eigty:   v = u;
      fv = fu;
      goto loop; // continue

//c  end of main loop
end:  // fmin = x;

      return x; // fmin;
}


