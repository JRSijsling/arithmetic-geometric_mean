"""
 *  Calculation of an approximation of the period matrix
 *
 *  Copyright (C) 2016  Nicolas Mascot (n.a.v.mascot@warwick.ac.uk),
 *                      J.R. Sijsling (sijsling@gmail.com)
 *
 *  Distributed under the terms of the GNU General License (GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc., 51
 *  Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

class IntegrationApproximation:
    def __init__(self, f, alphas, prec_approx = 20, mesh = 2):
        self.f = f
        self.alphas = alphas
        # Derived objects:
        self.g = (self.f.degree() - 1) // 2
        self.prec_approx = prec_approx
        self.mesh = mesh
        self.CC_approx = ComplexField(self.prec_approx * log(10) / log(2))
        self.alphas_approx = [ self.CC_approx(alpha) for alpha in alphas ]
        # For caching results:
        self.period_dict = { }

        # Loading Nicolas' code for integration:

        # Branch values:
        gp("RNEG=101; RPOS=102; INEG=103; IPOS=104;");

        # Print branch label given a branch value.
        # The label corresponds to the branch cut.
        gp("""
        printbranch(branch)=
        {
         if( branch == RNEG, print("RNEG"));
         if( branch == RPOS, print("RPOS"));
         if( branch == INEG, print("INEG"));
         if( branch == IPOS, print("IPOS"));
         return();
         error("Unknown sqrt branch");
        }
        """.replace('\n','')
        );

        # Square root determined by branch value, which should be thought of as a
        # continuation to a half-plane.
        gp(
        """
        rac2(z, branch)=
        {
         if( branch == RNEG, return( sqrt(z) ) );
         if( branch == RPOS, return( I * sqrt(-z) ) );
         if( branch == INEG, return( sqrt(I) * sqrt(-I*z) ) );
         if( branch == IPOS, return( sqrt(-I) * sqrt(I * z) ) );
         error("Unknown sqrt branch");
        }
        """.replace('\n','')
        );

        # Best branch is the one that stays farthest away from the cut corresponding to
        # it.
        gp("""
        best_branch(z)=
        {
         my( x = real(z), y = imag(z) );
         if( abs(y) > abs(x),
          if( y > 0, return(INEG), return(IPOS))
         ,
          if( x > 0, return(RNEG), return(RPOS))
         );
        }
        """.replace('\n',''));

        # Evaluation:
        gp("EvPol(f, x) = subst(f, variable(f), x);");

        # Heuristic zero test:
        gp("IsApprox0(x) = abs(x) < 10^(-default(realprecision) / 2);");

        # Return real positive roots of a complex polynomial
        gp(
        """
        {
        polposroots( f ) =
         if( poldegree( f ) <= 0, return( [] ) );
         my( Z = polroots( f ) );
         RZ = List();
         for(i = 1,#Z,
          if( IsApprox0( imag( Z[i] ) ) && real( Z[i] ) > 10^(-default(realprecision) / 2), listput( RZ, real(Z[i]) ))
         );
         Vec(RZ);
        }
        """.replace('\n',''));

        # Return distinct roots in [0, 1] of a complex polynomial
        gp("""{
        pol01roots( f )=
         if(poldegree( f )<=0,return([]));
         my( Z = polroots(f), i, z, eps = 10^(-default(realprecision) / 2));
         RZ = List();
         for(i = 1,#Z,
          if( IsApprox0( imag(Z[i]) ),
           z = real( Z[i] );
           if( z > eps && z < 1-eps && (#RZ==0 || abs( z - RZ[#RZ]) > eps),
            listput(RZ, z)
           )
          )
         );
         Vec(RZ);
        }""".replace('\n',''));

        # Normalises a polynomial by removing the null leading terms (should fix the bug reported by Jeroen)
        gp("""{
        normalise_lc_pol(f)=
         my(d=poldegree(f),d2=-1,var=variable(f),eps=10^(-default(realprecision)/2));
         for(i=0,d,
          if(abs(polcoeff(f,i))>eps,d2=i)
         );
         sum(i=0,d2,polcoeff(f,i)*var^i);
        }""".replace('\n',''));

        # Returns values in [0, 1] where f crosses one of the branch cuts
        gp("""{
        polcrossaxes01(f)=
         my(u=real(f),v=imag(f),Z,Zf,eps=10^(-default(realprecision)/2));
         Z=pol01roots(f);
         u=normalise_lc_pol(u);
         v=normalise_lc_pol(v);
         if(u==0 || v==0, return(Z));
         Zf=prod(i=1,#Z,'t-Z[i]);
         Z=concat(Z,pol01roots(u\Zf)); Z=concat(Z,pol01roots(v\Zf));
         Z=vecsort(Z); /* TODO fix can be multiple roots there */
        }""".replace('\n',''));

        # Integrates division by sqrt(f) from a to b (complexes), assuming no branch
        # cut and f square free. Detects zeros at ends and simplifies them out by changing
        # variable.
        gp("""{
        AJ_fin_nocut(f, a, b, branch, mesh)=

         /*print("AJfin");*/
         my(t, g, h, S, z0, z1);

         g = EvPol(f, a + (b - a) * 't);

         z0 = IsApprox0( polcoeff(g, 0) );
         z1 = IsApprox0( EvPol(g, 1) );

         if(z0,

          if(z1,
           /*print("we have zeros at both ends");*/
           h = g\('t-'t^2);
           S = 6 * intnum( t = 0, 1, [1, a + (b - a) * t^2 * (3 - 2*t)] / rac2((3 - 2*t)*(2*t + 1) * EvPol(h, t^2 * (3 - 2*t)) ,branch),mesh)
          ,
           /*print("a is a zero, and b is not");*/
           h = g\ 't;
           S = 2 * intnum( t = 0, 1, [1 , a + (b - a) * t^2 ] / rac2( EvPol( h, t^2 ), branch), mesh)
          )
         ,
          if(z1,
           /*print("b is zero, and a is not");*/
           h = g\(1-'t);
           S = 2 * intnum(t = 0, 1, [1, b + (a - b) * t^2 ] / rac2( EvPol( h, 1 - t^2 ), branch), mesh)
          ,
           /*print("neither a or b are zeros");*/
           S = intnum(t = 0, 1, [1, a + (b - a) * t ] / rac2(EvPol(g, t), branch), mesh)
          )
         );
         (b - a) * S;
        }""".replace('\n',''));





        # Abelian integral from x=a straight to x=b, starting at y=ya. Also returns yb.
        # Technique: Keep track of a branch switch throughout.
        # (The correct signed root is chosen by explicit continuation.)
        gp(
        """
        AJ_fin(f, a, b, ya, mesh) =
        {
         my(g, stops, S = [0, 0], s, yprev, z, i, w, branch, ynew);

         /*print(["Int",a,b]);*/
         g = EvPol(f, a + (b - a) * 't);

         /* Checks for which 0<t<1 the image by f of the integration path crosses the axes */
         stops = polcrossaxes01( g );

         /*print(Str(#stops, "obstacles"));*/
         /* Turn these t into z */
         stops = apply( t-> a + t * (b-a), stops);
         /* Add endpoints */
         stops = concat([a], concat(stops,[b]));
         /* print(stops); */

         /* Keep track of the sheet */
         s = 1;
         yprev = ya;
         z = a;

         for(i = 1,#stops - 1,

          w = EvPol(f,z);
          /* w is in the area we're starting from now. */
          if( IsApprox0(w), w = EvPol(deriv(f),z));
          branch = best_branch(w);

          /* Choose a sqrt accordingly */
          /* print(Str("New value ",w)); */
          /* printbranch(branch); */
          if( !IsApprox0( EvPol(f, z) ),
          /* adjust sign so that we stay in the same sheet */
           ynew = rac2( EvPol(f, z), branch);
           if( s * real(ynew/yprev) < 0, s = -s );
          );
          S += s * AJ_fin_nocut(f, stops[i], stops[i+1], branch, mesh);

          /* integrate 1 step */

          /* compute new start and new y */
          z = stops[i + 1];
          yprev = s * rac2(EvPol(f, z) , branch);
         );
         [S, yprev];
        }
        """.replace('\n',''));

        # Image of P+Q-OO by Abel-Jacobi
        gp(
        """
        {
        AJ2(f, P, Q, mesh) =
         my();
         Z = polroots(f);
         /* Integrate P+Q-2(z,0) instead, it is the same up to a period */
         z = Z[1];
         AJ_fin(f, P[1], z, P[2], mesh)[1] + AJ_fin(f, Q[1], z, Q[2], mesh)[1];
        }
        """.replace('\n',''));

        # Image of P-OO by Abel-Jacobi
        gp(
        """
        {
        AJ1(f,P,mesh)=
         my();
         Z = polroots(f);
         z = Z[1]; /* Integrate P-(z,0) instead, it is the same up to a period */
         AJ_fin(f,P[1],z,P[2],mesh)[1];
        }
        """.replace('\n',''));

        gp(
        """
        {
        whatZ(f)=
         Z = polroots(f);
         z = Z[1];
         z;
        }
        """.replace('\n',''));

        # Actual initialization:
        gp.set_default('realprecision', self.prec_approx)
        gp("my_mesh=intnuminit(0,1,{})".format(str(self.mesh)));

    def __repr__(self):
        return "The class of integration functionality for the hyperelliptic curve with equation y^2 = {}".format(str(self.f))

    def change_precision(self, prec_approx):
        self.prec_approx = prec_approx
        self.CC_approx = ComplexField(self.prec_approx * log(10) / log(2))
        self.alphas_approx = [ self.CC_approx(alpha) for alpha in alphas ]
        gp.set_default('realprecision', self.prec_approx)

    def change_mesh(self, mesh):
        self.mesh = mesh
        gp("my_mesh=intnuminit(0,1,{})".format(str(self.mesh)));

    def change_roots(self, alphas):
        self.alphas = alphas
        # NOTE: Next line is very unlikely to matter.
        self.alphas_approx = [ self.CC_approx(alpha) for alpha in alphas ]

    def abel_jacobi(self, a, b, ya):
        # Wrapper of AJ_fin: calculates the integral corresponding to the path
        # between the points a and b for the give value ya of y at a.
        return gp("AJ_fin({},{},{},{},my_mesh)".format(str(self.f).replace('x', 't'), str(a), str(b), str(ya)))

    def circle_intersections(self, pair, center, radius):
        # Determines the intersection of the line determined by the pair of
        # points pair with the circle with center center and radius radius.
        z1, z2 = pair
        u = (z2 - z1)/abs(z2 - z1)
        im = (-(center - z1)/u).imag_part()
        re = sqrt(radius^2 - im^2)
        return [ center + u*self.CC_approx([re, im]), center + u*self.CC_approx([-re, im]) ]

    def base_point(self):
        # Finds a base point from which the line segments to the individual
        # roots do not come too close to the other roots.
        center = sum(self.alphas_approx) / len(self.alphas_approx)
        radius = max([ abs(alpha - center) for alpha in self.alphas_approx ])
        n = len(self.alphas_approx)
        pairs = [ [ self.alphas_approx[i], self.alphas_approx[j] ] for i in range(n) for j in range(n) if i < j ]
        intersections = reduce(lambda x, y : x + y, [ self.circle_intersections(pair, center, 2*radius) for pair in pairs ])
        logs = sorted([ (intersection - center).log() for intersection in intersections ], key = lambda z : z.imag_part())
        # Take the largest gap in the circle segments determined by the
        # intersection of lines determined by roots:
        logs.append(logs[0] + 2*self.CC_approx.pi()*self.CC_approx([0,1]))
        # TODO: Remove this line and check that things still work
        #logs = [ log for log in logs if log.imag_part() >= -0.01 and log.imag_part() <= 2*self.CC_approx.pi() - 0.01 ]
        diffs = [ abs(logs[i + 1] - logs[i]) for i in range(0, len(logs) - 1) ]
        index = diffs.index(max(diffs))
        log = (logs[index + 1] + logs[index])/2
        x0 = center + exp(log)
        return x0

    def reorder_roots(self, x0):
        # This is an order on the roots that corresponds to sweeping out from
        # the base point counterclockwise.
        # The algorithm depends on log being defined with a branch cut on the
        # negative imaginary axis!
        center = sum(self.alphas_approx) / len(self.alphas_approx)
        alphas_with_args = [ [ self.alphas[index], self.alphas_approx[index], ((self.alphas_approx[index] - x0) / (center - x0)).argument() ] for index in range(len(self.alphas)) ]
        alphas_with_args = sorted(alphas_with_args, key = lambda x : x[2])
        self.alphas = [ alphas_with_arg[0] for alphas_with_arg in alphas_with_args ]
        self.alphas_approx = [ alphas_with_arg[1] for alphas_with_arg in alphas_with_args ]

    def mumford_tuples(self):
        # Returns the combinations of the loops around the individual branch
        # points that are used in the integrals considered by Mumford.
        tups1 = [ [ 2*n, 2*n + 1 ] for n in range(self.g) ]
        tups2 = [ [ (2*n + 1)..(2*self.g) ] for n in range(self.g) ]
        return tups1 + tups2

    def paths(self, i, P0):
        # The path around the ith branch point for a loop starting at P0 with
        # radius radius.
        I = self.CC_approx([0, 1])
        x0, y0 = P0
        radius = min([ abs(self.alphas_approx[i] - self.alphas_approx[j]) for j in range(len(self.alphas_approx)) if j != i ]) / 2
        u = x0 - self.alphas_approx[i]
        u /= abs(u)
        path_to = [ x0, self.alphas_approx[i] + radius * u ]
        paths_around = [ [ self.alphas_approx[i] + radius * u * I^k, self.alphas_approx[i] + radius * u * I^(k + 1) ] for k in range(4) ]
        path_from = [ self.alphas_approx[i] + radius * u, x0 ]
        paths = [ path_to ] + paths_around + [ path_from ]
        return paths

    def cycle_integrals(self, i, P0):
        # Integrals of the canonical base of differentials over a loop around
        # the ith branch point when starting with a point P0 on the curve
        # (which in particular includes some information on the y-coordinate).
        paths = self.paths(i, P0)
        cycle_ints = vector([ self.CC_approx(0) for n in range(self.g) ])
        y0 = P0[1]
        for path in paths:
            path_int, y0 = self.abel_jacobi(path[0], path[1], y0)
            cycle_ints += vector([ self.CC_approx(comp) for comp in path_int ])
        return cycle_ints

    def calculate_periods(self):
        # Calculates the periods for a Mumford basis after reordering the
        # roots.
        # We could approximate alpha in what follows, but it does not seem
        # worth the trouble.
        self.index = tuple([ self.prec_approx, self.mesh ])
        if not self.index in self.period_dict:
            x0 = self.base_point()
            # NOTE: Side effect of reordering the roots occurs here and is
            # crucial.
            self.reorder_roots(x0)
            y0 = sqrt(f(x0))
            P0 = [ x0, y0 ]
            P0c = [ x0, -y0 ]
            cycle_intss_pairs = [ [ self.cycle_integrals(i, P0), self.cycle_integrals(i, P0c) ] for i in range(2*self.g + 1) ]
            tups = self.mumford_tuples()
            Omega = Matrix([ sum([ cycle_intss_pairs[tup[i]][i % 2] for i in range(len(tup)) ]) for tup in tups ])
            Omega1 = Omega.submatrix(row = 0, col = 0, nrows = self.g, ncols = self.g)
            Omega2 = Omega.submatrix(row = self.g, col = 0, nrows = self.g, ncols = self.g)
            self.period_dict[self.index] = [ Omega, Omega2 * Omega1^(-1) ]

    def period_matrix(self):
        self.calculate_periods()
        return self.period_dict[self.index][0]

    def tau(self):
        self.calculate_periods()
        return self.period_dict[self.index][1]
