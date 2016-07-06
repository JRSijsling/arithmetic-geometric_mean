"""
 *  Example calculations of period matrices to high precision
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

load("Initialize.sage")

CC.<I> = ComplexField()
alphas = [ 1, 2, 3, 4, 5, 7 ]
alphas = [ 1, 2 + I, 2 - I, 5 + I, 5 - I, 7 ]
alphas = [ -1, 0, 1, 2, 3, 7 ]
alphas = [ -1, 1, 2, 3, 7 ]
alphas = [ CC(alpha) for alpha in alphas ]
R.<x> = PolynomialRing(CC)
f = prod([ x - alpha for alpha in alphas ])

PeriodInfo = PeriodInformation(f)
print PeriodInfo.tau()
