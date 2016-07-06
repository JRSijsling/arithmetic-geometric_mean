"""
 *  Initialization of the package
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

import os
# The following line is a bad solution:
if not '__agmdir__' in globals():
    __agmdir__ = os.getenv("PWD") + "/"

load(__agmdir__ + "Characteristics.sage")
load(__agmdir__ + "SymplecticAction.sage")
load(__agmdir__ + "IntegrationApproximation.sage")
load(__agmdir__ + "ThomaeValues.sage")
load(__agmdir__ + "BorchardtIterator.sage")
load(__agmdir__ + "Thetas.sage")
load(__agmdir__ + "PeriodInformation.sage")
