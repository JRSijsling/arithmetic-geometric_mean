"""
 *  Calculation of Thomae values for a given set of roots
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

class ThomaeValues:
    def __init__(self, alphas, chars):
        self.alphas = alphas
        self.chars = chars
        self.g = self.chars.g
        # For caching:
        self.thomae_from_set_dict = { }
        self.thomae_from_eta_dict = { }

    def __repr__(self):
        return "The class of Thomae value functionality for the hyperelliptic curve corresponding to the roots {}".format(str(self.alphas))

    def change_roots(self, alphas):
        # NOTE: For now this clears the cache.
        self.alphas = alphas
        self.thomae_from_set_dict = { }
        self.thomae_from_eta_dict = { }

    def thomae_from_set(self, S):
        # Thomae value for a subset S subset of {0..2g-1}, where alphas is a
        # set of roots. Even cardinality assumed.
        # Dupont page 192.
        if not S in self.thomae_from_set_dict:
            U = Set([ 2*i for i in range(self.g + 1) ])
            S = Set(S)
            # To deal with a branch point at infinity:
            if 2*self.g + 1 in S:
                S = Set(range(len(self.alphas))).difference(S)
            SU = S.symmetric_difference(U)
            SUc = Set(range(len(self.alphas))).difference(SU)
            s = (-1)^(len(S.intersection(U)))
            if len(SU) != self.g + 1:
                return 0
            self.thomae_from_set_dict[S] = s * prod([ 1/(self.alphas[i] - self.alphas[j]) for i in SU for j in SUc ])
        return self.thomae_from_set_dict[S]

    def thomae_from_eta(self, eta):
        if not eta in self.thomae_from_eta_dict:
            self.thomae_from_eta_dict[eta] = self.thomae_from_set(self.chars.set_from_eta(eta))
        return self.thomae_from_eta_dict[eta]

    def fill_thomae_from_eta_dictionary(self, etas):
        # Calculates the Thomae values for the provided characteristics.
        for eta in etas:
            if not eta in self.thomae_from_eta_dict.keys():
                # NOTE: In fact evaluating the next line suffices:
                self.thomae_from_eta_dict[eta] = self.thomae_from_eta(eta)

    def thomae_from_set_dictionary(self):
        return self.thomae_from_set_dict

    def thomae_from_eta_dictionary(self):
        return self.thomae_from_eta_dict
