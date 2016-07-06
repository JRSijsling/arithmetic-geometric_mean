"""
 *  Correspondence between characteristics and sets of branch points
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

class Characteristics:
    def __init__(self, g):
        self.g = g
        # Derived objects:
        self.F = FiniteField(2)
        self.V = VectorSpace(self.F, 2*self.g)
        self.eta0 = tuple(self.V(0))
        # For caching:
        self.eta_from_set_dict = { }
        self.set_from_eta_dict = { }
        self.fundamental_added = False

    def __repr__(self):
        return "The class of characteristic functionality in genus {}".format(self.g)

    def eta_from_branch_point(self, n):
        # Conversion of a single branch point to a characteristic.
        # Dupont page 191.
        # First we deal with the two final points, which Dupont and Mumford do
        # not mention:
        if n == 2*self.g + 1:
            return self.V([ self.F(0) for i in range(2*self.g) ])
        if n == 2*self.g:
            return self.V([ self.F(0) for i in range(self.g) ] + [ self.F(1) for i in range(self.g) ])
        # The other points are as prescribed in loc. cit.:
        m = n // 2
        s = n % 2
        a = [ self.F(kronecker_delta(i, m)) for i in range(self.g) ]
        b = [ self.F(0) for i in range(self.g) ]
        for i in range(m + s):
            b[i] = self.F(1)
        return self.V(a + b)

    def eta_from_set(self, S):
        # Conversion of a set (of labels of branch points) to a characteristic.
        # Caches the result.
        if not S in self.eta_from_set_dict.keys():
            self.eta_from_set_dict[S] = self.V(sum([ self.eta_from_branch_point(n) for n in S ]))
        return self.eta_from_set_dict[S]

    def set_from_eta(self, eta):
        # Finds subset of even cardinality corresponding to fundamental
        # characteristic.
        # Dupont page 191.
        # Caches the result.
        # Sets for the basis elements:
        if not eta in self.set_from_eta_dict.keys():
            ess1 = [ Set([ 0..2*i ]) + Set([ 2*self.g + 1 ]) for i in range(self.g) ]
            ess2 = [ Set([ 2*i, 2*i + 1 ]) for i in range(self.g) ]
            ess = ess1 + ess2
            S = Set([ ])
            for i in range(2*self.g):
                if self.F(eta[i]) != 0:
                    S = S.symmetric_difference(ess[i])
            # Avoid final branch point (which can correspond to infinity) by
            # taking a complement if necessary:
            if 2*self.g + 1 in S:
                S = Set([0..2*self.g + 1]).difference(S)
            self.set_from_eta_dict[eta] = S
        return self.set_from_eta_dict[eta]

    def create_second_component_dictionary(self, g):
        # Creates a dictionary with subsets of even cardinality that correspond to
        # the fundamental characteristics (those starting with 0).
        if g == 1:
            V = VectorSpace(self.F, 1)
            return { tuple(V([0])): Set([ ]), tuple(V([1])) : Set([ 0, 1 ]) }
        # Recursive creation:
        es = Set([ 2*g - 2, 2*g - 1 ])
        dictionary_old = self.create_second_component_dictionary(g = g - 1)
        dictionary = { }
        for tup in dictionary_old.keys():
            S = dictionary_old[tup]
            dictionary[tup + tuple([self.F(0)])] = S
            dictionary[tup + tuple([self.F(1)])] = S.symmetric_difference(es)
        return dictionary

    def fill_set_from_eta_dictionary_fundamental(self):
        # Fills the dictionary with the values for the fundamental
        # characteristics.
        if not self.fundamental_added:
            second_part_dict = self.create_second_component_dictionary(self.g)
            dictionary_to_add = { tuple(self.V([ 0 for i in range(self.g)] + list(second_part))) : second_part_dict[second_part] for second_part in second_part_dict.keys() }
            self.set_from_eta_dict.update(dictionary_to_add)
            # NOTE: Side effect introduced here and used later.
            self.etas_fund = dictionary_to_add.keys()
            self.fundamental_added = True

    def fill_set_from_eta_dictionary_all(self):
        # Same as above, but for all theta characteristics.
        # Main sets used:
        ess1 = [ Set([ 2*i, 2*i + 1 ]) for i in range(self.g) ][::-1]
        ess2 = [ Set([ 0..2*i ]) + Set([ 2*self.g + 1 ]) for i in range(self.g) ][::-1]
        ess = ess1 + ess2
        # "Reverse dictionary":
        self.set_from_eta_dict = { tuple(self.V(0)) : Set([ ]) }
        for es in ess:
            eta = self.eta_from_set(es)
            for key in self.set_from_eta_dict.keys():
                # Again avoid branch point a infinity:
                S = self.set_from_eta_dict[key].symmetric_difference(es)
                if 2*self.g + 1 in S:
                    S = Set([0..2*self.g + 1]).difference(S)
                self.set_from_eta_dict[tuple(self.V(key) + eta)] = S

    def eta_from_set_dictionary(self):
        return self.eta_from_set_dict

    def set_from_eta_dictionary(self):
        return self.set_from_eta_dict

    def fundamental_etas(self):
        if not hasattr(self, "etas_fund"):
            self.etas_fund = [ ]
            for n in range(2^self.g):
                m = n
                eta = [ ]
                for index in range(2*self.g):
                    eta.insert(0, self.F(m))
                    m //= 2
                self.etas_fund.append(tuple(eta))
        return self.etas_fund
