"""
 *  Class for iterating Thomae values after Borchardt
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

class BorchardtIterator:
    def __init__(self, av_dict):
        # Save original input:
        self.av_dict_orig = av_dict
        # Dictionary for the actual integration:
        self.av_dict = av_dict
        # Recover vector space over which we work:
        self.V = vector(list(self.av_dict.keys()[0])).parent()
        self.eta0 = tuple(self.V(0))

    def __repr__(self):
        return "The class of Borchardt functionality for the dictionary {}".format(str(self.av_dict_orig))

    def iterates(self, only_good = False):
        # Input av_dict is a dictionary with as keys the elements of the base
        # vector space self.V.
        # Dupont page 155.
        etas = self.av_dict.keys()
        if only_good:
            # In this case the first root determines the others.
            # TODO: Treat boundary cases (not quite specified by Dupont).
            root0 = sqrt(self.av_dict[self.eta0])
            alphav_dict = { self.eta0 : sqrt(self.av_dict[self.eta0]) }
            for eta in etas:
                root = sqrt(self.av_dict[eta])
                if abs(root - root0) < abs(root + root0):
                    alphav_dict[eta] = root
                else:
                    alphav_dict[eta] = -root
            alphav_dicts = [ alphav_dict ]
        else:
            alphav_dict0 = { eta : sqrt(self.av_dict[eta]) for eta in etas }
            # Possible sign choices (the range runs from 1 because changing it everywhere has no effect)
            cp = cartesian_product_iterator([ [1, -1] for i in range(1, len(alphav_dict0)) ])
            alphav_dicts = [ ]
            for signs_end in cp:
                signs = [ 1 ] + list(signs_end)
                etas_signs_zip = zip(etas, signs)
                alphav_dict = { etas_signs[0] : etas_signs[1] * alphav_dict0[etas_signs[0]] for etas_signs in etas_signs_zip }
                alphav_dicts.append(alphav_dict)
        av_dicts = [ ]
        for alphav_dict in alphav_dicts:
            av_dict = { eta : 0 for eta in etas }
            for eta1 in etas:
                for eta2 in etas:
                    index = tuple(self.V(eta1) + self.V(eta2))
                    av_dict[index] += alphav_dict[eta1] * alphav_dict[eta2]
            av_dict = { eta : av_dict[eta] / len(etas) for eta in etas }
            av_dicts.append(av_dict)
        return av_dicts

    def iterate(self, only_good = True, iterate_approx = "not set"):
        # Choose the correct iterate. In case only_good is set, we get the
        # "good" iterate in the sense of Dupont page 155. Otherwise a
        # comparison with an approximation is made and the closest result is
        # returned.
        if not only_good and iterate_approx == "not set":
            raise ValueError("Please specify a value for iterate_approx to obtain the closest iterate, or set only_good to True to obtain the good iterate")
        if only_good:
            self.av_dict = self.iterates(only_good = True)[0]
        else:
            iterates = self.iterates(only_good = False)
            max_diff_abss = [ max([ abs(iterate[eta] - iterate_approx[eta]) for eta in iterate.keys() ]) for iterate in iterates ]
            self.av_dict = iterates[max_diff_abss.index(min(max_diff_abss))]

    def current_values(self):
        return self.av_dict
