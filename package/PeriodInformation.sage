"""
 *  Main functionality
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

class PeriodInformation:
    def __init__(self, f, prec = 300, prec_buffer = 20, prec_approx = 20, mesh = 2, term_bound = 2, agm_stopping_precision = 2^(-20)):
        self.f = f
        self.prec = prec
        self.prec_buffer = prec_buffer
        # Derived objects:
        self.g = (self.f.degree() - 1) // 2
        CC_work.<I> = ComplexField((self.prec + self.prec_buffer) * log(10) / log(2))
        self.CC_work = CC_work
        CC_return.<I> = ComplexField(self.prec * log(10) / log(2))
        self.CC_return = CC_return
        alphas = [ tup[0] for tup in self.f.roots(self.CC_work) ]
        # Classes needed for the calculations:
        # NOTE: We only need theta around but keep these for debugging.
        self.integration_approx = IntegrationApproximation(self.f, alphas, prec_approx = prec_approx, mesh = mesh)
        self.chars = Characteristics(self.g)
        self.symplectic_action = SymplecticAction(self.chars)
        self.thomae_values = ThomaeValues(alphas, self.chars)
        self.thetas = Thetas(self.integration_approx, self.symplectic_action, self.thomae_values, term_bound = term_bound, agm_stopping_precision = agm_stopping_precision)
        # The final results are kept track of in a dictionary:
        self.tau_dict = { }

    def __repr__(self):
        return "The period matrix information of the hyperelliptic curve defined by y^2 = {}".format(str(self.f))

    def change_precision(self, prec):
        # Changes main precision, so the precision of the period matrix returned in the end
        # TODO: This is wrong because other objects inherit. Modify other updates too.
        self.prec = prec
        CC_work.<I> = ComplexField((self.prec + self.prec_buffer) * log(10) / log(2))
        self.CC_work = CC_work
        CC_return.<I> = ComplexField(self.prec * log(10) / log(2))
        self.CC_return = CC_return
        alphas = self.f.roots(self.CC_work)
        self.thetas.integration_approx.change_roots(alphas)
        self.thetas.thomae_values.change_roots(alphas)

    def change_precision_buffer(self, prec_buffer):
        # Changes main precision, so the precision of the period matrix returned in the end
        self.prec_buffer = prec_buffer
        CC_work.<I> = ComplexField((self.prec + self.prec_buffer) * log(10) / log(2))
        self.CC_work = CC_work
        alphas = self.f.roots(self.CC_work)
        self.thetas.integration_approx.change_roots(alphas)
        self.thetas.thomae_values.change_roots(alphas)

    def change_approximation_precision(self, prec_approx):
        # Changes precision of approximation used
        self.integration_approx.change_precision(prec_approx)

    def change_approximation_mesh(self, mesh):
        # Changes mesh used in approximation
        self.integration_approx.change_mesh(mesh)

    def change_term_bound(self, term_bound):
        # Changes number of terms used when evaluating theta functions
        self.thetas.change_term_bound(term_bound)

    def change_agm_stopping_precision(self, agm_stopping_precision):
        # Determines from which point on only good roots of the AGM are used.
        # TODO: Find a mathematical way to deal with this.
        self.thetas.change_agm_stopping_precision(agm_stopping_precision)

    def calculate_periods(self):
        # The main function.
        index = tuple([self.prec, self.prec_buffer, self.thetas.integration_approx.prec_approx, self.thetas.integration_approx.mesh, self.thetas.term_bound, self.thetas.agm_stopping_precision])
        if not index in self.tau_dict:
            # Calculating approximate value and passing on reordering of roots.
            tau_approx = self.thetas.integration_approx.tau()
            print tau_approx
            self.thetas.thomae_values.change_roots(self.thetas.integration_approx.alphas)
            # Create special matrices and their actions, along with squares of
            # theta values:
            self.thetas.thomae_values.chars.fill_set_from_eta_dictionary_fundamental()
            special_action_info = self.thetas.symplectic_action.special_action_information()
            etas_special = self.thetas.special_etas()
            # NOTE: In the next line we can get away with a smaller dictionary.
            # (Also mentioned in Thetas.sage.)
            theta_square_dict = self.thetas.theta_square_dictionary(etas_special)
            # Initialize final outcome as zero:
            tau = Matrix(self.CC_return, self.g, self.g, 0)
            for i in range(self.g):
                # Backwards here to make antidiagonal calculation possible:
                for j in range(i + 1)[::-1]:
                    action_dict = special_action_info[i,j]
                    theta0_transf_square = self.thetas.theta0_transformed_square(action_dict)
                    # TODO: Change structure a bit to avoid these unclear
                    # indices.
                    eta_transf = action_dict[1][self.thetas.eta0][0]
                    # TODO: Modify following line (wrong)
                    quot = theta0_transf_square / theta_square_dict[eta_transf]
                    if i == j:
                        entry = quot
                    else:
                        entry = sqrt(quot + tau[i,i]*tau[j,j])
                    # Take correct fourth root:
                    entries = [ self.CC_work([0, 1])^e * entry for e in range(4) ]
                    diffs = [ abs(entry - tau_approx[i,j]) for entry in entries ]
                    tau[i,j] = tau[j,i] = self.CC_return(entries[diffs.index(min(diffs))])
            self.tau_dict[index] = tau

    def tau(self):
        self.calculate_periods()
        return self.tau_dict[self.prec, self.prec_buffer, self.integration_approx.prec_approx, self.integration_approx.mesh, self.thetas.term_bound, self.thetas.agm_stopping_precision]
