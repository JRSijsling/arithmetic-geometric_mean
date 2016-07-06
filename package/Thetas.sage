"""
 *  Calculation of dictionaries of theta-values to high precision
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

class Thetas:
    def __init__(self, integration_approx, symplectic_action, thomae_values, term_bound = 2, agm_stopping_precision = 2^(-20)):
        self.integration_approx = integration_approx
        self.symplectic_action = symplectic_action
        self.thomae_values = thomae_values
        self.term_bound = term_bound
        self.agm_stopping_precision = agm_stopping_precision
        # Derived objects:
        self.eta0 = self.symplectic_action.chars.eta0
        self.g = (len(self.integration_approx.alphas) - 1) // 2
        self.CC_work = self.integration_approx.alphas[0].parent()
        self.CC_approx = self.integration_approx.CC_approx
        # For caching:
        # TODO: This class could have a lot more of these, and (see below) we
        # should recycle old values when the precisions change. Not now.
        self.theta_square_dict = { }
        self.theta0_transf_square_dict = { }

    def __repr__(self):
        return "The class of theta functionality for the roots {}".format(str(self.integration_approx.alphas))

    def change_term_bound(self, term_bound):
        # Modifies the number of terms used.
        self.term_bound = term_bound
        # TODO: We could cache the old results. Not for now, since it is done
        # one level up; instead we clear the dictionaries.
        self.theta_square_dict = { }
        self.theta0_transf_square_dict = { }

    def change_agm_stopping_precision(self, agm_stopping_precision):
        # Determines from which point on only good roots of the AGM are used.
        # TODO: Find a mathematical way to deal with this.
        self.agm_stopping_precision = agm_stopping_precision
        # TODO: We could cache the old results. Not for now, since it is done
        # one level up; instead we clear the dictionaries.
        self.theta_square_dict = { }
        self.theta0_transf_square_dict = { }

    def theta_approximation_from_eta(self, eta, tau_approx):
        # An approximation of the values of the theta function with
        # characteristic eta. Note that the internal bound term_bound that
        # controls the number of terms is used.
        # TODO: Cache values, especially for the main tau_approx.
        # TODO: Could make creation of ZZg global.
        ZZg = cartesian_product_iterator([ [-self.term_bound..self.term_bound] for i in range(self.g) ])
        a = eta[0:self.g]
        b = eta[self.g:2*self.g]
        am = Matrix(self.CC_approx, [ [ Integers()(x) for x in a ] ])
        bm = Matrix(self.CC_approx, [ [ Integers()(x) for x in b ] ])
        theta_approx = self.CC_approx(0)
        for n in ZZg:
            nm = Matrix(self.CC_approx, [ [ x for x in n ] ])
            nam = nm + (1/2)*am
            e = (nam * tau_approx * nam.transpose() + nam * bm.transpose())[0, 0]
            theta_approx += exp(e * self.CC_approx.pi() * self.CC_approx([0,1]))
        return theta_approx

    def theta_approximation_dictionary(self, etas, tau_approx):
        # TODO: Cache values, especially for the main tau_approx.
        return { eta : self.theta_approximation_from_eta(eta, tau_approx) for eta in etas }

    def calculate_special_etas(self):
        # Minimal set of etas needed to calculate the entries of the period matrix.
        # TODO: We could also construct the even smaller set of transformations
        # of the zero element and use those later. Again, not for now.
        if not hasattr(self.symplectic_action.chars, "etas_special"):
            special_action_info = self.symplectic_action.special_action_information()
            self.symplectic_action.chars.etas_special = reduce(lambda x, y : x.union(y), [ Set([ tuple(eta_epsilon[0]) for eta_epsilon in action_dict[1].values() ]) for action_dict in special_action_info.values() ])

    def special_etas(self):
        self.calculate_special_etas()
        return self.symplectic_action.chars.etas_special

    def theta0_borchardt_start(self, etas):
        # Start of iteration for theta_{0,0} (tau) via the arithmetic-geometric
        # mean. It is slightly artificial to distinguish this from the other
        # cases, but it seems better for readability.
        # Dupont page 201.
        # NOTE: Can add more etas for more quotients. We usually add the
        # special ones only. At any rate, the function below is for the
        # transformation is note quite a generalization of this one (as there
        # the etas are not specified).
        # TODO: Add caching (needed later, but indexing the results will be
        # more difficult). The value of theta0_square can then be calculated
        # separately.
        self.thomae_values.fill_thomae_from_eta_dictionary(etas)
        theta_approx_dict = self.theta_approximation_dictionary(etas, self.integration_approx.tau())
        theta0_approx_square = theta_approx_dict[self.eta0]^2
        # Quotients of the theta approximations will later determine the
        # correct roots for the roots of the Thomae quotients (which will have
        # higher precision):
        theta_square_quot_approx_dict = { eta : theta_approx_dict[eta]^2 / theta0_approx_square for eta in etas }
        # Now we create the roots of the Thomae quotients. First we create the
        # Thomae values themselves:
        theta_square_quot_dict = { eta : sqrt(self.thomae_values.thomae_from_eta_dict[eta]) for eta in etas }
        # Then we quotient out by the value at the zero characteristic, just as
        # we did for the approximations of the theta values:
        theta_square_quot_dict = { eta : theta_square_quot_dict[eta] / theta_square_quot_dict[self.eta0] for eta in etas }
        # We choose minus signs by using the approximation:
        for eta in etas:
            if abs(theta_square_quot_dict[eta] - theta_square_quot_approx_dict[eta]) > abs(theta_square_quot_dict[eta] + theta_square_quot_approx_dict[eta]):
                theta_square_quot_dict[eta] *= -1
        return [ theta_square_quot_dict, theta0_approx_square ]

    def fill_theta_square_dictionary(self, etas):
        # Use iterator to get theta0 to high precision.
        # NOTE: At least all fundamental entries should be in the input
        # dictionary etas for now.
        # TODO: Add more caching (right now it still has to recalculate a lot
        # even if one entry fails).
        if not all([ self.theta_square_dict.has_key(eta) for eta in etas ]):
            theta_square_quot_dict, theta0_square_approx = self.theta0_borchardt_start(etas)
            etas_fund = self.symplectic_action.chars.fundamental_etas()
            borchardt_iterator = BorchardtIterator({ eta : theta_square_quot_dict[eta] for eta in etas_fund })
            tau_approx_fold = self.integration_approx.tau()
            only_good = False
            while True:
                if not only_good:
                    # First use approximations:
                    tau_approx_fold *= 2
                    theta_approx_dict = self.theta_approximation_dictionary(etas_fund, tau_approx_fold)
                    iterate_approx = { eta : theta_approx_dict[eta]^2 / theta0_square_approx for eta in etas_fund }
                    borchardt_iterator.iterate(only_good = False, iterate_approx = iterate_approx)
                else:
                    # Then switch to good iterates:
                    borchardt_iterator.iterate(only_good = True)
                avs = borchardt_iterator.av_dict.values()
                maximum_difference = max([ abs(av - avs[0]) for av in avs ])
                # TODO: Justify the 3.
                if abs(maximum_difference) < self.CC_work(2^(-self.CC_work.precision() + 3)):
                    break
                only_good = maximum_difference < self.CC_work(self.agm_stopping_precision)
            theta0_square = 1 / avs[0]
            theta_square_dict_new = { eta : theta0_square * theta_square_quot_dict[eta] for eta in etas }
            self.theta_square_dict.update(theta_square_dict_new)

    def theta0_square(self):
        self.fill_theta_square_dictionary(self.symplectic_action.chars.fundamental_etas())
        return self.theta_square_dict[self.eta0]

    def theta_square(self, eta):
        self.fill_theta_square_dictionary(Set(self.symplectic_action.chars.fundamental_etas()).union(Set([eta])))
        return self.theta_square_dict[eta]

    def theta_square_dictionary(self, etas):
        self.fill_theta_square_dictionary(Set(self.symplectic_action.chars.fundamental_etas()).union(Set(etas)))
        return { eta : self.theta_square_dict[eta] for eta in etas }

    def theta0_transformed_borchardt_start(self, action_dict):
        # Start of iteration for theta_{0,0} (gamma tau) via the
        # arithmetic-geometric mean. The input is a list consisting of a matrix
        # and its action on characteristics along with the corresponding roots
        # of unity.
        # Dupont page 201.
        gamma, eta_epsilon_dict = action_dict
        etas_fund = self.symplectic_action.chars.fundamental_etas()
        etas_transf = [ eta_epsilon_dict[eta][0] for eta in etas_fund ]
        self.thomae_values.fill_thomae_from_eta_dictionary(etas_transf)
        tau_transf_approx = self.symplectic_action.action_on_tau(gamma, self.integration_approx.tau())
        theta_transf_approx_dict = self.theta_approximation_dictionary(etas_fund, tau_transf_approx)
        theta0_transf_approx_square = theta_transf_approx_dict[self.eta0]^2
        # Quotients of the theta approximations will later determine the
        # correct roots for the roots of the Thomae quotients (which will have
        # higher precision):
        theta_transf_square_quot_approx_dict = { eta : theta_transf_approx_dict[eta]^2 / theta0_transf_approx_square for eta in etas_fund }
        # Now we create the roots of the Thomae quotients. First we create the
        # Thomae values themselves after finding their quasi-transformations:
        theta_transf_square_quot_dict = { eta : self.CC_work([0,1])^(eta_epsilon_dict[eta][1]) * self.thomae_values.thomae_from_eta_dict[tuple(eta_epsilon_dict[eta][0])] for eta in etas_fund }
        theta_transf_square_quot_dict = { eta : sqrt(theta_transf_square_quot_dict[eta]) for eta in etas_fund }
        # Then we quotient out by the value at the zero characteristic, just as
        # we did for the approximations of the theta values:
        theta_transf_square_quot_dict = { eta : theta_transf_square_quot_dict[eta] / theta_transf_square_quot_dict[self.eta0] for eta in etas_fund }
        # We choose minus signs by using the approximation:
        for eta in etas_fund:
            if abs(theta_transf_square_quot_dict[eta] - theta_transf_square_quot_approx_dict[eta]) > abs(theta_transf_square_quot_dict[eta] + theta_transf_square_quot_approx_dict[eta]):
                theta_transf_square_quot_dict[eta] *= -1
        return [ theta_transf_square_quot_dict, theta0_transf_approx_square ]

    def theta0_transformed_square(self, action_dict):
        # Use iterator to get the transformed theta0-value to high precision.
        gamma = action_dict[0]
        theta_transf_square_quot_dict, theta0_transf_square_approx = self.theta0_transformed_borchardt_start(action_dict)
        etas_fund = self.symplectic_action.chars.fundamental_etas()
        borchardt_iterator = BorchardtIterator({ eta : theta_transf_square_quot_dict[eta] for eta in etas_fund })
        tau_transf_approx_fold = self.symplectic_action.action_on_tau(gamma, self.integration_approx.tau())
        only_good = False
        while True:
            if not only_good:
                # First use approximations:
                tau_transf_approx_fold *= 2
                theta_transf_approx_dict = self.theta_approximation_dictionary(etas_fund, tau_transf_approx_fold)
                iterate_approx = { eta : theta_transf_approx_dict[eta]^2 / theta0_transf_square_approx for eta in theta_transf_approx_dict.keys() }
                borchardt_iterator.iterate(only_good = False, iterate_approx = iterate_approx)
            else:
                # Then switch to good iterates:
                borchardt_iterator.iterate(only_good = True)
            avs = borchardt_iterator.av_dict.values()
            maximum_difference = max([ abs(av - avs[0]) for av in avs ])
            if abs(maximum_difference) < self.CC_work(2^(-self.CC_work.precision())):
                break
            only_good = maximum_difference < self.CC_work(self.agm_stopping_precision)
        return 1 / avs[0]
