"""
 *  Creation of special symplectic matrices for period calculation
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

class SymplecticAction:
    def __init__(self, chars):
        self.chars = chars
        self.g = self.chars.g
        # Free ZZ-module needed in what follows:
        self.M = FreeModule(Integers(), 2*self.g)
        # For caching:
        # TODO: That currently fails
        self.action_on_etas_fund_dict = { }

    def __repr__(self):
        return "The class of symplectic actions in genus {}".format(str(self.g))

    def calculate_special_matrix_dictionary(self):
        # Special matrices in symplectic group.
        # Dupont page 201.
        if not hasattr(self, "gamma_dict"):
            S = MatrixSpace(Integers(), self.g, self.g)
            J = block_matrix([ [ S(0), S(-1) ], [ S(1), S(0) ] ])
            A = S(1)
            B = S(0)
            C = S(0)
            D = S(1)
            self.gamma_dict = { }
            for i in range(self.g):
                for j in range(i + 1):
                    B = S(0)
                    B[i,j] = B[j,i] = 1
                    M = block_matrix([ [A, B], [C, D] ], subdivide = False)
                    self.gamma_dict[i,j] = (J * M)^2

    def special_matrix_dictionary(self):
        if not hasattr(self, "gamma_dict"):
            self.calculate_special_matrix_dictionary()
        return self.gamma_dict

    def split_matrix(self, gamma):
        A = gamma.submatrix(0, 0, self.g, self.g)
        B = gamma.submatrix(0, self.g, self.g, self.g)
        C = gamma.submatrix(self.g, 0, self.g, self.g)
        D = gamma.submatrix(self.g, self.g, self.g, self.g)
        return A, B, C, D

    def action_on_tau(self, gamma, tau):
        # The fractional linear transformation of tau by gamma.
        A, B, C, D = self.split_matrix(gamma)
        return (A * tau + B) * (C * tau + D)^(-1)

    def action_on_etas_fundamental(self, gamma):
        # Action on the characteristics used in transformation behavior of
        # theta functions.
        # Dupont page 128.
        # TODO: Upcoming trick fails because for some a reason a tuple is still
        # not hashable.
        #if not index in self.action_on_etas_fund_dict.keys():
        A, B, C, D = self.split_matrix(gamma)
        m00 = (C * D.transpose()).diagonal()
        m01 = (A * B.transpose()).diagonal()
        m0 = self.M(m00 + m01)
        m0_transf = gamma.transpose() * m0
        eta_epsilon_dict = { }
        for eta in self.chars.fundamental_etas():
            alpha_beta = self.M(gamma.transpose() * self.M(eta) - m0_transf)
            eta_transf = tuple(self.chars.V(alpha_beta))
            alpha = Matrix(Integers(), self.g, 1, alpha_beta[0:self.g].list())
            beta = Matrix(Integers(), self.g, 1, alpha_beta[self.g:2*self.g].list())
            diag = Matrix(Integers(), self.g, 1, (A*B.transpose()).diagonal())
            rt_exp = 2 * (B*alpha).transpose() * C*beta\
                    - (D*alpha).transpose() * B*alpha\
                    - (C*beta).transpose() * A*beta\
                    + 2 * diag.transpose() * (D*alpha - C*beta)
            eta_epsilon_dict[eta] = [ eta_transf, rt_exp[0, 0] % 4 ]
        return eta_epsilon_dict
        #    self.action_on_etas_fund_dict[index] = eta_epsilon_dict
        #return self.action_on_etas_fund_dict[index]

    def actions_on_etas_fundamental(self, gamma_dict):
        # Dictionary of the above for a fixed dictionary of transformations
        # gamma_dict. Also returns all the transformed characteristics.
        return { indices : [ gamma_dict[indices], self.action_on_etas_fundamental(gamma_dict[indices]) ] for indices in gamma_dict.keys() }

    def calculate_special_action_information(self):
        if not hasattr(self, "special_action_info"):
            self.special_action_info = self.actions_on_etas_fundamental(self.special_matrix_dictionary())

    def special_action_information(self):
        self.calculate_special_action_information()
        return self.special_action_info
