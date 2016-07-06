"""
 *  File for debugging
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
"""
print PeriodInfo.tau()
"""

"""
Chars = PeriodInfo.chars
print Chars.fundamental_etas()
Chars.fill_set_from_eta_dictionary_fundamental()
print Chars.set_from_eta_dictionary()
print all([ tuple(Chars.eta_from_set(Chars.set_from_eta_dict[eta])) == eta for eta in Chars.set_from_eta_dict.keys() ])
"""

"""
SymplecticAction = PeriodInfo.symplectic_action
print SymplecticAction.special_matrix_dictionary()
print SymplecticAction.special_action_information()
"""

"""
IntegrationApproximation = PeriodInfo.integration_approx
print IntegrationApproximation.period_matrix()
print IntegrationApproximation.tau()
print IntegrationApproximation.alphas
print IntegrationApproximation.alphas_approx
"""

"""
# Use "fake" dictionary:
ThomaeValues = PeriodInfo.thomae_values
ThomaeValues.chars.fill_set_from_eta_dictionary_fundamental()
etas_fund = ThomaeValues.chars.fundamental_etas()
print etas_fund
ThomaeValues.fill_thomae_from_eta_dictionary(etas_fund)
av_dict = { eta : ThomaeValues.thomae_from_eta_dict[eta] for eta in etas_fund }
print av_dict
iterator = BorchardtIterator(av_dict)
iterator.iterate()
print iterator.av_dict
print iterator.iterates()
"""

"""
Thetas = PeriodInfo.thetas
etas_fund = Thetas.symplectic_action.chars.fundamental_etas()
tau_approx = Thetas.integration_approx.tau()
print tau_approx
print Thetas.theta_approximation_dictionary(etas_fund, tau_approx)
print Thetas.special_etas()
print Thetas.theta0_borchardt_start(etas_fund)
print Thetas.theta_square_dictionary(etas_fund)
special_info = Thetas.symplectic_action.special_action_information()
action_dict = special_info.values().pop()
print action_dict
print Thetas.theta0_transformed_borchardt_start(action_dict)
print Thetas.theta0_transformed_square(action_dict)
"""
