#!/usr/bin/env python

import numpy as N
import matplotlib as P


def net_rate(n_e,n_xi,n_xip1,n_xim1,alpha_rec, alpha_recm1, I_rat, I_ratm1):
	# Time dependence of ionization balance (dn_Xi/dt)
	# electron density, relative densities of ions (i, i+1, i-1), recombination coefficients, ionization rates

	# Where to get these from? Cloudy? Xstar? <- some sort of utility exists for taking rates etc from xstar, might be useful

	recomb    = - n_xi * n_e * alpha_recm1 # recombination of Xi -> Xi-1
	recomb_p1 = n_xip1 * n_e * alpha_rec   # recombination of Xi+1 -> Xi
	ioniz     = - n_xi * I_rat             # ionization of Xi -> Xi+1
	ioniz_m1  = n_xim1 * I_ratm1           # ionization of Xi-1 -> Xi
	return recomb + recomb_p1 + ioniz + ioniz_m1

def t_eq(n_e, n_xi, n_xip1 , alpha_rec, alpha_recm1):
	return (alpha_rec * n_e)**-1 * ((alpha_recm1/alpha_rec)+(n_xip1/n_xi))**-1

if __name__ == '__main__':
	# Test functions. Write main code to call.
	
	# Should I implement two cases? Long and short equilibrium timescales?
	
	# Can consider only a few ions. Probably helps lots. Focus on Si XIV / S XVI first? No detectable lines nearby
	
	# Use Stingray to do the timing.

	# May need to optimise the code reasonably well.