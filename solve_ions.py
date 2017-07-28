#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *

elements=['Fe']

# Load lightcurve
test_lightcurve=lightcurve('example_lightcurve.lc')
test_lightcurve.cut_interval(5.856e8+128000,5.856e8+140000)
test_lightcurve.rebin(10)
test_lightcurve.filter_null()

# Mean ionization - 10^3 seems reasonable...
mean_countrate=test_lightcurve.mean
mean_xi=1000.
densities=[1.e10]

initial_countrate=test_lightcurve.countrate[0]
initial_xi=calc_xi_from_countrate(initial_countrate,mean_countrate,mean_xi)
print '\nInitial count rate:',initial_countrate
print 'Initial ionization:',initial_xi

# Load rates and equilibrium concentrations
rates=pion_rates()
concentrations=pion_concentrations()

for density in densities:
	print '\nRunning calculations for density',density
	for element in elements:
		print '\nCalculating time-dependent ion concentrations for',element

		# Find ions
		ions=rates.get_ions(element)
		
		# Get intial ion concentrations
		initial_concs=concentrations.get_concentrations(element, ions, initial_xi)
		current_concs=initial_concs
		temp_concs=np.zeros(current_concs.shape

		for time, countrate in zip(test_lightcurve.time, test_lightcurve.countrate):

			#### Some sort of convergence needed here?
			for index,ion in enumerate(ions):
				print 'Ion',ion
				if ion==min(ions):
					# if no ions in current state or state above, set no ions for next step
					if current_concs[index]==0. and current_concs[index+1]==0.:
						temp_concs[index]=0.

				if ion==max(ions):
					# if no ions in current state or state below, set no ions for next step
					if current_concs[index]==0. and current_concs[index-1]==0.:
						temp_concs[index]=0
				else:
					# if no ions in state and adjacent states, set no ions for next step
					if current_concs[index]==0. and current_concs[index-1]==0. and current_concs[index+1]==0:
						temp_concs[index]=0