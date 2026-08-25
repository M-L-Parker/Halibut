#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
import roman
import os



# Load lightcurve
lc_filename='example_lightcurve.lc'
test_lightcurve=lightcurve(lc_filename)
test_lightcurve.cut_interval(5.856e8+128000,5.856e8+140000)
test_lightcurve.rebin(10) # higher number = bigger time bins
test_lightcurve.filter_null()

# Define time resolution, resample lightcurve
resample_factor=1 # higher number = smaller time bins 
lc_spline=test_lightcurve.spline()
manual_times=np.linspace(min(test_lightcurve.time),max(test_lightcurve.time),resample_factor*len(test_lightcurve.time))
time_resolution=float(max(manual_times)-min(manual_times))/float(len(manual_times))

# Mean ionization - 10^3 seems reasonable...
mean_countrate=test_lightcurve.mean
mean_xi=5000.
element='Fe'


### Note that these are spex densities! hden is in units of 10^20/m^3.
### Silva et al. (2016) use 10^7/cm^3 as representative value, corresponds to 10^-7 in these units.
densities=[1.e-7]


toolbar_width=60

for density in densities:

	initial_countrate=test_lightcurve.countrate[0]
	initial_xi=calc_xi_from_countrate(initial_countrate,mean_countrate,mean_xi)
	print '\nInitial count rate:',initial_countrate
	print 'Initial ionization:',initial_xi

	# Load rates and equilibrium concentrations
	rates=pion_rates(density=density)
	ions=rates.get_ions(element)
	concentrations=pion_concentrations(density=density)

	true_density=density*1.e20
	print '\nCalculating time-dependent ion concentrations for:'
	print '\tDensity:',true_density,'m^-3'
	print '\tElement:',element
	print 'Sampling every',time_resolution,'seconds:'

	# Get intial ion concentrations
	initial_concs=concentrations.get_concentrations(element, ions, initial_xi)
	current_concs=initial_concs

	time_dependent_concentrations=[]
	delta_ion=0

	i_rates=[]
	r_rates=[]

	i_rates_m1=[]
	r_rates_p1=[]

	xi_values=[]


	for t_step, countrate in zip(range(0,len(manual_times)), lc_spline(manual_times)):
		time=manual_times[t_step]

		### This needs modifying to not run every step. I can't be bothered now.
		toolbar_update(float(t_step)/float(len(manual_times)),toolbar_width)

		if t_step != len(manual_times)-1:
			delta_t=manual_times[t_step+1]-time

			current_xi=calc_xi_from_countrate(countrate, mean_countrate, mean_xi)
			xi_values.append(current_xi)

			net_rates, temp_i_rates, temp_r_rates = rates.get_net_rates(element, np.log10(current_xi), ions, current_concs)

			current_concs=current_concs+net_rates*delta_t

			i_rates.append(temp_i_rates[24])
			r_rates.append(temp_r_rates[24])
			i_rate_m1=temp_i_rates[23]
			r_rate_p1=temp_r_rates[25]

			i_rates_m1.append(i_rate_m1)
			r_rates_p1.append(r_rate_p1)

			time_dependent_concentrations.append(current_concs)
	toolbar_update(1,toolbar_width)

	print '\nDone.'

	# Final time-dependent ion concentrations array for element. Axis 0 is time, axis 1 ion number
	time_dependent_concentrations=np.array(time_dependent_concentrations)

	print '\nPlotting test figure:'

	fig=pl.figure(figsize=(6,12))

	ax1=pl.subplot(411)
	ax1.set_ylabel(r'$\xi$')
	pl.plot(manual_times[:-1],xi_values)


	ax2=pl.subplot(412)
	ax2.set_ylabel(r'$\mathrm{Rate\ (s^{-1})}$')
	pl.plot(manual_times[:-1],i_rates,label='Ionization')
	pl.plot(manual_times[:-1],r_rates,label='Recombination')
	pl.plot(manual_times[:-1],[x+y for x,y in zip(i_rates, r_rates)],label='Net loss')
	pl.legend()

	ax3=pl.subplot(413)
	ax3.set_ylabel(r'$\mathrm{Rate\ (s^{-1})}$')
	pl.plot(manual_times[:-1],i_rates_m1,label='Fe XXIV Ionization')
	pl.plot(manual_times[:-1],r_rates_p1,label='Fe XXVI Recombination')
	pl.plot(manual_times[:-1],[x+y for x,y in zip(i_rates_m1, r_rates_p1)],label='Net gain')
	pl.legend()


	ax4=pl.subplot(414)
	pl.plot(manual_times[:-1], time_dependent_concentrations[:,24-1]/np.sum(time_dependent_concentrations[-1,:]),label='Fe XXIV')
	pl.plot(manual_times[:-1], time_dependent_concentrations[:,25-1]/np.sum(time_dependent_concentrations[-1,:]),label='Fe XXV')
	pl.plot(manual_times[:-1], time_dependent_concentrations[:,26-1]/np.sum(time_dependent_concentrations[-1,:]),label='Fe XXVI')
	pl.legend()


	# pl.show()
	pl.savefig('Fe_XXV_density%s_test.pdf' % str(density) ,bbox_inches='tight')
	print '\tSaved as Fe_XXV_density%s_test.pdf' % str(density)

	print '\nSaving ion concentrations:'
	output_dir='time_dependent_ions'
	lightcurve_stem=''.join(lc_filename.split('.')[:-1])
	outfilename='ion_concs_'+lightcurve_stem+'_'+element+'_'+str(density)+'.npz'
	if not os.path.exists(output_dir):
		print '\tPath',output_dir,'does not exist, making folder'
		os.mkdir(output_dir)
	if os.path.exists(output_dir+'/'+outfilename):
		print '\tFile',outfilename,'already exists, deleting'
		os.remove(output_dir+'/'+outfilename)
	np.savez(output_dir+'/'+outfilename, times=manual_times, concentrations='time_dependent_concentrations',\
			ionizations=xi_values, lightcurve=lc_spline(manual_times))
	print '\tSaved as',outfilename

