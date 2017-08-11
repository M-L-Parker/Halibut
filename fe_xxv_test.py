#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
import roman



# Load lightcurve
test_lightcurve=lightcurve('example_lightcurve.lc')
test_lightcurve.cut_interval(5.856e8+128000,5.856e8+140000)
test_lightcurve.rebin(10)
test_lightcurve.filter_null()

# Define time resolution, resample lightcurve
resample_factor=10
lc_spline=test_lightcurve.spline()
# print len(test_lightcurve.time)
# exit()
manual_times=np.linspace(min(test_lightcurve.time),max(test_lightcurve.time),resample_factor*len(test_lightcurve.time))

# Mean ionization - 10^3 seems reasonable...
mean_countrate=test_lightcurve.mean
mean_xi=5000.
element='Fe'

### Silva et al. 2016 says constant e density, but what should it be??
### Should be HII dominated, so ne ~ nH

densities=[1.e-6]

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


	# Get intial ion concentrations
	initial_concs=concentrations.get_concentrations(element, ions, mean_xi)
	current_concs=initial_concs#*density

	# print sum(current_concs)
	# exit()

	# temp_concs=np.zeros(current_concs.shape)

	# print sum(initial_concs)
	# print sum(current_concs)
	# exit()

	time_dependent_concentrations=[]
	delta_ion=0

	i_rates=[]
	r_rates=[]

	i_rates_m1=[]
	r_rates_p1=[]

	xi_values=[]

	# for ion in ions:
	for t_step, countrate in zip(range(0,len(manual_times)), lc_spline(manual_times)):
		time=manual_times[t_step]
		# temp_concs=np.copy(current_concs)

		if t_step != len(manual_times)-1:
			delta_t=manual_times[t_step+1]-time
			# print float(time-min(manual_times))/float(max(manual_times)-min(manual_times))

			current_xi=calc_xi_from_countrate(countrate, mean_countrate, mean_xi)
			xi_values.append(current_xi)

			net_rates, temp_i_rates, temp_r_rates = rates.get_net_rates(element, np.log10(current_xi), ions, current_concs)

			current_concs=current_concs+net_rates*delta_t

			# print 'ionization:',current_xi
			#### Some sort of convergence needed here?
			# print 'Ion',element, roman.toRoman(ion)



			# # if no ions in state and adjacent states, set no ions for next step
			# if current_concs[24]==0. and current_concs[24-1]==0. and current_concs[24+1]==0:
			# 	temp_concs[24]=0
			# 	i_rate=0
			# 	r_rate=0
			# else:
			# 	ion_rates=rates.get_ion_rates(element,ion,np.log10(current_xi))
			# 	i_rate=ion_rates[0]*current_concs[24]
			# 	r_rate=ion_rates[1]*current_concs[24]
			i_rates.append(temp_i_rates[24])
			r_rates.append(temp_r_rates[24])

			# 	ion_rates_m1=rates.get_ion_rates(element,ion-1,np.log10(current_xi))
			# i_rate_m1=ion_rates_m1[0]*current_concs[24-1]
			i_rate_m1=temp_i_rates[23]
			# 	ion_rates_p1=rates.get_ion_rates(element,ion+1,np.log10(current_xi))
			# i_rate_p1=ion_rates_p1[0]*current_concs[24+1]
			r_rate_p1=temp_r_rates[25]

			i_rates_m1.append(i_rate_m1)
			r_rates_p1.append(r_rate_p1)

			# 	delta_ion = i_rate_m1+r_rate_p1-i_rate-r_rate
			# 	temp_concs[24]=current_concs[24]+delta_ion*delta_t
			# 	# print current_concs[24], delta_ion*delta_t

			# print 'Rates:', i_rate, r_rate
			# print 'Current concentration:', current_concs[24]

		# exit()
			# current_concs=temp_concs
			# print current_concs
			# raw_input()

			time_dependent_concentrations.append(current_concs)
	time_dependent_concentrations=np.array(time_dependent_concentrations)
	# print time_dependent_concentrations.shape
	# print time_dependent_concentrations.shape
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
	# ax4.set_ylim(0,1)
	# ax2.set_yscale('log')
	# pl.plot(manual_times[:-1], time_dependent_concentrations[:,26])
	# print time_dependent_concentrations[:,25-1]
	# pl.plot(manual_times[:-1], time_dependent_concentrations[:,24-1])
	pl.plot(manual_times[:-1], time_dependent_concentrations[:,24-1]/np.sum(time_dependent_concentrations[-1,:]),label='Fe XXIV')
	pl.plot(manual_times[:-1], time_dependent_concentrations[:,25-1]/np.sum(time_dependent_concentrations[-1,:]),label='Fe XXV')
	pl.plot(manual_times[:-1], time_dependent_concentrations[:,26-1]/np.sum(time_dependent_concentrations[-1,:]),label='Fe XXVI')
	# pl.plot(manual_times[:-1], time_dependent_concentrations[:,26-1])
	pl.legend()


	# pl.show()
	pl.savefig('Fe_XXV_density%s_test.pdf' % str(density) ,bbox_inches='tight')

