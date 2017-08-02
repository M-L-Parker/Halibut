#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
import roman

elements=['Fe']

# Load lightcurve
test_lightcurve=lightcurve('example_lightcurve.lc')
test_lightcurve.cut_interval(5.856e8+128000,5.856e8+140000)
test_lightcurve.rebin(10)
test_lightcurve.filter_null()

# Mean ionization - 10^3 seems reasonable...
mean_countrate=test_lightcurve.mean
mean_xi=1000.
densities=[1.e15]

### Silva et al. 2016 says constant e density, but what should it be??
### Should be HII dominated, so ne ~ nH

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
		initial_concs=concentrations.get_concentrations(element, ions, mean_xi)
		current_concs=initial_concs*density
		temp_concs=np.zeros(current_concs.shape)

		# print sum(initial_concs)
		# print sum(current_concs)
		# exit()

		time_dependent_concentrations=[]
		delta_ion=0

		for t_step, countrate in zip(range(0,len(test_lightcurve.time)), test_lightcurve.countrate):
			time=test_lightcurve.time[t_step]
			temp_concs=np.zeros(current_concs.shape)

			if t_step != len(test_lightcurve.time)-1:
				delta_t=test_lightcurve.time[t_step+1]-time
				print float(time-min(test_lightcurve.time))/float(max(test_lightcurve.time)-min(test_lightcurve.time))

				current_xi=calc_xi_from_countrate(countrate, mean_countrate, mean_xi)
				print 'ionization:',current_xi
				#### Some sort of convergence needed here?
				for index,ion in enumerate(ions):
					print 'Ion',element, roman.toRoman(ion)
					# print ion, min(ions)
					if True:#ion==25:
						# print current_concs
						# exit()
						if ion==min(ions):
							# if no ions in current state or state above, set no ions for next step
							if current_concs[index]==0. and current_concs[index+1]==0.:
								temp_concs[index]=0.
								i_rate=0
								r_rate=0
							else:
								# find ionization rate of ion, recombination rate of ion above
								ion_rates=rates.get_ion_rates(element,ion,current_xi)
								i_rate=ion_rates[0]*current_concs[index]
								r_rate=ion_rates[1]*current_concs[index]
								ion_rates_p1=rates.get_ion_rates(element,ion+1,current_xi)
								i_rate_p1=ion_rates_p1[0]*current_concs[index+1]
								r_rate_p1=ion_rates_p1[1]*current_concs[index+1]

								delta_ion = net_rate(density,current_concs[index],current_concs[index+1],0,r_rate_p1, 0, i_rate, 0)
								temp_concs[index]=current_concs[index]+delta_ion*delta_t
								# print current_concs[index], delta_ion*delta_t
								# print i_rate, r_rate
								# exit()


						elif ion==max(ions):
							# if no ions in current state or state below, set no ions for next step
							if current_concs[index]==0. and current_concs[index-1]==0.:
								temp_concs[index]=0
								i_rate=0
								r_rate=0
							else:
								ion_rates=rates.get_ion_rates(element,ion,current_xi)
								i_rate=ion_rates[0]*current_concs[index]
								r_rate=ion_rates[1]*current_concs[index]
								ion_rates_m1=rates.get_ion_rates(element,ion-1,current_xi)
								i_rate_m1=ion_rates_m1[0]*current_concs[index-1]
								r_rate_m1=ion_rates_m1[1]*current_concs[index-1]


								delta_ion = net_rate(density,current_concs[index],0,current_concs[index-1],0, r_rate, 0, i_rate_m1)
								temp_concs[index]=current_concs[index]+delta_ion*delta_t
								# print current_concs[index], delta_ion*delta_t
						else:
							# if no ions in state and adjacent states, set no ions for next step
							if current_concs[index]==0. and current_concs[index-1]==0. and current_concs[index+1]==0:
								temp_concs[index]=0
								i_rate=0
								r_rate=0
							else:
								ion_rates=rates.get_ion_rates(element,ion,current_xi)
								i_rate=ion_rates[0]*current_concs[index]
								r_rate=ion_rates[1]*current_concs[index]
								ion_rates_m1=rates.get_ion_rates(element,ion-1,current_xi)
								i_rate_m1=ion_rates_m1[0]*current_concs[index-1]
								r_rate_m1=ion_rates_m1[1]*current_concs[index-1]
								ion_rates_p1=rates.get_ion_rates(element,ion+1,current_xi)
								i_rate_p1=ion_rates_p1[0]*current_concs[index+1]
								r_rate_p1=ion_rates_p1[1]*current_concs[index+1]

								delta_ion = net_rate(density,current_concs[index],current_concs[index+1],current_concs[index-1],r_rate_m1, r_rate, i_rate, i_rate_m1)
								temp_concs[index]=current_concs[index]+delta_ion*delta_t
							# print current_concs[index], delta_ion*delta_t

						print 'Rates:', i_rate, r_rate
						print 'Current concentration:', current_concs[index]
				
			# exit()
				current_concs=temp_concs
				# print current_concs
				# raw_input()

				time_dependent_concentrations.append(temp_concs)
		time_dependent_concentrations=np.array(time_dependent_concentrations)
		# print time_dependent_concentrations.shape
		# print time_dependent_concentrations.shape
		ax1=pl.subplot(211)
		pl.plot(test_lightcurve.time,test_lightcurve.countrate)

		ax2=pl.subplot(212)
		ax2.set_yscale('log')
		# pl.plot(test_lightcurve.time[:-1], time_dependent_concentrations[:,26])
		print time_dependent_concentrations[:,25-1]
		# pl.plot(test_lightcurve.time[:-1], time_dependent_concentrations[:,24-1])
		pl.plot(test_lightcurve.time[:-1], time_dependent_concentrations[:,25-1])
		# pl.plot(test_lightcurve.time[:-1], time_dependent_concentrations[:,26-1])
		pl.show()


			# exit()




