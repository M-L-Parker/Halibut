#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *

if __name__ == '__main__':

	run_settings=settings('halibut_settings.txt')
	elements=run_settings.elements
	densities=run_settings.densities
	column=run_settings.column
	lc_name=run_settings.lightcurve


	print 'Loading lightcurve with appropriate settings. Move this to a function.'
	test_lightcurve=lightcurve(lc_name)
	for interval in run_settings.cut_intervals:
		test_lightcurve.cut_interval(interval[0],interval[1])
	if run_settings.rebin is not None:
		test_lightcurve.rebin(run_settings.rebin)
	test_lightcurve.filter_null()

	# Define time resolution, resample lightcurve
	resample_factor=run_settings.resample # higher number = smaller time bins 
	lc_spline=test_lightcurve.spline()
	manual_times=np.linspace(min(test_lightcurve.time),max(test_lightcurve.time),resample_factor*len(test_lightcurve.time))
	zero_times=manual_times-min(manual_times)
	# time_resolution=float(max(manual_times)-min(manual_times))/float(len(manual_times))
	
	for element in elements:
		print 'Generating',element,'plot'
		fig=pl.figure(element+' ion abundances',figsize=(6,12))
		ax1=pl.subplot(411)
		ax2=pl.subplot(412)
		ax3=pl.subplot(413)
		ax4=pl.subplot(414)
		ax1.plot(zero_times,lc_spline(manual_times),color='k',lw=1)
		for density in densities[:-1][::-1]:
			print density
			element_filename='time_dependent_ions/ion_concs_%s_%s_%s.npz' % (''.join(lc_name.split('.')[:-1]), element, str(density))
			concentrations=np.load(element_filename)['concentrations'] #ion concentrations
			print concentrations.shape[1],'ions, plotting last 3'

			ax2.plot(zero_times[:-1],concentrations[:,-1]/np.sum(concentrations[0,:]),lw=1)
			ax3.plot(zero_times[:-1],concentrations[:,-2]/np.sum(concentrations[0,:]),lw=1)
			ax4.plot(zero_times[:-1],concentrations[:,-3]/np.sum(concentrations[0,:]),lw=1)
		# pl.show()

			# ax1.plot(concentrations[:,-1],label=concentrations.shape[1])
			# ax2.plot(concentrations[:,-2],label=concentrations.shape[1]-1)
			# ax3.plot(concentrations[:,-3],label=concentrations.shape[1]-2)
		# pl.show()
		# exit()
		pl.savefig(element+'_ion_abundances.pdf',bbox_inches='tight')