#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from ufo_functions import *
import roman
from matplotlib import gridspec

colours=['dodgerblue','forestgreen','goldenrod','orangered','red'][::-1]


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
		gs=gridspec.GridSpec(4,1,hspace=0.025)
		ax1=pl.subplot(gs[0,0])
		ax2=pl.subplot(gs[1,0])
		ax3=pl.subplot(gs[2,0],sharey=ax2)
		ax4=pl.subplot(gs[3,0],sharey=ax2)

		ax1.set_xlim(0,max(zero_times))
		ax2.set_xlim(0,max(zero_times))
		ax3.set_xlim(0,max(zero_times))
		ax4.set_xlim(0,max(zero_times))

		ax1.set_ylabel(r'$\mathrm{Count\ rate\ (s^{-1})}$')
		ax4.set_xlabel('Time (s)')

		labels = [item.get_text() for item in ax1.get_xticklabels()]
		empty_string_labels = ['']*len(labels)
		ax1.set_xticklabels(empty_string_labels)
		ax2.set_xticklabels(empty_string_labels)
		ax3.set_xticklabels(empty_string_labels)
		# ax4.set_xticklabels(labels)
		ax1.get_xaxis().set_tick_params(which='both', direction='in')
		ax2.get_xaxis().set_tick_params(which='both', direction='in')
		ax3.get_xaxis().set_tick_params(which='both', direction='in')
		ax4.get_xaxis().set_tick_params(which='both', direction='inout')




		ax1.plot(zero_times,lc_spline(manual_times),color='k',lw=1)
		for i,density in enumerate(densities):
			print density
			element_filename='time_dependent_ions/ion_concs_%s_%s_%s.npz' % (''.join(lc_name.split('.')[:-1]), element, str(density))
			concentrations=np.load(element_filename)['concentrations'] #ion concentrations
			print concentrations.shape[1],'ions, plotting last 3'

			ax2.plot(zero_times[:-1],concentrations[:,-1]/np.sum(concentrations[0,:]),lw=1,color=colours[i],label=density*1.e7)
			ax3.plot(zero_times[:-1],concentrations[:,-2]/np.sum(concentrations[0,:]),lw=1,color=colours[i])
			ax4.plot(zero_times[:-1],concentrations[:,-3]/np.sum(concentrations[0,:]),lw=1,color=colours[i])

			ax2.legend(loc='best', frameon=False, ncol=2, title=r'$\mathrm{Density}\ (\times10^7\ \mathrm{cm^{-3}})$')

		# pl.show()
		ax2.set_ylabel(element+' '+roman.toRoman(concentrations.shape[1])+' abundance')
		ax3.set_ylabel(element+' '+roman.toRoman(concentrations.shape[1]-1)+' abundance')
		ax4.set_ylabel(element+' '+roman.toRoman(concentrations.shape[1]-2)+' abundance')

		pl.savefig(element+'_ion_abundances.pdf',bbox_inches='tight')