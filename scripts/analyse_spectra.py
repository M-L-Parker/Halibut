#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from ufo_functions import *
from glob import glob
from scipy.interpolate import UnivariateSpline as spline
import stingray
import h5py
from dcfirr import dcfirr

colours=['orchid','dodgerblue','forestgreen','goldenrod','orangered','red'][::-1]
ratio=True
calc_dcf=False ### Douglas wrote the DCF script and it's crazy slow and I hate him.
stem='O_'
segment_size=1000




def main():
	run_settings=settings('halibut_settings.txt')
	elements=run_settings.elements
	densities=run_settings.densities
	lc_name=run_settings.lightcurve
	spectra_dir=run_settings.spectra_dir

	temp_filename='time_dependent_ions/ion_concs_%s_%s_%s.npz' % (''.join(lc_name.split('.')[:-1]), elements[0], str(densities[0]))

	print '\nReading lightcurve...'
	temp_file=np.load(temp_filename)
	times=temp_file['times']
	times=times-min(times)
	lightcurve=temp_file['lightcurve']
	mean_cr=np.mean(lightcurve)

	fig1=pl.figure()
	f1ax1=pl.subplot(211)
	f1ax2=pl.subplot(212)
	f1ax2.set_xlabel('time (s)')
	f1ax2.set_ylabel(r'ratio (test/reference)')
	f1ax1.set_ylabel(r'counts s$^{-1}$')
	f1ax1.set_xlim(0,max(times))
	f1ax2.set_xlim(0,max(times))

	fig2=pl.figure()
	f2ax1=pl.subplot(111)
	f2ax1.set_xscale('log')
	# f2ax
	f2ax1.set_xlabel(r'$f\ (\mathrm{Hz})$')
	f2ax1.set_ylabel('Lag (s)')

	if calc_dcf:
		fig3=pl.figure()
		f3ax1=pl.subplot(111)
		f3ax1.set_xlim(-10000,10000)
		f3ax1.set_xlabel('Lag (s)')
		f3ax1.set_ylabel('Correlation')
	

	cgs_densities=[]
	lag_peaks=[]
	for i,density in enumerate(densities):


		print '\nLoading',spectra_dir+'spectral_data_density_%s.hdf5' % str(density)
		cgs_densities.append(spex2cgs(density))
		infile=h5py.File(spectra_dir+'spectral_data_density_%s.hdf5' % str(density),'r')
		spectra_dset=infile['spectra']

		energies_dset=infile['energies']

		energy_bins=np.logspace(np.log10(0.3),np.log10(10),101)


		# print times
		test_lc=get_binned_lightcurve(energies_dset[:], (6.65,6.75), spectra_dset,times)
		# print test_lc.counts
		length=test_lc.counts.shape[0]
		# print times[1:length+1],lightcurve[:length]
		reference_lc=stingray.lightcurve.Lightcurve(times[1:length+1],lightcurve[:length])
		if ratio:
			ratio_lc=stingray.lightcurve.Lightcurve(times[1:length+1],test_lc.counts/reference_lc.counts,input_counts=False)


		if i==0:
			f1ax1.plot(reference_lc.time,reference_lc.counts,color='k',lw=1)
			np.savetxt('reference_lc',reference_lc.counts)
			np.savetxt('times',reference_lc.time)
		if ratio:
			f1ax2.plot(ratio_lc.time,ratio_lc.counts,lw=1,color=colours[i],label=density*1.e7)
		else:
			f1ax2.plot(test_lc.time,test_lc.counts,lw=1,color=colours[i],label=density*1.e7)

		# test_lc.write('test_lc_%s' % str(i+1),format='ascii')
		# ratio_lc.write('ratio_lc_%s' % str(i+1),format='ascii')
		np.savetxt('ratio_lc_%s' % str(i+1),ratio_lc.counts)

		segsize=10000
		print '\tCalculating time lag with',segsize,'s segments'
		if ratio:
			cross_spec=stingray.crossspectrum.AveragedCrossspectrum(reference_lc,ratio_lc,segment_size=segsize)
		else:
			cross_spec=stingray.crossspectrum.AveragedCrossspectrum(reference_lc,test_lc,segment_size=segsize)
		cross_spec=cross_spec.rebin_log(f=0.001)
		timelag=cross_spec.time_lag()
		f2ax1.plot(cross_spec.freq,timelag[0],lw=1,color=colours[i],label=density*1.e7)


		if calc_dcf:
			print '\tCalculating DCF using dcfirr.py'
			print "\tThis is really slow, this is Douglas' fault. We hate Douglas."
			if ratio:
				lag,cor,numf,indices,indfinal=dcfirr(ratio_lc.time, ratio_lc.counts, reference_lc.time, reference_lc.counts, minpt=100000,minlag=-10000,maxlag=10000)
			else:
				lag,cor,numf,indices,indfinal=dcfirr(test_lc.time, test_lc.counts, reference_lc.time, reference_lc.counts, minpt=100000,minlag=-10000,maxlag=10000)
			
			splinefit=spline(lag,cor*1,s=0.0001) # s is the smoothing condition,, but I don't really understand it. This value is based on trial and error.

			dummylags=np.linspace(-10000,10000,1.e6)
			dummycors=splinefit(dummylags)
			f3ax1.plot(lag,cor,lw=1,color=colours[i],label=density*1.e7)
			# f3ax1.plot(dummylags,splinefit(dummylags),color='k',lw=0.5)

			lag_peak=dummylags[dummycors==max(dummycors)]
			lag_peaks.append(lag_peak)

		#### Make covariance spectrum...
		#make_covariance_spectrum(energies_dset[:],spectra_dset,energy_bins,times,reference_lc)

		lc_list=[]
		print '\tCalculating lightcurves for covariance spectrum'
		for b1, b2 in zip(energy_bins[:-1],energy_bins[1:]):

			temp_lc=get_binned_lightcurve(energies_dset[:], (b1,b2), spectra_dset,times)
			lc_list.append(temp_lc)
			# exit()

		print '\tCalculating covariance spectrum'
		# print '\tNB. This is extremely slow, probably because of the long lightcurves.'
		cov_spec=stingray.Covariancespectrum(lc_list)

		print '\tCalculating average covariance spectrum'
		average_cov_spec=stingray.AveragedCovariancespectrum(lc_list,segment_size)
		# np.savetxt('energy_covar_%s.txt' % str(i),cov_spec.energy_covar)
		# np.savetxt('covar_density_%s.txt' % str(density),cov_spec.covar)
		np.savetxt('av_covar_density_%s_%ss.txt' % (str(density),str(segment_size)),average_cov_spec.covar)



	if calc_dcf:
		pl.legend(loc='best', title=r'$\mathrm{Density}\ (\times10^7\ \mathrm{cm^{-3}})$', frameon=False)

	# pl.show()
	print 'Saving lightcurves to '+stem+'analysed_lightcurves.pdf'
	fig1.savefig(stem+'analysed_lightcurves.pdf',bbox_inches='tight')
	print 'Saving lag-frequency to '+stem+'analysed_lagfreq.pdf'
	fig2.savefig(stem+'analysed_lagfreq.pdf',bbox_inches='tight')
	if calc_dcf:
		fig3.savefig(stem+'analysed_dcf.pdf',bbox_inches='tight')

		print 'Saving DCF peaks to '+stem+'lag_values.txt'
		np.savetxt(stem+'lag_values.txt',lag_peaks)


	pass


if __name__ == '__main__':
	main()
