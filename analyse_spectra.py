#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
from glob import glob
import stingray
import h5py
from dcfirr import dcfirr

colours=['dodgerblue','forestgreen','goldenrod','orangered','red'][::-1]
ratio=True


def get_binned_lightcurve(energies, eband, spectra_dset,times):
	'''Returns the lightcurve (in stingray format) for a given energy band. Eband should be an iterable object with a lower and upper limit limit'''
	### Find indices of matching energies. This is horribly inelegant.
	# print eband[0],eband[1]
	lowindex=np.where(energies>eband[0])[0][0]
	highindex=np.where(energies<eband[1])[-1][-1]
	# print lowindex, highindex
	counts=np.sum(spectra_dset[:,lowindex:highindex+1],axis=1)
	# pl.plot(counts)
	return stingray.lightcurve.Lightcurve(times[1:counts.shape[0]+1],counts,input_counts=True)




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

	# print times.shape


	fig1=pl.figure()
	f1ax1=pl.subplot(211)
	f1ax2=pl.subplot(212)
	f1ax2.set_xlabel('time (s)')
	f1ax2.set_ylabel(r'ratio (test/reference)')
	f1ax1.set_ylabel(r'counts s$^{-1}$')
	f1ax1.set_xlim(0,max(times))
	f1ax2.set_xlim(0,max(times))
	# f1ax2.set_ylim(0.42,0.55)

	fig2=pl.figure()
	f2ax1=pl.subplot(111)
	f2ax1.set_xscale('log')
	# f2ax
	f2ax1.set_xlabel(r'$f\ (\mathrm{Hz})$')
	f2ax1.set_ylabel('Lag (s)')

	fig3=pl.figure()
	f3ax1=pl.subplot(111)
	f3ax1.set_xlim(-10000,10000)
	f3ax1.set_xlabel('Lag (s)')
	f3ax1.set_ylabel('Correlation')
	# f3ax1.set_xscale('log')
	# f3ax1.set_xlim(1)

	for i,density in enumerate(densities):


		print '\nLoading',spectra_dir+'spectral_data_density_%s.hdf5' % str(density)
		infile=h5py.File(spectra_dir+'spectral_data_density_%s.hdf5' % str(density),'r')
		spectra_dset=infile['spectra']
		energies_dset=infile['energies']

		energy_bins=np.logspace(np.log10(0.3),np.log10(10),101)


		test_lc=get_binned_lightcurve(energies_dset[:], (6.6,6.8), spectra_dset,times)
		length=test_lc.counts.shape[0]
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

		segsize=20000
		if ratio:
			cross_spec=stingray.crossspectrum.AveragedCrossspectrum(ratio_lc,reference_lc,segment_size=20000)
		else:
			cross_spec=stingray.crossspectrum.AveragedCrossspectrum(test_lc,reference_lc,segment_size=20000)
		# cross_spec.rebin_log(f=0.1)
		timelag=cross_spec.time_lag()
		f2ax1.errorbar(cross_spec.freq,timelag[0],timelag[1],lw=1,color=colours[i],label=density*1.e7)




		#### Plot reference lightcurve
		# reference_lc.plot()
		# pl.show()

		#### Make image of spectra over time. TBC.
		# pl.imshow(spectra.T,aspect=0.1)
		# pl.show()
		print '\nCalculating DCF using dcfirr.py'
		print "This is really slow, this is Douglas' fault. We hate Douglas."
		# if i==0:
		if ratio:
			lag,cor,numf,indices,indfinal=dcfirr(ratio_lc.time, ratio_lc.counts, reference_lc.time, reference_lc.counts, minpt=100000,minlag=-10000,maxlag=10000)
		else:
			lag,cor,numf,indices,indfinal=dcfirr(test_lc.time, test_lc.counts, reference_lc.time, reference_lc.counts, minpt=100000,minlag=-10000,maxlag=10000)
		# else:
			# lag,cor,numf,indices,indfinal=dcfirr(ratio_lc.time, ratio_lc.counts, reference_lc.time, reference_lc.counts, minpt=100000,minlag=-30000,maxlag=30000,indices=indices,indfinal=indfinal,numf=numf)
		f3ax1.plot(lag,cor,lw=1,color=colours[i],label=density*1.e7)

	pl.legend(loc='best', title=r'$\mathrm{Density}\ (\times10^7\ \mathrm{cm^{-3}})$', frameon=False)

	# pl.show()
	fig1.savefig('analysed_lightcurves.pdf',bbox_inches='tight')
	fig2.savefig('analysed_lagfreq.pdf',bbox_inches='tight')
	fig3.savefig('analysed_dcf.pdf',bbox_inches='tight')


	pass


if __name__ == '__main__':
	main()