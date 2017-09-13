#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
from glob import glob
import stingray


def load_spectra(spectra_stem='spectrum_*.qdp'):
	toolbar_width=60
	print '\nLoading spectra...'
	print "This is real slow, but I can't be bothered to fix it. Sorry."

	spectra_raw=glob('model_spectra/'+spectra_stem)
	indices=[int(x.split('_')[-1].split('.')[0]) for x in spectra_raw]
	sorted_spectra=[y for x,y in sorted(zip(indices,spectra_raw))]
	
	spec_list=[]
	n_spec=len(sorted_spectra)
	for i,spectrum in enumerate(sorted_spectra):
		toolbar_update(float(i)/float(n_spec),toolbar_width)
		indata=np.loadtxt(spectrum,skiprows=1)
		if i==0:
			E=10.**indata[:,0]
		F=indata[:,3]
		spec_list.append(F)
	toolbar_update(1,toolbar_width)
	print '\nDone.',n_spec,'spectra loaded'
	spec_array=np.array(spec_list)
	return E, spec_array


def get_binned_lightcurve(energies, eband, spectra,times):
	'''Returns the lightcurve (in stingray format) for a given energy band. Eband should be an iterable object with a lower and upper limit limit'''
	### Find indices of matching energies. This is horribly inelegant.
	# print eband[0],eband[1]
	lowindex=np.where(energies>eband[0])[0][0]
	highindex=np.where(energies<eband[1])[-1][-1]
	# print lowindex, highindex
	counts=np.sum(spectra[:,lowindex:highindex+1],axis=1)
	# pl.plot(counts)
	return stingray.lightcurve.Lightcurve(times[:-1],counts,input_counts=False)




def main():
	run_settings=settings('halibut_settings.txt')
	elements=run_settings.elements
	densities=run_settings.densities
	lc_name=run_settings.lightcurve

	temp_filename='time_dependent_ions/ion_concs_%s_%s_%s.npz' % (''.join(lc_name.split('.')[:-1]), elements[0], str(densities[0]))

	print '\nReading lightcurve...'
	temp_file=np.load(temp_filename)
	times=temp_file['times']
	times=times-min(times)
	lightcurve=temp_file['lightcurve']
	mean_cr=np.mean(lightcurve)


	# energies, spectra=load_spectra()
	print '\nLoading spectra, energies from text files. Remember to turn this off and load properly.'
	energies=np.loadtxt('energies_temp.txt')
	spectra=np.loadtxt('spectra_temp.txt')

	#### Save to files so I can skip loading for dev purposes.
	# np.savetxt('energies_temp.txt',energies)
	# np.savetxt('spectra_temp.txt',spectra)



	reference_lc=stingray.lightcurve.Lightcurve(times[:-1],np.sum(spectra,axis=1),input_counts=False) #reference (total) lightcurve

	energy_bins=np.logspace(np.log10(0.3),np.log10(10),101)

	print '\nCalculating binned lightcurves...'
	print 'This needs to be added to the settings. Fix it.'
	lightcurves=[]
	for elow,ehigh in zip(energy_bins[:-1],energy_bins[1:]):
		lightcurves.append(get_binned_lightcurve(energies,(elow,ehigh),spectra,times))
	print 'Done.'


	print '\nCalculating cross spectra...'
	cross_spectra=[]
	ax=pl.subplot(111)
	ax.set_xscale('log')
	for lightcurve in lightcurves:
		cross_spectra.append(stingray.crossspectrum.Crossspectrum(lightcurve,reference_lc))
		pl.plot(cross_spectra[-1].freq,cross_spectra[-1].time_lag())
	print 'Done.'


	# print '\nCalculating lags...'
	# lags=[x.time_lag() for x in cross_spectra]
	# print 'Done.'


	pl.show()
	#### Plot reference lightcurve
	# reference_lc.plot()
	# pl.show()

	#### Make image of spectra over time. TBC.
	# pl.imshow(spectra.T,aspect=0.1)
	# pl.show()



	pass


if __name__ == '__main__':
	main()