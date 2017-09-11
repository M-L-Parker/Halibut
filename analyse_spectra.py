#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
from glob import glob
import stingray


def load_spectra(spectra_stem='spectrum_*.qdp'):
	toolbar_width=60
	print '\nLoading spectra...'

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


	energies, spectra=load_spectra()
	reference_lc=np.mean(spectra,axis=1) #reference (average) lightcurve

	# pl.plot(reference_lc)
	# pl.show()

	#### Make image of spectra over time. TBC.
	# pl.imshow(spectra.T,aspect=0.1)
	# pl.show()



	pass


if __name__ == '__main__':
	main()