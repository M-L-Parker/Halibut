#!/usr/bin/env python

import numpy as np
import h5py
from ufo_functions import *

def load_spectra(spectra_stem='spectrum_*.qdp'):
	toolbar_width=60
	print '\nLoading spectra...'
	print "This is real slow, but I can't be bothered to fix it. Sorry."

	spectra_raw=glob(spectra_stem)
	# print spectra_raw
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

if __name__ == '__main__':
	print '\nProcessing spectra files into data tables.'


	run_settings=settings('halibut_settings.txt')
	elements=run_settings.elements
	densities=run_settings.densities
	lc_name=run_settings.lightcurve
	spectra_dir=run_settings.spectra_dir

	for density in densities:

		print '\nDensity:',density

		energies, spectra=load_spectra(spectra_dir+'density_%s_spectrum_*.qdp' % str(density))

		print '\tWriting to file:',spectra_dir+'spectral_data_density_%s.hdf5' % str(density)
		outfilename=spectra_dir+'spectral_data_density_%s.hdf5' % str(density)
		outfile=h5py.File(outfilename,'w')

		spectra_dset=outfile.create_dataset("spectra",data=spectra)
		energies_dset=outfile.create_dataset("energies",data=energies)

		outfile.flush()
		outfile.close()


