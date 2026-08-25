#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from ufo_functions import *
from glob import glob
import stingray
import h5py



def main():
	timestep=10.
	print "WARNING: xs_variance.py is assuming delta-t=10s for all lightcurves, because I'm lazy"

	run_settings=settings('halibut_settings.txt')
	elements=run_settings.elements
	densities=run_settings.densities
	lc_name=run_settings.lightcurve
	spectra_dir=run_settings.spectra_dir

	ebins=np.logspace(np.log10(0.3),np.log10(10),1001)
	mean_es=[(x+y)/2. for x, y in zip(ebins[:-1],ebins[1:])]


	for density in densities:
		print density
		indata=h5py.File('/Users/mlparker/data/ufo_timing_data/spectral_data_density_%s.hdf5' % str(density))
		spectra_dset=indata['spectra']
		energies_dset=indata['energies']

		times=[10.* i for i in range(0,len(spectra_dset[:,0]))]

		ref_lc=stingray.lightcurve.Lightcurve(counts=get_binned_lightcurve(energies_dset[:],(0.3,10),spectra_dset),time=times)

		rms_spectrum=[]
		coherence_spectrum=[]

		for elow, ehigh in zip(ebins[:-1],ebins[1:]):
		    temp_lc=stingray.lightcurve.Lightcurve(counts=get_binned_lightcurve(energies_dset[:],(elow,ehigh),spectra_dset)*100,time=times)

		    # print temp_lc.time.shape
		    # print ref_lc.time.shape
		    # pl.plot(temp_lc.time,temp_lc.counts)
		    # pl.show()
		    # exit()
		    # print temp_lc.counts,temp_lc.time
		    xs_var=stingray.excess_variance(temp_lc,'fvar')
		    # print xs_var
		    # exit()
		    coherence=stingray.coherence(temp_lc,ref_lc)
		    rms_spectrum.append(xs_var[0])
		    coherence_spectrum.append(coherence)

		np.savetxt('rms_spectrum_density_%s.dat' % str(density), rms_spectrum)
		np.savetxt('coherence_spectrum_density_%s.dat' % str(density), coherence_spectrum)



	pass

if __name__ == '__main__':
	main()