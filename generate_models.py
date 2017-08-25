#!/usr/bin/env python

import numpy as np
import pylab as pl
from ufo_functions import *
from glob import glob
import os


##### Need to make a settings file to hold this stuff! This is impractical. Too many scripts.
elements=['Fe','Si','S','Ne', 'Mg', 'Ar','Ca']
densities=[1e-7] #Units of 10^20 m^-3
lc_name='example_lightcurve.lc'


def get_lightcurve():
	pass

def main():

	### load settings - to be implemented

	for density in densities:
		filenames=[]
		for element in elements:
			
			### Find 
			element_filename='time_dependent_ions/ion_concs_%s_%s_%s.npz' % (''.join(lc_name.split('.')[:-1]), element, str(density))
			if not os.path.exists(element_filename):
				print 'ERROR: file',element_filename,'not found. Exiting.'
				exit()
			filenames.append(element_filename)

		filenames=np.array(filenames)

		temp_file=np.load(filenames[0])
		times=temp_file['times']
		times=times-min(times)
		lightcurve=temp_file['lightcurve']
		### Might want to downsample again to some output resolution...


		concentrations=[np.load(x)['concentrations'] for x in filenames] # Can't convert to array because different sizes. Slightly irritating.
		
		# print files

		for i,t in enumerate(times[:-1]):
			countrate=lightcurve[i]
			for j,element in enumerate(elements):

				# Concnetrations of ions for each element at each timestep.
				elemental_concentrations=concentrations[j][i]
				# Need to apply to model, save spectrum (don't simulate, waste of time)





if __name__ == '__main__':
	main()