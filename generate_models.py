#!/usr/bin/env python

import numpy as np
import pylab as pl
import os
from shutil import move
from ufo_functions import *
from glob import glob
from subprocess import call

def generate_spex_string(commands,normalization,outfilename):
	command_str='spex<<EOF\ncom slab\ncom po\ncom rel 2 1\n'
	command_str+='par 1 2 norm v '+str(normalization)+'\n'
	for c in commands:
		command_str+=c+'\n'
	command_str+='c\npar sh\npar sh free\nmo sh\n'
	plot_str='p de xs\np ty mo\np ux ke\np uy a\np rx 2 10\np ry 0 0.3\np x  log\np y  lin\np fil dis f\np\n'
	command_str+=plot_str
	out_str='plot adum '+str(outfilename)+' over\n'
	command_str+=out_str+'q\nEOF\n'
	return command_str
	
def tidy_spectra(stem,model_dir):
	# if not os.path.exists('model_spectra'):
		# os.mkdir('model_spectra')
	for filename in glob(stem+'.qdp'):
		print filename,model_dir+filename
		move(filename,model_dir+filename)

def main():
	run_settings=settings('halibut_settings.txt')
	elements=run_settings.elements
	densities=run_settings.densities
	column=run_settings.column
	lc_name=run_settings.lightcurve
	# spectra_dir=

	for density in densities:
		filenames=[]
		for element in elements:
			
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

		mean_cr=np.mean(lightcurve)
		### Might want to downsample again to some output resolution...


		concentrations=[np.load(x)['concentrations'] for x in filenames] # Can't convert to array because different sizes. Slightly irritating.
		

		for i,t in enumerate(times[:-1]):
			countrate=lightcurve[i]
			pow_norm=countrate/mean_cr*12.566 # I have no idea where this number came from, Ciro wrote it down once.
			commands=[]

			if not os.path.exists(run_settings.spectra_dir+'density_%s_spectrum_%s.qdp' % (str(density),str(i))) or run_settings.clobber:
				for j,element in enumerate(elements):

					# Concnetrations of ions for each element at each timestep.
					elemental_concentrations=concentrations[j][i]

					if len(elemental_concentrations)>=9:
						for k,conc in enumerate(elemental_concentrations):
							if conc>0:
								# print conc, column,float(conc)*float(column)
								commands.append('par 1 1 '+element.lower()+'0'*(2-len(str(k+1)))+str(k+1)+' v '+str(np.log10(float(conc)*float(column))))
					else:	
						for k,conc in enumerate(elemental_concentrations):
							if conc>0:
								commands.append('par 1 1 '+element.lower()+str(k+1)+' v '+str(np.log10(float(conc)*float(column))))

				spex_string=generate_spex_string(commands,pow_norm,'density_%s_spectrum_%s' % (str(density),str(i)))
				call(spex_string,shell=True)
				tidy_spectra('density_*_spectrum_*',run_settings.spectra_dir)
			else:
				print '\tFile:',run_settings.spectra_dir+'density_%s_spectrum_%s.qdp' % (str(density),str(i)),'already exists, skipping...'



if __name__ == '__main__':
	main()