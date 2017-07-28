#!/usr/bin/env python

import numpy as np
import pylab as pl
from glob import glob
from astropy.io import fits

def calc_xi(n_total,r,flux):
	return flux/(n_total*r**2)

def net_rate(n_e,n_xi,n_xip1,n_xim1,alpha_rec, alpha_recm1, I_rat, I_ratm1):
	# Time dependence of ionization balance (dn_Xi/dt)
	# electron density, relative densities of ions (i, i+1, i-1), recombination coefficients, ionization rates

	# Where to get these from? Cloudy? Xstar? <- some sort of utility exists for taking rates etc from xstar, might be useful
	# SPEX has a rec_time tool. Ask CP about it.w

	recomb    = - n_xi * n_e * alpha_recm1 # recombination of Xi -> Xi-1
	recomb_p1 = n_xip1 * n_e * alpha_rec   # recombination of Xi+1 -> Xi
	ioniz     = - n_xi * I_rat             # ionization of Xi -> Xi+1
	ioniz_m1  = n_xim1 * I_ratm1           # ionization of Xi-1 -> Xi
	return recomb + recomb_p1 + ioniz + ioniz_m1

def t_eq(n_e, n_xi, n_xip1 , alpha_rec, alpha_recm1):
	return (alpha_rec * n_e)**-1 * ((alpha_recm1/alpha_rec)+(n_xip1/n_xi))**-1


class pion_rates:
	"""Big stupid class for holding data. To be tidied up nicely later, when I can be bothered."""

	def __init__(self, source='pion_rates'):
		print '\nReading ionization/recombination rates from '+source
		file_names=glob(source+'/*.asc')
		xis=[float('.'.join(x.split('.')[-3:-1]).split('_')[-1]) for x in file_names]
		self.file_names = [y for x,y in sorted(zip(xis,file_names))]
		self.xis=sorted(xis)
		self.data_stack=[]
		for file_name in self.file_names:
			self.data_stack.append(np.loadtxt(file_name,dtype='str',skiprows=3))

		### 1000 xi values, 207 ions, 6 columns (element, ion, ionization rate, recomb rate, junk, junk)
		self.data_stack=np.array(self.data_stack)
		self.elements=list(set(self.data_stack[0,:,0]))

	def filter_element(self,element):
		print '\tFiltering rates table for',element
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		return filtered_data

	def filter_ion(self,element,ion):
		print '\tFiltering rates table for',element,ion
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		ions=filtered_data[0,:,1]
		filtered_data=filtered_data[:,ions==ion,:]
		return filtered_data




class pion_concentrations:
	"""Big stupid class for holding data. To be tidied up nicely later, when I can be bothered."""

	def __init__(self, source='pion_concs'):
		print '\nReading equilibrium concentrations from '+source
		file_names=glob(source+'/*.asc')
		xis=[float('.'.join(x.split('.')[-3:-1]).split('_')[-1]) for x in file_names]
		self.file_names = [y for x,y in sorted(zip(xis,file_names))]
		self.xis=sorted(xis)
		self.data_stack=[]
		for file_name in self.file_names:
			###HAVE TO MANUALLY READ THESE FILES. ANNOYING.
			# self.data_stack.append(np.loadtxt(file_name,dtype='str',skiprows=4))
			i=0
			temp_array=[]
			for row in open(file_name,'r'):
				if i>3:
					# print row.strip()
					element= row[0:3].strip()
					# print element
					ion=row[3:10].strip()
					iontype=row[11:20].strip()
					remaining=row[21:].strip().split()
					charge, fraction, absolute = remaining[0],remaining[1],remaining[2]
					# print element,ion,iontype,charge,fraction, absolute
					temp_array.append([element,ion,iontype,charge,fraction, absolute])
				i+=1
			# exit()
			temp_array=np.array(temp_array)
			self.data_stack.append(temp_array)


		### 1000 xi values, 207 ions, 6 columns (element, ion, ionization rate, recomb rate, junk, junk)
		self.data_stack=np.array(self.data_stack)
		self.calc_elements()
		# print self.data_stack.shape
		# exit()
	def calc_elements(self):
		self.elements=[]
		for sub_array in self.data_stack:
			self.elements=self.elements+list(set(sub_array[:,0]))
			# print temp_elements
			# exit()
		self.elements=list(set(self.elements))
		# print self.elements
		# exit()
		# self.elements=list(set(self.data_stack[0,:,0]))

	def filter_element(self,element):
		print '\tFiltering concentrations table for',element
		filtered_data=[]
		for sub_array in self.data_stack:
			filtered_sub_array=sub_array[sub_array[:,0]==element]
			filtered_data.append(filtered_sub_array)
		# filtered_data=self.data_stack[:,elements==element,:]
		return np.array(filtered_data)

	def filter_ion(self,element,ion):
		print '\tFiltering concentrations table for',element,ion
		filtered_data=[]
		for sub_array in self.data_stack:
			filtered_sub_array=sub_array[sub_array[:,0]==element]
			filtered_sub_array=filtered_sub_array[filtered_sub_array[:,1]==ion]
			if len(filtered_sub_array)==1:
				filtered_data.append(filtered_sub_array[0][-1])
			else:
				filtered_data.append(0)


		# elements=self.data_stack[0,:,0]
		# filtered_data=self.data_stack[:,elements==element,:]
		# ions=filtered_data[0,:,1]
		# filtered_data=filtered_data[:,ions==ion,:]
		return np.array(filtered_data)


class lightcurve:
	"""Lightcurve class. Self explanatory. Assumes standard fits file format."""
	def __init__(self, filename):
		print '\nInitializing lightcurve:',filename
		self.filename=filename
		rate_table=fits.open(filename)['RATE'].data
		self.countrate=rate_table['RATE']
		self.error=rate_table['ERROR']
		self.time=rate_table['TIME']

	def rebin(self,factor):
		"""Function to resample the lightcurve. Probably necessary."""
		### This is very crude and could use re-writing.
		print '\tResampling lightcurve ('+self.filename+') by a factor of',factor
		self.time=self.time[::factor]
		n_timesteps=len(self.countrate)
		remainder=n_timesteps%factor
		new_countrate=[]
		for i in range(0,n_timesteps/int(factor)):
			new_countrate.append(np.sum(self.countrate[factor*i:factor*i+factor])/float(factor))
		if remainder !=0:
			new_countrate.append(np.sum(self.countrate[factor*(i+1):factor*(i+1)])/float(remainder))
		self.countrate=np.array(new_countrate)
		pass

	def filter_null(self):
		"""Remove all zeros from lightcurve"""
		print '\tRemoving zeros from lightcurve. Caution! Only run this AFTER rebinning, not before.'
		self.time=self.time[self.countrate>0]
		self.countrate=self.countrate[self.countrate>0]
		pass


		







if __name__ == '__main__':
	# Test functions. Write main code to call.
	
	# Should I implement two cases? Long and short equilibrium timescales?
	
	# Can consider only a few ions. Probably helps lots. Focus on Si XIV / S XVI first? No detectable lines nearby
	
	# Use Stingray to do the timing.

	# May need to optimise the code reasonably well.

	#print "Doesn't actually do anything yet."

	rates=pion_rates()
	concentrations=pion_concentrations()

	xi_vals=rates.xis
	fexxv_rates=rates.filter_ion('Fe','XXV')
	fexxvi_rates=rates.filter_ion('Fe','XXVI')

	# print concentrations.elements
	fexxv_concs=concentrations.filter_ion('Fe','XXV')
	fexxvi_concs=concentrations.filter_ion('Fe','XXVI')
	sixiv_concs=concentrations.filter_ion('Si','XIV')
	sxvi_concs=concentrations.filter_ion('S','XVI')
	arxviii_concs=concentrations.filter_ion('Ar','XVIII')
	caxx_concs=concentrations.filter_ion('Ca','XX')
	fexxvi_concs=concentrations.filter_ion('Fe','XXVI')
	# print fexxv_concs.shape
	# exit()


	concentrations_fig=pl.figure('Equilibrium concentrations')
	ax=pl.subplot(111)
	pl.plot(xi_vals,fexxv_concs, label='Fe XXV')
	pl.plot(xi_vals, fexxvi_concs, label='Fe XXVI')
	pl.plot(xi_vals, sixiv_concs, label='Si XIV')
	pl.plot(xi_vals, sxvi_concs, label='S XVI')
	pl.plot(xi_vals, arxviii_concs, label='Ar XVIII')
	pl.plot(xi_vals, caxx_concs, label='Ca XX')
	pl.legend()
	pl.savefig('concentrations_test.pdf',bbox_inches='tight')



	test_lightcurve=lightcurve('src_10s.lc')

	lightcurve_fig=pl.figure('Lightcurve')
	ax=pl.subplot(111)
	pl.plot(test_lightcurve.time,test_lightcurve.countrate,color='k',label='Raw lightcurve')
	test_lightcurve.rebin(10)
	test_lightcurve.filter_null()
	pl.plot(test_lightcurve.time,test_lightcurve.countrate,color='r',label='Filtered, resampled')
	pl.legend()
	pl.savefig('lightcurve_test.pdf',bbox_inches='tight')
	
	# pl.show()

	# pl.plot(xi_vals,fexxv_rates[:,0,2],label='Fe XXV Ionization')
	# pl.plot(xi_vals,fexxv_rates[:,0,3],label='Fe XXV Recombination')
	# pl.plot(xi_vals,fexxvi_rates[:,0,2],label='Fe XXVI Ionization')
	# pl.plot(xi_vals,fexxvi_rates[:,0,3],label='Fe XXVI Recombination')

	# pl.show()