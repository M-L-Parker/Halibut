#!/usr/bin/env python

import numpy as np
import pylab as pl
from glob import glob
from astropy.io import fits
import roman
from scipy.interpolate import UnivariateSpline as spline


def calc_xi(n_total,r,flux):
	return flux/(n_total*r**2)

def calc_xi_from_countrate(countrate,meanrate,meanxi):
	# print countrate, meanrate, meanxi
	return meanxi*countrate/meanrate


# def net_rate(n_e,n_xi,n_xip1,n_xim1,alpha_rec, alpha_recm1, I_rat, I_ratm1):
# 	# Time dependence of ionization balance (dn_Xi/dt)
# 	# electron density, relative densities of ions (i, i+1, i-1), recombination coefficients, ionization rates
# 	recomb    = - n_xi * alpha_recm1 # recombination of Xi -> Xi-1
# 	recomb_p1 = n_xip1 * alpha_rec   # recombination of Xi+1 -> Xi
# 	ioniz     = - n_xi * I_rat             # ionization of Xi -> Xi+1
# 	ioniz_m1  = n_xim1 * I_ratm1           # ionization of Xi-1 -> Xi
# 	# print 'recombination rate:',recomb
# 	# print 'recombination rate of above ion:', recomb_p1
# 	# print 'ionization rate:', ioniz
# 	# print 'ionization rate of below ion:', ioniz_m1
# 	# exit()
# 	return recomb + recomb_p1 + ioniz + ioniz_m1

def t_eq(n_e, n_xi, n_xip1 , alpha_rec, alpha_recm1):
	return (alpha_rec * n_e)**-1 * ((alpha_recm1/alpha_rec)+(n_xip1/n_xi))**-1

class pion_rates:
	"""Big stupid class for holding data. To be tidied up nicely later, when I can be bothered."""

	def __init__(self, source='pion_rates', density=1.e-10):
		print '\nReading ionization/recombination rates from '+source
		file_names=glob(source+'/*%s*.asc' % str(density))
		xis=[float('.'.join(x.split('.')[-3:-1]).split('_')[-1]) for x in file_names]

		self.file_names = [y for x,y in sorted(zip(xis,file_names))]
		self.xis=sorted(xis)
		self.data_stack=[]
		for file_name in self.file_names:
			self.data_stack.append(np.loadtxt(file_name,dtype='str',skiprows=3))

		### 1000 xi values, 207 ions, 6 columns (element, ion, ionization rate, recomb rate, junk, junk)
		self.data_stack=np.array(self.data_stack)
		self.elements=list(set(self.data_stack[0,:,0]))
		self.ions=list(set(self.data_stack[0,:,1]))
		# self.define_spline_array()

	def filter_element(self,element):
		# print '\tFiltering rates table for',element
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		return filtered_data

	def filter_ion(self,element,ion):
		# print '\tFiltering rates table for',element,ion
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		ions=filtered_data[0,:,1]
		filtered_data=filtered_data[:,ions==ion,:]
		return filtered_data

	def get_ions(self,element):
		print '\tCalculating ions of',element
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		ions=sorted([roman.fromRoman(x) for x in list(set(filtered_data[0,:,1]))])
		return np.array(ions)

	def get_splines(self, element, ion):
		# print element, ion
		filtered_sub_array = self.filter_ion(element,roman.toRoman(ion))
		ionization_spline = spline(self.xis,filtered_sub_array[:,0,2],s=0)
		recomb_spline = spline(self.xis,filtered_sub_array[:,0,3],s=0)
		return ionization_spline,recomb_spline

	def get_ion_rates(self,element,ion,xi):
		i_spline, r_spline=self.get_splines(element,ion)
		return i_spline(xi), r_spline(xi)

	def get_net_rates(self,element,xi,ions,concentrations):
		r_rates=[]
		i_rates=[]
		net_rates=[]
		for ion,concentration in zip(ions,np.nan_to_num(concentrations)):
			temp_rates=self.get_ion_rates(element,ion,xi)
			# print temp_rates[1], concentration
			r_rates.append(temp_rates[1]*max(concentration,0.))
			i_rates.append(temp_rates[0]*max(concentration,0.))
		for i,ion in enumerate(ions):
			if i==0:
				net_rate=r_rates[i+1]-r_rates[i]-i_rates[i]
			elif ion==ions[-1]:
				net_rate=i_rates[i-1]-r_rates[i]-i_rates[i]
			else:
				net_rate=i_rates[i-1]+r_rates[i+1]-r_rates[i]-i_rates[i]
			net_rates.append(net_rate)
			# print net_rate, concentration
		return np.array(net_rates), np.array(i_rates), np.array(r_rates)

class pion_concentrations:
	"""Big stupid class for holding data. To be tidied up nicely later, when I can be bothered."""

	def __init__(self, source='pion_concs', density=1.e-10):
		print '\nReading equilibrium concentrations from '+source
		file_names=glob(source+'/*%s*.asc' % str(density))
		xis=[float('.'.join(x.split('.')[-3:-1]).split('_')[-1]) for x in file_names]
		self.file_names = [y for x,y in sorted(zip(xis,file_names))]
		self.xis=sorted(xis)
		self.data_stack=[]
		for file_name in self.file_names:
			###HAVE TO MANUALLY READ THESE FILES. ANNOYING.
			###The array format is different for concentrations than rates - 
			###only includes ions that are present, whereas rates includes all
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
					temp_array.append([element,ion,iontype,charge,fraction, absolute])
				i+=1
			temp_array=np.array(temp_array)
			self.data_stack.append(temp_array)

		self.data_stack=np.array(self.data_stack)
		self.calc_elements()

	def calc_elements(self):
		self.elements=[]
		for sub_array in self.data_stack:
			self.elements=self.elements+list(set(sub_array[:,0]))
		self.elements=list(set(self.elements))

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

		return np.array(filtered_data,dtype=float)

	def get_concentrations(self, element, ions, xi):
		if not hasattr(ions,'__iter__'):
			ions=[ions]
		concentrations=[]
		for ion in ions:
			print '\nGetting equilibrium concentration of',element,roman.toRoman(ion)
			filtered_array=self.filter_ion(element,roman.toRoman(ion))
			# print filtered_array
			if sum(filtered_array)==0:
				concentrations.append(0.)
			else:
				ion_spline=spline(self.xis,filtered_array,s=0)
				estimate=float(ion_spline(np.log10(xi)))
				if estimate>1.e-20:
					concentrations.append(estimate)
				else:
					concentrations.append(0.)

		return np.array(concentrations)




class lightcurve:
	"""Lightcurve class. Self explanatory. Assumes standard fits file format."""
	def __init__(self, filename):
		print '\nInitializing lightcurve:',filename
		self.filename=filename
		rate_table=fits.open(filename)['RATE'].data
		self.countrate=rate_table['RATE']
		self.error=rate_table['ERROR']
		self.time=rate_table['TIME']
		self.mean=np.mean(self.countrate)
		self.filtered=False

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
		self.mean=np.mean(self.countrate)
		pass

	def filter_null(self):
		"""Remove all zeros from lightcurve"""
		print '\tRemoving zeros from lightcurve. Caution! Only run this AFTER rebinning, not before.'
		self.time=self.time[self.countrate>0]
		self.countrate=self.countrate[self.countrate>0]
		self.filtered=True
		pass

	def spline(self):
		"""Fit a splint to lightcurve"""
		if not self.filtered:
			print 'WARNING: Lightcurve should be filtered for zeros (lighcurve.filter_null()) before fitting spline'
		lc_spline=spline(self.time,self.countrate,s=0)
		return lc_spline

	def cut_interval(self,tmin,tmax):
		# find indices of intervals to keep
		print '\tCutting times between',5.856e8+130000,'and',5.856e8+140000
		indices=[]
		for i,t in enumerate(self.time):
			if t<tmin or t>tmax:
				indices.append(i)
		self.time=self.time[indices]
		self.countrate=self.countrate[indices]
		self.error=self.error[indices]
		self.mean=np.mean(self.countrate)




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
	sxvi_rates=rates.filter_ion('S','XVI')
	sixiv_rates=rates.filter_ion('Si','XIV')
	rates_fig=pl.figure('Ionization rates')
	ax=pl.subplot(111)
	ax.set_yscale('log')
	pl.plot(xi_vals,fexxvi_rates[:,0,2],color='dodgerblue',label='Fe XXVI')
	pl.plot(xi_vals,fexxvi_rates[:,0,3],color='dodgerblue',ls='--')
	pl.plot(xi_vals,fexxv_rates[:,0,2],color='red', label='Fe XXV')
	pl.plot(xi_vals,fexxv_rates[:,0,3],color='red',ls='--')
	pl.plot(xi_vals,sxvi_rates[:,0,2],color='goldenrod', label='S XVI')
	pl.plot(xi_vals,sxvi_rates[:,0,3],color='goldenrod',ls='--')
	pl.plot(xi_vals,sixiv_rates[:,0,2],color='forestgreen', label='Si XIV')
	pl.plot(xi_vals,sixiv_rates[:,0,3],color='forestgreen',ls='--')
	pl.legend()
	pl.savefig('rates_test.pdf',bbox_inches='tight')



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
	pl.plot(xi_vals,fexxv_concs, label='Fe XXV', color='r')
	pl.plot(xi_vals, fexxvi_concs, label='Fe XXVI', color='dodgerblue')
	pl.plot(xi_vals, sixiv_concs, label='Si XIV', color='forestgreen')
	pl.plot(xi_vals, sxvi_concs, label='S XVI', color='goldenrod')
	pl.plot(xi_vals, arxviii_concs, label='Ar XVIII', color='orchid')
	pl.plot(xi_vals, caxx_concs, label='Ca XX', color='navy')
	pl.legend()
	pl.savefig('concentrations_test.pdf',bbox_inches='tight')



	test_lightcurve=lightcurve('example_lightcurve.lc')

	lightcurve_fig=pl.figure('Lightcurve',figsize=(10,5))
	ax=pl.subplot(111)
	ax.set_xlim(min(test_lightcurve.time),max(test_lightcurve.time))
	pl.plot(test_lightcurve.time,test_lightcurve.countrate,color='k',label='Raw lightcurve')
	test_lightcurve.cut_interval(5.856e8+128000,5.856e8+140000)
	test_lightcurve.rebin(10)
	test_lightcurve.filter_null()
	pl.plot(test_lightcurve.time,test_lightcurve.countrate,color='r',label='Filtered, resampled')
	spline_test=test_lightcurve.spline()
	dummy_times=np.linspace(min(test_lightcurve.time),max(test_lightcurve.time),100000)
	pl.plot(dummy_times,spline_test(dummy_times),color='dodgerblue',label='Spline')
	# pl.plot()
	pl.legend()
	pl.savefig('lightcurve_test.pdf',bbox_inches='tight')
	
	# pl.show()

	# pl.plot(xi_vals,fexxv_rates[:,0,2],label='Fe XXV Ionization')
	# pl.plot(xi_vals,fexxv_rates[:,0,3],label='Fe XXV Recombination')
	# pl.plot(xi_vals,fexxvi_rates[:,0,2],label='Fe XXVI Ionization')
	# pl.plot(xi_vals,fexxvi_rates[:,0,3],label='Fe XXVI Recombination')

	# pl.show()