#!/usr/bin/env python

import numpy as np
import pylab as pl
from glob import glob

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

	def __init__(self, source='/home/mlparker/computing/ufo_timing/pion_rates'):
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
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		return filtered_data

	def filter_ion(self,element,ion):
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		ions=filtered_data[0,:,1]
		filtered_data=filtered_data[:,ions==ion,:]
		return filtered_data




class pion_concentrations:
	"""Big stupid class for holding data. To be tidied up nicely later, when I can be bothered."""

	def __init__(self, source='/home/mlparker/computing/ufo_timing/pion_concs'):
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
					element= row[0:2].strip()
					ion=row[3:10].strip()
					iontype=row[11:20].strip()
					remaining=row[21:].strip().split()
					charge, fraction, absolute = remaining[0],remaining[1],remaining[2]
					# print element,ion,iontype,charge,fraction, absolute
					temp_array.append([element,ion,iontype,charge,fraction, absolute])
				i+=1
			temp_array=np.array(temp_array)
			print temp_array.shape
			self.data_stack.append(temp_array)


		### 1000 xi values, 207 ions, 6 columns (element, ion, ionization rate, recomb rate, junk, junk)
		self.data_stack=np.array(self.data_stack)
		print self.data_stack.shape
		exit()
		self.elements=list(set(self.data_stack[0,:,0]))

	def filter_element(self,element):
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		return filtered_data

	def filter_ion(self,element,ion):
		elements=self.data_stack[0,:,0]
		filtered_data=self.data_stack[:,elements==element,:]
		ions=filtered_data[0,:,1]
		filtered_data=filtered_data[:,ions==ion,:]
		return filtered_data










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

	fexxv_concs=concentrations.filter_ion('Fe','XXV')

	ax=pl.subplot(111)
	ax.set_yscale('log')

	pl.plot(xi_vals,fexxv_concs[:,0,-1])

	# pl.plot(xi_vals,fexxv_rates[:,0,2],label='Fe XXV Ionization')
	# pl.plot(xi_vals,fexxv_rates[:,0,3],label='Fe XXV Recombination')
	# pl.plot(xi_vals,fexxvi_rates[:,0,2],label='Fe XXVI Ionization')
	# pl.plot(xi_vals,fexxvi_rates[:,0,3],label='Fe XXVI Recombination')

	pl.legend()
	pl.show()