#!/usr/bin/env python

import numpy as np
import pylab as pl
from subprocess import call
from shutil import move
from ufo_functions import *
import os

# xi_range=np.linspace(2,5,1000) ##### THIS NEEDS UPDATING TO USE SETTINGS

# densities=np.logspace(-10,-5,6)
run_settings=settings('halibut_settings.txt')
elements=run_settings.elements
densities=run_settings.densities
column=run_settings.column
lc_name=run_settings.lightcurve
xi_range=run_settings.xi_range

if not os.path.exists('pion_rates'):
	os.mkdir('pion_rates')
if not os.path.exists('pion_concs'):
	os.mkdir('pion_concs')

for density in densities:
	spex_str='spex <<EOF\ncom po\ncom pion\ncom rel 1 2\npar 1 2 hden v %s\n' % str(density)

	for xi in xi_range:
		if not os.path.exists('pion_rates/'+'rates_%s_%s.asc' % (str(density),str(xi))) or not os.path.exists('pion_concs/'+'concs_%s_%s.asc' % (str(density),str(xi))):
			spex_str+='par 1 2 xi v %s\ncalc\nasc file rates_%s_%s 1 2 rate\nasc file concs_%s_%s 1 2 icon\n' % (str(xi),str(density), str(xi),str(density), str(xi))

	spex_str+='q\nEOF\n'
	call(spex_str,shell=True)

	for xi in xi_range:
		if os.path.exists('rates_%s_%s.asc' % (str(density),str(xi))):
			move('rates_%s_%s.asc' % (str(density),str(xi)),'pion_rates/'+'rates_%s_%s.asc' % (str(density),str(xi)))
		if os.path.exists('concs_%s_%s.asc' % (str(density),str(xi))):
			move('concs_%s_%s.asc' % (str(density),str(xi)),'pion_concs/'+'concs_%s_%s.asc' % (str(density),str(xi)))
			