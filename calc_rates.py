#!/usr/bin/env python

import numpy as np
import pylab as pl
from subprocess import call


xi_range=np.linspace(2,5,1000)
densities=np.logspace(-10,-5,6)

for density in densities:
	for xi in xi_range:
		call('spex <<EOF\n\
				com po\n\
				com pion\n\
				com rel 1 2\n\
				par 1 2 xi v %s\n\
				par 1 2 hden v %s\n\
				calc\n\
				asc file rates_%s_%s 1 2 rate\n\
				asc file concs_%s_%s 1 2 icon\n\
				' % (str(xi),str(density),str(density), str(xi),str(density), str(xi)), shell=True)