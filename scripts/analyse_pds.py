#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from dcfirr import *
import matplotlib.pyplot as pl


obsids=['A','B','C','D','E']
cors=[]
for obs in obsids:

	reference_fits=fits.open('obs%s_2000_10000_src.fits' % obs)
	test_fits=fits.open('obs%s_8000_10000_src.fits' % obs)

	nonzero=reference_fits['RATE'].data['RATE']>0.

	ref_times=reference_fits['RATE'].data['TIME'][nonzero]-min(reference_fits['RATE'].data['TIME'][nonzero])
	ref_counts=reference_fits['RATE'].data['RATE'][nonzero]*300
	test_counts=test_fits['RATE'].data['RATE'][nonzero]*300

	ratio_lc=test_fits['RATE'].data['RATE'][nonzero]/reference_fits['RATE'].data['RATE'][nonzero]*300

	lag,cor,numf,indices,indfinal=dcfirr(ref_times,\
	                                     ratio_lc,\
	                                     ref_times,\
	                                     ref_counts,\
	                                     minpt=2000,\
	                                     minlag=-10000,\
	                                     maxlag=10000)

	# fig2=pl.figure()

	cors.append(list(cor))
	pl.plot(lag,cor)

# cors=np.array(cors)
# print cors

# pl.plot(np.mean(cors,axis=0),color='k')
pl.show()