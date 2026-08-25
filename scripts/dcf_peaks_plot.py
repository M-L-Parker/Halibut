#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from ufo_functions import *
from scipy.stats import linregress
from matplotlib.ticker import *

stem='Fe_'
lags=np.loadtxt(stem+'lag_values.txt')

run_settings=settings('halibut_settings.txt')
densities=[spex2cgs(x) for x in run_settings.densities]

ax=pl.subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
ax.yaxis.set_major_locator(FixedLocator([10 ,30,100,300,1000,3000,10000]))
# ax.yaxis.set_major_locator(FixedLocator([0.03,0.1,0.3,1.,3.,10.]))
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.set_xlim(5.e5,1.e9)
ax.set_ylim(3.e1,1.e4)
ax.set_xlabel(r'$\mathrm{Density}\ (\mathrm{cm^{-3}})$')
ax.set_ylabel('DCF peak lag (s)')

logx=np.log10(densities)
logy=np.log10(lags)

m, c, r, p, std_err = linregress(logx,logy)
print 'log(t) =',m,'* log(nH) +',c
print 't =',10.**c,'*','nH^'+str(m)

dummy_xs=np.logspace(min(logx)-1.,max(logx)+1.)
new_ys=[10.**(m*np.log10(x)+c) for x in dummy_xs]

pl.plot(dummy_xs,new_ys,color='k',lw=1.)
pl.plot(densities,lags,marker='o',color='k',ls='none')

pl.vlines([6.66e9/2.,6.66e8/2.],1.e1,1.e4,color='r',linestyles='dashed')
pl.text(6.66e8/2., 6000, r'IRAS 13224-3809, 200$R_\mathrm{G}$',rotation=270,color='r')


pl.vlines([4.6e8/2.,4.6e7/2.],1.e1,1.e4,color='dodgerblue',linestyles='dashed')
pl.text(4.6e8/2., 6000, r'PDS 456, 20$R_\mathrm{G}$',rotation=270,color='dodgerblue')
pl.text(4.6e7/2., 6000, r'PDS 456, 200$R_\mathrm{G}$',rotation=270,color='dodgerblue')





pl.savefig(stem+'DCF_density_plot.pdf',bbox_inches='tight')