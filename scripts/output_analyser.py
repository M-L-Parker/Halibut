#!/usr/bin/env python

import numpy as np
import pylab as pl
import argparse
import roman



def parse():
	parser=argparse.ArgumentParser()
	parser.add_argument("filename", help="Name of output npz file from solve_ions.py")
	parser.add_argument("ion", help='Ion to examine')
	args=parser.parse_args()
	return args

def lc_plot(infile,plot_number=111):
	print '\tGenerating lightcurve...'
	lightcurve=infile['lightcurve']
	times=infile['times']
	times=times-min(times)
	lc_ax=pl.subplot(plot_number)
	lc_ax.set_xlabel(r'$t\ \mathrm{(s)}$')
	lc_ax.set_ylabel(r'$\mathrm{Count\ rate\ (s^{-1})}$')
	pl.plot(times,lightcurve)

def concs_plot(infile,ion,plot_number=111):
	concs_ax=pl.subplot(plot_number)
	concs_ax.set_ylabel(r'$\mathrm{Ion\ abundance}$')
	concs_ax.set_xlabel(r'$t\ \mathrm{(s)}$')
	times=infile['times']
	times=times-min(times)
	concs=infile['concentrations']
	ion_concs=concs[:,ion-1]
	ion_sum=np.sum(concs[0,:])
	if ion != 1:
		ion_concs_m1=concs[:,ion-2]
		pl.plot(times[:-1],ion_concs_m1/ion_sum,label=roman.toRoman(ion-1))
	pl.plot(times[:-1],ion_concs/ion_sum,label=roman.toRoman(ion))
	if ion != concs.shape[1]:
		ion_concs_p1=concs[:,ion]
		pl.plot(times[:-1],ion_concs_p1/ion_sum,label=roman.toRoman(ion+1))
	pl.legend(loc='upper left')


def ionization_plot(infile,ion,plot_number=111):
	ionization_ax=pl.subplot(plot_number)
	ionization_ax.set_ylabel(r'$\mathrm{Ionization\ rate\ (s^{-1})}$')
	ionization_ax.set_xlabel(r'$t\ \mathrm{(s)}$')
	times=infile['times']
	times=times-min(times)
	rates=infile['ionization_rates']
	ion_rates=rates[:,ion-1]
	# ion_sum=np.sum(rates[0,:])
	if ion != 1:
		ion_rates_m1=rates[:,ion-2]
		pl.plot(times[:-1],ion_rates_m1,label=roman.toRoman(ion-1))
	pl.plot(times[:-1],ion_rates,label=roman.toRoman(ion))
	if ion < rates.shape[1]-1:
		ion_rates_p1=rates[:,ion]
		pl.plot(times[:-1],ion_rates_p1,label=roman.toRoman(ion+1))
	pl.legend(loc='upper left')

def recombination_plot(infile,ion,plot_number=111):
	recombination_ax=pl.subplot(plot_number)
	recombination_ax.set_ylabel(r'$\mathrm{Recombination\ rate\ (s^{-1})}$')
	recombination_ax.set_xlabel(r'$t\ \mathrm{(s)}$')
	times=infile['times']
	times=times-min(times)
	rates=infile['recombination_rates']
	ion_rates=rates[:,ion-1]
	# ion_sum=np.sum(rates[0,:])
	if ion > 2:
		ion_rates_m1=rates[:,ion-2]
		pl.plot(times[:-1],ion_rates_m1,label=roman.toRoman(ion-1))
	pl.plot(times[:-1],ion_rates,label=roman.toRoman(ion))
	if ion != rates.shape[1]:
		ion_rates_p1=rates[:,ion]
		pl.plot(times[:-1],ion_rates_p1,label=roman.toRoman(ion+1))
	pl.legend(loc='upper left')

def main():
	args=parse()
	infile=np.load(args.filename)
	ion=int(args.ion)

	element = args.filename.split('.')[0].split('_')[-2]
	print '\nGenerating plots for',element,roman.toRoman(ion)+':'

	fig=pl.figure(figsize=(6,12))

	lc_plot(infile,plot_number=411)
	concs_plot(infile,ion,plot_number=412)
	ionization_plot(infile,ion,plot_number=413)
	recombination_plot(infile,ion,plot_number=414)


	# fig.suptitle(element+' '+roman.toRoman(ion))

	outfilename='output_plot_'+element+roman.toRoman(ion)+'.pdf'
	print '\tSaving as',outfilename
	fig.tight_layout()
	# fig.subplots_adjust(top=0.95)
	pl.savefig(outfilename)





if __name__ == '__main__':
	main()




