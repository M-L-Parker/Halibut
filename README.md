# Halibut
### Because X-ray software is fish now.

This is a series of scripts for estimating the X-ray timing/variability properties of ionized gas responding to an X-ray source. 

In brief, it takes a lightcurve and calculates the ion abundances and corresponding rates at each time step (allowing the lightcurve to be sampled at an arbitrary frequency). These numbers are then used to calculate the spectrum of absorbing gas at each step, which can then in turn be used as input to timing/variability analysis.

Based on the method described by Silva, Uttley & Constantini (2016).

Ion concentrations and ionization/recombination rates are calculated using SPEX (Kaastra, Mewe & Nieuwenhuijzen, 1996).

## Requirements:
* Python 2.6/7
* SPEX: https://www.sron.nl/astrophysics-spex
* astropy: http://www.astropy.org/
* Stingray: http://stingraysoftware.github.io/
* h5py: http://www.h5py.org/

## Sequence:
1. Run calc_rates.py to generate tables of equilibrium concentrations and ionization rates. This just calls SPEX a bunch of times, very badly.
	Run fe_xxv_test.py to check these are working. This runs the calculations for Fe only, and produces some diagnostic plots.
2. Run the full ion solver, solve_ions.py. This calculates the time-dependent ionization rates and ion concentrations for a given lightcurve, and saves them to .npz files
	Use output_analyser.py to view these output files, check that everything looks reasonable. At higher densities it may be necessary to sample the lightcurve more frequently (use the resample_factor parameter in solve_ions.py) to increase the time resolution.
3. Run generate_models.py to simulate spectra based on the output ion concentrations.
4. Get lag products, using Stingray [NOT YET IMPLEMENTED]

## To Do:
* ~~Implement simulator to reconstruct spectra based on the time-dependent ion concentrations~~
* Implement time lag code to analyse the reconstructed spectra
* ~~Add some sort of settings file, to be used by the different scripts~~
* Come up with a better name.
* ~~Correct calc_rates.py to tidy up after itself.~~
* Write a paper about this.
