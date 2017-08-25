# ufo-timing



## Sequence:
1. Run calc_rates.py to generate tables of equilibrium concentrations and ionization rates. This just calls SPEX a bunch of times, very badly.
  * Run fe_xxv_test.py to check these are working. This runs the caluclations for Fe only, and produces some diagnostic plots.
2. Run the full ion solver, solve_ions.py. This calculates the time-dependent ionization rates and ion concentrations for a given lightcurve, and saves them to .npz files
  * Use output_analyser.py to view these output files, check that everything looks reasonable. At higher densities it may be necessary to sample the lightcurve more frequently (use the resample_factor parameter in solve_ions.py) to increase the time resolution.
3. Run generate_models.py to simulate spectra based on the output ion concentrations. [NOT YET IMPLEMENTED]

## To Do:
* Implement simulator to reconstruct spectra based on the time-dependent ion concentrations
* Implement time lag code to analyse the reconstructed spectra
* Add some sort of settings file, to be used by the different scripts
