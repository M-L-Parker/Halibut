#!/usr/bin/env python

from ufo_functions import *
from scipy.interpolate import UnivariateSpline as spline
import os


def main():
    run_settings = settings('halibut_settings.txt')
    elements = run_settings.elements
    densities = run_settings.densities
    lc_name = run_settings.lightcurve
    spectra_dir = run_settings.spectra_dir

    temp_filename = 'time_dependent_ions/ion_concs_%s_%s_%s.npz' % (
        ''.join(lc_name.split('.')[:-1]), elements[0], str(densities[0]))

    print '\nReading lightcurve...'
    temp_file = np.load(temp_filename)
    times = temp_file['times']
    times = times - min(times)
    lightcurve = temp_file['lightcurve']
    mean_cr = np.mean(lightcurve)
    sorted_lc = sorted(lightcurve)

    print 'Calculating cumulative flux distribution...'
    cumulative = [sum(sorted_lc[:i]) / sum(lightcurve)
                  for i in range(0, len(lightcurve))]
    cumulative_spline = spline(sorted_lc, cumulative, s=0.1)

    root_splines = [spline(sorted_lc, cumulative - x, s=0.1)
                    for x in np.linspace(0.1, 0.9, 9)]
    roots = [s.roots()[0] for s in root_splines]

    fig1 = pl.figure(figsize=(12, 6))

    ax1 = pl.subplot(121)
    ax1.set_xlim(0, max(times) / 1000)
    ax1.set_ylim(0, max(lightcurve))
    ax1.plot(times / 1000, lightcurve)
    ax1.hlines(roots, 0, max(times) / 1000, lw=0.5)
    # print roots

    ax2 = pl.subplot(122)
    ax2.plot(cumulative, sorted_lc)
    ax2.set_ylim(0, max(lightcurve))
    ax2.set_xlim(0, 1)
    ax2.vlines(np.linspace(0.1, 0.9, 9), 0, roots, lw=0.5)
    ax2.hlines(roots, 0, np.linspace(0.1, 0.9, 9), lw=0.5)

    pl.savefig('figures/spectra_slicing.pdf', bbox_inches='tight')

    print 'Getting flux-sliced intervals...'

    times_metalist = []
    for i, root in enumerate(roots):
        # print i
        temp_times = []
        if i == 0:
            temp_times = times[lightcurve < root]
        elif root == roots[-1]:

            temp_times = times[lightcurve < root]
            temp_lc = lightcurve[lightcurve < root]
            temp_times = temp_times[temp_lc > roots[i - 1]]
            times_metalist.append(temp_times)

            temp_times = times[lightcurve > root]
        else:
            temp_times = times[lightcurve < root]
            temp_lc = lightcurve[lightcurve < root]
            temp_times = temp_times[temp_lc > roots[i - 1]]
        times_metalist.append(temp_times)
    # print sum([len(x) for x in times_metalist])
    # Find the indexes of each timestep
    index_metalist = [np.array([np.where(times == x)[0][0]
                                for x in y]) for y in times_metalist]
    # exit()

    pl.figure()
    ax = pl.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')

    for density in densities:
        print 'Loading spectra for density', density

        infile = h5py.File(
            spectra_dir + 'spectral_data_density_%s.hdf5' % str(density), 'r')
        spectra_dset = infile['spectra']
        energies_dset = infile['energies']

        for i, indices in enumerate(index_metalist):

            mean_spectrum = np.mean(
                spectra_dset[indices[indices < spectra_dset.shape[0]], :], axis=0)
            # print indices[indices<spectra_dset.shape[0]]
            # print spectra_dset.shape[0]
            # print indices
            pl.plot(energies_dset, mean_spectrum)
            np.savetxt('flux_resolved_spectra/density_%s_flux_%s.dat' %
                       (str(density), str(i + 1)), mean_spectrum)
        pl.savefig('figures/spectra_test_%s.pdf' % str(density))
        # exit()
    np.savetxt('flux_resolved_spectra/ebins.dat',energies_dset)


if __name__ == '__main__':
    if not os.path.exists('figures'):
        os.mkdir('figures')
    if not os.path.exists('flux_resolved_spectra'):
        os.mkdir('flux_resolved_spectra')
    main()
