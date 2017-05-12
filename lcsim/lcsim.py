"""
Simulate the noise and flux smearing for a light curve with the perfect
simulared flux based on the observation logs obtained from SIMLIB files
"""
import numpy as np


class LCSim():
    """
    Simulate light curve noise and smearing based on SIMLIBReader output

    Properties
    ----------
    seed : int, optional
        Seed for random number generators used thoughout the class
    """
    def __init__(self, seed=0):
        self.seed = seed

    def compute_statistical_error(self):
        """
        Compute the statistical uncertaintly from a set of observing logs

        Returns
        -------
        snr : np.ndarray
            signal-to-noise ratio for the a given observation
        """
        zptdif_pe = (self.simlib['zpt'] +
                     2.5*np.log10(self.simlib['gain']) - 27.5)
        self.pe_over_fluxcal = 10**(0.4*zptdif_pe)

        self.flux_adu = 10.0**(0.4 * (self.simlib['zps'] - self.genmag))
        flux_pe = self.flux_adu * self.simlib['gain']

        area_bg = 4.0 * np.pi * self.simlib['psf1']**2

        sqskyerr_pe = area_bg*(self.simlib['skysigs'] * self.simlib['gain'])**2
        sqccderr_pe = area_bg * self.simlib['noise']**2

        template_sqskyerr_pe = area_bg * (self.simlib['skysigs'] *
                                          self.simlib['gain'])**2
        zfac = 10.0**(0.8*(self.simlib['zps'] - self.simlib['zpt']))
        template_sqskyerr_pe *= zfac
        self.template_pe_err = np.sqrt(template_sqskyerr_pe)

        sqsum = flux_pe + sqskyerr_pe + sqccderr_pe
        self.flux_pe_err = np.sqrt(sqsum)

        return flux_pe / self.flux_pe_err

    def scale_fluxerr_model(self):
        SPol = {
            'g': [4.53181, -4.33347,  1.962600, -0.3855380, 0.0280264],
            'r': [3.70768, -3.44778,  1.612420, -0.3268320, 0.0245651],
            'i': [2.31471, -1.67548,  0.773386, -0.1527000, 0.0112755],
            'z': [1.60241, -0.829311, 0.417534, -0.0892309, 0.00715445]
        }

        DPol = {
            'g': [1.08856, 0.02900760, 0, 0, 0],
            'r': [1.06293, 0.02090900, 0, 0, 0],
            'i': [1.03998, 0.04424530, 0, 0, 0],
            'z': [1.17307, 0.00649704, 0, 0, 0]
        }

        field = self.simlib['field'].values[0]
        band = self.simlib['band'].values[0]
        if field == 'C3' or field == 'X3':
            Pol = DPol
        else:
            Pol = SPol

        scale = 0.0
        for i in range(5):
            scale += Pol[band][i] * self.simlib['psf1']**i
        return scale

    def smear_generated_mag(self):
        """
        Apply random smearing to the simulated photometry based on the
        observing logs
        """
        flux_adu_errS = self.flux_pe_err / self.simlib['gain']
        template_adu_err = self.template_pe_err / self.simlib['gain']

        relerr = 10.0**(0.4*self.simlib['zps']) - 1.0
        err1 = flux_adu_errS
        err2 = self.flux_adu * relerr
        flux_adu_errSZ = np.sqrt(err1**2 + err2**2)

        gaussian_search = np.random.normal(size=self.flux_adu.size)
        gaussian_template = np.random.normal(size=self.flux_adu.size)

        flux_adu_err_real = flux_adu_errSZ
        flux_obs_adu = (self.flux_adu +
                        flux_adu_err_real * gaussian_search +
                        template_adu_err * gaussian_template)

        mask = (flux_obs_adu >= 0).astype(int)
        sqerr_ran = (flux_obs_adu * mask - self.flux_adu) / self.simlib['gain']
        sqsum = flux_adu_errSZ**2 + template_adu_err**2 + sqerr_ran
        flux_adu_errSZT = np.sqrt(sqsum)

        sqsum = flux_adu_errS**2 + template_adu_err**2 + sqerr_ran
        errstat = np.sqrt(sqsum)

        scale_flux_err = self.scale_fluxerr_model()
        flux_adu_errSZT *= scale_flux_err
        errstat *= scale_flux_err

        zp_scale = 10**(-0.4*(self.simlib['zps'] - 31.4))
        fluxcal = flux_obs_adu * zp_scale
        fluxcal_err = errstat * zp_scale

        relerr = 10**(0.4*self.simlib['sigzps']) - 1.0
        tmperr = fluxcal * relerr
        fluxcal_err = np.sqrt(fluxcal_err**2 + tmperr**2)

        return fluxcal, fluxcal_err

    def simulate(self, generated_mag, simlib):
        """

        """
        self.genmag = generated_mag
        self.simlib = simlib

        self.compute_statistical_error()
