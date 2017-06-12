"""
Simulate the noise and flux smearing for a light curve with the perfect
simulared flux based on the observation logs obtained from SIMLIB files
"""
import numpy as np


class LCSim():
    """
    Simulate light curve noise and smearing based on SIMLIBReader output

    Parameters
    ----------
    seed : int, optional
        Seed for random number generators used thoughout the class
    """
    def __init__(self, seed=None):
        if seed is not None:
            np.random.seed(seed)

    def psf_area(self):
        """
        Calculate the area of the effective area of the PSF
        Properties
        ----------
        psf1 : float
            First axis of the PSF

        psf2 : float
            Second axis of the PSF

        psfratio : float
            Ratio between the two psf axis

        Returns
        -------
        area_bg : float
            Area of the background
        """
        if ((self.simlib['psfratio'][0] < 1e-5) or
           (self.simlib['psf2'][0] < 1e-4)):
            return 4.0 * np.pi * self.simlib['psf1']**2

        else:
            tmp = (self.simlib['psfratio'] * (self.simlib['psf2']**2) /
                   (self.simlib['psf1']**2))
            a1 = 1.0 / (1.0 + tmp)
            a2 = 1.0 - a1
            tmp = ((a1 * self.simlib['psf2'] / self.simlib['psf1']) +
                   (a2 * self.simlib['psf1'] / self.simlib['psf2']))
            return 4.0 * np.pi * (((self.simlib['psf2']**2) *
                                  (self.simlib['psf1']**2)) / (1.0 + tmp**2))

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

        flux_pe = self.flux_adu * self.simlib['gain']

        area_bg = self.psf_area()

        sqskyerr_pe = area_bg*(self.simlib['skysigs'] * self.simlib['gain'])**2
        sqccderr_pe = area_bg * self.simlib['noise']**2

        self.template_pe_err = np.zeros_like(self.simlib['skysigt'])
        if self.simlib.survey == 'DES':
            template_sqskyerr_pe = area_bg * (self.simlib['skysigt'] *
                                              self.simlib['gain'])**2
            zfac = 10.0**(0.8*(self.simlib['zps'] - self.simlib['zpt']))
            template_sqskyerr_pe *= zfac
            self.template_pe_err = np.sqrt(template_sqskyerr_pe)

        sqsum = np.abs(flux_pe) + sqskyerr_pe + sqccderr_pe
        self.flux_pe_err = np.sqrt(sqsum)

        return flux_pe, self.flux_pe_err

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
        if field == 'C3' or field == 'X3':
            Pol = DPol
        else:
            Pol = SPol

        poly_index = self.simlib['flt'].map(Pol)
        poly_index = np.array([np.array(poly_row) for poly_row in poly_index])
        scale = np.zeros(poly_index.shape[0])
        for i in range(5):
            scale += poly_index[:, i] * self.simlib['psf1']**i

        return scale

    def smear_generated_mag(self):
        """
        Apply random smearing to the simulated photometry based on the
        observing logs
        """
        flux_adu_errS = self.flux_pe_err / self.simlib['gain']
        template_adu_err = self.template_pe_err / self.simlib['gain']

        relerr = 10.0**(0.4*self.simlib['sigzps']) - 1.0
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
        sqerr_ran = np.abs((flux_obs_adu * mask - self.flux_adu) /
                           self.simlib['gain'])
        sqsum = flux_adu_errSZ**2 + template_adu_err**2 + sqerr_ran
        flux_adu_errSZT = np.sqrt(sqsum)

        sqsum = flux_adu_errS**2 + template_adu_err**2 + sqerr_ran
        errstat = np.sqrt(sqsum)

        if self.simlib.survey == 'DES':
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

    def simulate(self, generated_flux, simlib, unit='nJy'):
        """
        Simulate observed flux and flux error for a given simulated magnitude
        or flux

        Parameters
        ----------
        generated_flux : ndarray
            Array of generated fluxes or magnitude this should not have any
            scatter other than model uncertainties, must have the same size as
            `simlib`. generated_flux must be in the same order as the MJDs
            in simlib.

        simlib : pandas.DataFrame
            DataFrame of observing logs created by SIMLIBReader.get_obs() file.

        unit : str, optional
            Unit of the input generated_flux. Currently supports:
            `mag` (AB magnitudes), `nJy` (nano Jansky) or kessler (10^-11 Jy)

        Returns
        -------
        fluxcal : ndarray
            Flux array with image noise smear applied.

        fluxcal_err : ndarray
            Flux error array.
        """
        self.simlib = simlib

        if unit == 'kessler':
            self.flux_adu = generated_flux

        elif unit == 'nJy':
            zp_scale = 10**(-0.4*(31.4 - self.simlib['zps']))
            self.flux_adu = generated_flux * zp_scale

        elif unit == 'mag' or unit == 'ab':
            self.flux_adu = 10.0**(0.4 * (self.simlib['zps'] - generated_flux))

        self.compute_statistical_error()
        return self.smear_generated_mag()
