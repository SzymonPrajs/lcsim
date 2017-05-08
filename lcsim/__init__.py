"""
lcsim: Python module for working with SNANA's SIMLIB and HOSTLIB files

This module is designed to simplify the process of simulating light curves in
astronomical surveys (e.g. DES, LSST etc) by separating the popularly used
*.SIMLIB files produced by SNANA (Kessler'09) from the rest of the package.

Example workflow using lcsim would be as follows:
-> Load *.SIMLIB file
-> Generate fake observing logs e.g. lists of MJDs per field and band
-> *Externally* compute perfect light curve fluxes or magnitudes
-> Pass fluxes corresponding to the observing logs back to lcsim
-> Simulate uncertainties and flux smearing for the light curves
"""
