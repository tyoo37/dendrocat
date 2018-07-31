from astropy.table import Table
from astropy.io import fits
import dendrocat
import astropy.units as u 
import numpy as np
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from collections import OrderedDict



irs2 = Table.read('/users/bmcclell/nrao/cat/w51IRS2_photometered.dat', format='ascii')

n1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51n_cont_briggsSC_tclean.image.fits'))
n2 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51North_QbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal16_final.image.tt0.pbcor.fits'))

n2.nu = 45.0*u.GHz
n2.freq_id = '45.0GHz'
n2.set_metadata()

mc = dendrocat.MasterCatalog(n1, n2, catalog=irs2)
accepted = mc.catalog[mc.catalog['rejected']==0]

for i in range(len(accepted)):
    nus, fluxes, errs = mc.plotsedgrid(accepted[i], alphas=[1, 2, 3], path='/users/bmcclell/nrao/documentation/SEDS_IRS2/')
    print('   nus', nus)
    print('fluxes', fluxes)
    print('errors', errs)
