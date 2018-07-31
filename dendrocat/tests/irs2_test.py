from astropy.io import fits
import dendrocat
from astropy.table import Table, Column, MaskedColumn, vstack, hstack
import astropy.units as u
n1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51n_cont_briggsSC_tclean.image.fits'))
n2 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51North_QbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_selfcal16_final.image.tt0.pbcor.fits'))
n2.nu = 45*u.GHz
n2.freq_id = '45.0GHz'
n2.set_metadata()
n1.autoreject()
n2.min_value = 1e-4
n2.min_delta = 1.4*1e-4
n2.autoreject()
n1.reject([226003, 226004, 226026])
n1.accept([226047, 226075, 226051])
mc = dendrocat.match(n1, n2)
mc.photometer(dendrocat.aperture.ellipse)
mc.photometer(dendrocat.aperture.annulus)

