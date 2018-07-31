from astropy.io import fits
from astropy.table import Table
import dendrocat
import astropy.units as u

t = Table(Table.read('/users/bmcclell/nrao/cat/45-93-226GHz_photometered_adjustedRADEC.dat', format='ascii'), masked=True)

rs1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs2 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2_band3_93GHz_adjustedRADEC.fits'))
rs3 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/W51e2w_QbandAarray_adjustedRADEC.fits'))

rs3.nu = 45*u.GHz
rs3.freq_id = '45GHz'
rs3.set_metadata()

mc = dendrocat.MasterCatalog(rs1, rs2, rs3, catalog=t)

ext = Table.read('/users/bmcclell/Data/EVLA_VLA_PointSourcePhotometry.ipac', format='ipac')

ra='gracen'
dec='gdeccen'
freq='Frequency'
flux_sum='aperture_flux'
flux_peak='peak_flux'
err='local_rms_noise'
shape='ellipse'
tolerance=0.1*u.arcsec
skip_rejects=True

me = mc.match_external(ext, ra=ra, dec=dec, freq=freq, flux_sum=flux_sum, 
                       flux_peak=flux_peak, err=err, shape=shape, 
                       tolerance=tolerance)
me.write('/users/bmcclell/nrao/cat/w51e_external_photometered.dat', format='ascii')
