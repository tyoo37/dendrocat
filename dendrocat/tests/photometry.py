import dendrocat
from astropy.io import fits
import astropy.units as u

rs1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs2 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz'))
rs3 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2w_QbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_incrementalselfcal8.image.tt0.pbcor.icrs.fits'))

rs3.nu = 45*u.GHz
rs3.freq_id = '45Ghz'
rs3.set_metadata()

rs1.autoreject()
rs2.autoreject()


rs1.reject([])



mc = dendrocat.MasterCatalog(rs1, rs2, rs3, catalog=dendrocat.match(rs1, rs2).catalog)
mc.photometer(dendrocat.ellipse)
mc.photometer(dendrocat.annulus)


