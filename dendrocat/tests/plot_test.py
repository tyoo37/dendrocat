import dendrocat
from astropy.io import fits
import astropy.units as u

rs1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs2 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz'))
rs3 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2w_QbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_incrementalselfcal8.image.tt0.pbcor.icrs.fits'))

rs3.nu = 45*u.GHz
rs3.freq_id = '45Ghz'

def set_metadata(obj):
    obj.metadata = {
        'data_unit': u.Unit(obj.header['BUNIT']),
        'spatial_scale': obj.pixel_scale,
        'beam_major': obj.beam.major,
        'beam_minor': obj.beam.minor,
        'wavelength': obj.nu,
        'velocity_scale': u.km/u.s,
        'wcs': obj.wcs,
            }

set_metadata(rs3)

rs1.autoreject()
rs2.autoreject()

mc = dendrocat.match(rs1, rs2)
mc.add_objects(rs3)
mc.photometer(dendrocat.annulus)
mc.photometer(dendrocat.ellipse)
mc.catalog.write('45-93-226GHz_photometry_ellipse_annulus.dat', format='ascii')

dendrocat.utils.save_regions(mc.catalog, '226GHz_93GHz_matched_regions.reg')
#rs3.plot_grid(catalog=mc.catalog, outfile='plotgrid_45GHzdata_226and93GHzapertures.png')


