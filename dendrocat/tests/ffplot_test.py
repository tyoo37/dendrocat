import dendrocat
from astropy.io import fits
import matplotlib.pyplot as plt

rs1 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'))
rs2 = dendrocat.RadioSource(fits.open('/lustre/aoc/students/bmcclell/w51/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz'))

rs1.autoreject()
rs2.autoreject()

rs1.reject([0, 12, 59, 89, 98, 99, 106])
rs2.reject([8, 9, 14, 18, 21, 42])

mc = dendrocat.match(rs1, rs2)
dendrocat.ffplot(mc, rs2, rs1)
plt.show()
