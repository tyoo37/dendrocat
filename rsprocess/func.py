import numpy as np
from radio_beam import Beams
import astropy.units as u

def mask(reg, cutout):
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')
    
def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5

def check_units(quantity, unit=u.deg):
    if unit.is_equivalent(quantity):
        return quantity.to(unit)
    else:
        return quantity * u.Unit(unit)

def commonbeam(major1, minor1, pa1, major2, minor2, pa2):
    """
    Create a smallest bounding ellipse around two other ellipses. 
    Give ellipse dimensions as astropy units quantities.
    """
    major1 = check_units(major1, unit=u.deg)
    minor1 = check_units(minor1, unit=u.deg)
    pa1 = check_units(pa1, unit=u.deg)
    major2 = check_units(major2, unit=u.deg)
    minor2 = check_units(minor2, unit=u.deg)
    pa2 = check_units(pa2, unit=u.deg)
    
    somebeams = Beams([major1.to(u.arcsec), major2.to(u.arcsec)]*u.arcsec, 
                      [minor1.to(u.arcsec), minor2.to(u.arcsec)]*u.arcsec, 
                      [pa1, pa2]*u.deg)
                      
    common = somebeams.common_beam()
    new_major = common._major
    new_minor = common._minor
    new_pa = common._pa
    
    return new_major.to(u.deg), new_minor.to(u.deg), new_pa

def save_ds9_regions(catalog, outfile, hide_rejects=True):
    
        if hide_rejects:
            catalog = catalog[np.where(catalog['rejected'] == 0)]
        else:
            catalog = catalog
            
        with open(outfile, 'w') as fh:
            fh.write("icrs\n")
            for row in catalog:
                fh.write("ellipse({x_cen}, {y_cen}, {major_fwhm}, " \
                         "{minor_fwhm}, {position_angle}) # text={{{_idx}}}\n"
                         .format(**dict(zip(row.colnames, row))))
