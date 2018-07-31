import regions
import astropy.units as u
import numpy as np


def mask(reg, cutout):
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')
    
    
def ellipse(source, cutout, obj):

    center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

    pix_major = (2.*source['major_fwhm']*u.deg / obj.pixel_scale).value
    pix_minor = (2.*source['minor_fwhm']*u.deg / obj.pixel_scale).value
    pa = source['position_angle']*u.deg
    
    radius = source['major_fwhm'] * u.deg
    reg = regions.EllipsePixelRegion(center, pix_major, pix_minor, angle=pa)
                     
    ellipse_mask = mask(reg, cutout)
    return ellipse_mask


def annulus(source, cutout, obj):
    
    center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
    inner_r = source['major_fwhm']*u.deg + obj.annulus_padding
    outer_r = inner_r + obj.annulus_width
    
    innerann_reg = regions.CirclePixelRegion(center, (inner_r/obj.pixel_scale).value)
    outerann_reg = regions.CirclePixelRegion(center, (outer_r/obj.pixel_scale).value)
    
    annulus_mask = (mask(outerann_reg, cutout) - mask(innerann_reg, cutout))
    return annulus_mask


def circle(source, cutout, obj, radius=None):
        
    if radius is None:
        radius = source['major_fwhm'] * u.deg
    else:
        radius = check_units(radius, unit=u.deg)
    
    center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
    reg = regions.CirclePixelRegion(center, radius/obj.pixel_scale)
    
    circle_mask = mask(reg, cutout)
    return circle_mask
