import regions
import astropy.units as u
import numpy as np
import warnings

from .utils import ucheck

def mask(reg, cutout):
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')
    
class Aperture():

    # This class will replace the previous mechanic for grabbing pixels in
    # certain regions. Now, the user will be able to define custom apertures
    # in addition to the presets. Support for using astropy regions as 
    # apertures will also be added. Thus, it is important that this class 
    # remains similar in function to astropy regions, and that any further
    # development take place in other classes' methods.
    
    def __init__(self, center, major, minor, pa, unit=None):
        '''
        Create an elliptical aperture, defined in either pixel or sky 
        coordinates.
        
        Parameters
        ----------
        center : list or tuple
            x and y (ra and dec) coordinates for the center of the ellipse.
        major : scalar or astropy.units.quantity.Quantity
            Major axis of the ellipse (i.e., longest diameter)
        minor : scalar or astropy.units.quantity.Quantity
            Minor axis of the ellipse (i.e., shortest diameter)
        pa : scalar or astropy.units.quantity.Quantity
            If scalar, assumed to be given in degrees. The position angle of 
            the major axis of the ellipse, measured from the positive x-axis
            toward the positive y-axis. Defined in the range 0 < pa <= 180.
        '''
        
        if unit is None:
            try:
                unit = center[0].unit
                self.unit = unit
            except AttributeError:
                warnings.warn('Unit not specified. Please specify a unit.')
                return
        else:
            self.unit = unit
        
        self.center = ucheck(center, self.unit)
        self.x_cen = ucheck(center[0], self.unit)
        self.y_cen = ucheck(center[1], self.unit)
        self.major = ucheck(major, self.unit)
        self.minor = ucheck(minor, self.unit)
        self.pa = ucheck(pa, u.deg)
        
    
    def place(self, image, wcs=None):
        """
        Place the aperture on an image, either using pixel coordinates or the
        image's WCS.
        
        Parameters
        ----------
        image : array
            The image upon which to place the aperture.
        wcs : astropy.wcs.wcs.WCS object
            The world coordinate system for the image, used to match up 
            apertures given in degrees (presumably RA and DEC).
            
        Returns
        ----------
        numpy.ndarray
            A boolean mask for the aperture with the same dimensions as `image`
        """
        
        
        
        
    def from_region(region):
        """
        Make a dendrocat.aperture.Aperture object from an existing regions
        object.
        
        Parameters
        ----------
        region : astropy regions region
            
        """

        
def ellipse(source, cutout, obj):

    center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])

    pix_major = (2.*source['major_fwhm']*u.deg / obj.pixel_scale).value
    pix_minor = (2.*source['minor_fwhm']*u.deg / obj.pixel_scale).value
    pa = source['position_angle']*u.deg
    
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
