import regions
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import warnings

from .utils import ucheck

class NoUnitError(Exception):
    pass
    
class NoWCSError(Exception):
    pass

class Aperture():
    
    def __init__(self, center, major, minor, pa, unit=None, frame='icrs'):
        '''
        Create an elliptical aperture, defined in either pixel (x, y) or sky 
        (ra, dec) coordinates.
        
        Parameters
        ----------
        center : list or tuple, as scalar or astropy.units.quantity.Quantity
            x and y (ra and dec) coordinates for the center of the ellipse.
        major : scalar or astropy.units.quantity.Quantity
            Major axis of the ellipse (i.e., longest diameter)
        minor : scalar or astropy.units.quantity.Quantity
            Minor axis of the ellipse (i.e., shortest diameter)
        pa : scalar or astropy.units.quantity.Quantity
            If scalar, assumed to be given in degrees. The position angle of 
            the major axis of the ellipse, measured from the positive x-axis
            toward the positive y-axis. Defined in the range 0 < pa <= 180.
        unit : astropy.unit.Unit or str
            The unit in which all other arguments are specified. Usually u.pix
            or u.deg.
        frame : str, optional
            The coordinate frame in which (ra, dec) coordinates are specified.
            Default is 'icrs'.
        '''
        
        if unit is None:
            try:
                unit = major.unit
                self.unit = unit
            except AttributeError:
                raise NoUnitError('No unit was specified. Use the keyword \
                                  argument `unit=`, and give an astropy units\
                                  object.') 
        else:
            if type(unit) is str:
                unit = u.Unit(unit)
            self.unit = unit
        
        self.center = ucheck(center, self.unit)
        
        self.major = ucheck(major, self.unit)
        self.minor = ucheck(minor, self.unit)
        self.pa = ucheck(pa, u.deg)
        self.frame = frame
    
    def _refresh_xycen(self):
        if type(self.center) == SkyCoord:
            self.x_cen = ucheck(self.center._sky_coord_frame._data._lon.value, self.unit)
            self.y_cen = ucheck(self.center._sky_coord_frame._data._lat.value, self.unit)
        
        elif type(self.center) == regions.PixCoord:
            self.x_cen = ucheck(self.center.x, self.unit)
            self.y_cen = ucheck(self.center.y, self.unit)
        
        else:
            self.x_cen = ucheck(self.center[0], self.unit)
            self.y_cen = ucheck(self.center[1], self.unit)
        return self.x_cen, self.y_cen
        
    def place(self, image, wcs=None):
        """
        Place the aperture on an image.
        
        Parameters
        ----------
        image : array
            The image upon which to place the aperture.
        wcs : astropy.wcs.wcs.WCS object, optional
            The world coordinate system for the image, used for coordinate 
            transformations.
        
        Returns
        ----------
        numpy.ndarray
            A boolean mask for the aperture with the same dimensions as `image`
        """
        self._refresh_xycen()
        if self.unit.is_equivalent(u.deg) and wcs is not None:
            pixel_scale = (np.abs(wcs.pixel_scale_matrix.diagonal()
                                      .prod())**0.5 * u.deg/u.pix)
            center = np.array(SkyCoord(self.x_cen.to(u.deg), 
                                       self.y_cen.to(u.deg), 
                                       frame=self.frame, 
                                       unit=(u.deg, u.deg)).to_pixel(wcs))
            center = regions.PixCoord(center[0], center[1])
            major = (self.major.to(u.deg)/pixel_scale).value
            minor = (self.minor.to(u.deg)/pixel_scale).value
            reg = regions.EllipsePixelRegion(center, major, minor, angle=self.pa)
             
        elif self.unit.is_equivalent(u.pix):
            center = regions.PixCoord(self.x_cen.value, self.y_cen.value)
            reg = regions.EllipsePixelRegion(center, self.major.value, 
                                             self.minor.value, angle=self.pa)
        else:
            raise NoWCSError('No WCS given.')    
        
        m, n = image.shape
        mask = reg.to_mask(mode='center')
        return np.array(mask.to_image((m, n)), dtype='bool')
        
        
    def from_region(region):
        """
        Make a dendrocat.aperture.Aperture object from an existing regions
        object.
        
        Parameters
        ----------
        region : astropy regions region
            
        """
 
 
class Ellipse(Aperture):
     
    def __init__(self, center, major, minor, pa, unit=None, frame='icrs'):
       Aperture.__init__(self, center, major, minor, pa, unit=unit)
       self.__name__ = 'ellipse'
    def place(self, image, wcs=None):
       return Aperture.place(self, image, wcs=wcs)


class Annulus(Aperture):

    def __init__(self, center, inner, outer, unit=None, frame='icrs'):
        if unit is None:
            try:
                unit = inner.unit
                self.unit = unit
            except AttributeError:
                raise NoUnitError('No unit was specified. Use the keyword' 
                                  'argument `unit=`, and give an astropy'  
                                  'units object.')
        else:
            self.unit = unit
        self.aperture_inner = Aperture(center, inner, inner, 0, unit=unit)
        self.aperture_outer = Aperture(center, outer, outer, 0, unit=unit)
        self.center = ucheck(center, self.unit)
        self.x_cen = ucheck(center[0], self.unit)
        self.y_cen = ucheck(center[1], self.unit)
        self.inner = ucheck(inner, self.unit)
        self.outer = ucheck(outer, self.unit)
        self.__name__ = 'annulus'
        
    def place(self, image, wcs=None):
        return (self.aperture_outer.place(image, wcs=wcs)
                ^ self.aperture_inner.place(image, wcs=wcs))


class Circle(Aperture):
    
    def __init__(self, center, radius, unit=None, frame='icrs'):
        Aperture.__init__(self, center, radius, radius, 0, unit=unit)
        self.radius = ucheck(radius, self.unit)
        self.__name__ = 'circle'
        
    def place(self, image, wcs=None):
        return Aperture.place(self, image, wcs=wcs)
