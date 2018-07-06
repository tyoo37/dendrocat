from astropy.io import fits
from astropy import wcs
import radio_beam
import numpy as np
import astropy.units as u
from astrodendro import Dendrogram, pp_catalog
from sourcecatalog import SourceCatalog
import warnings
warnings.filterwarnings('ignore')


class Image:
    """
    An object to store radio image data.
    """

    def __init__(self, hdu, region_id=None, freq_id=None):
        """
        Parameters
        ----------
        hdu : `~astropy.io.fits.hdu.image.PrimaryHDU`
            An astropy FITS HDU object containing the radio image data and 
            header.
        region_id : str
            An identifier specifying what sky object the radio image contains.
        freq_id : str
            An identifier specifying the observation frequency (Ex: 226GHz). 
            If not specified, it will be generated from the FITS image header.
        """

        self.hdu = hdu
        self.hdr = hdu[0].header
        self.data = hdu[0].data.squeeze()
        self.region_id = region_id
        self.freq_id = freq_id
        
        self.wcs = wcs.WCS(self.hdr).celestial
        self.beam = radio_beam.Beam.from_fits_header(self.hdr)
        self.pixel_scale = (np.abs(self.wcs.pixel_scale_matrix.diagonal()
                            .prod())**0.5 * u.deg)
        self.ppbeam = (self.beam.sr/(self.pixel_scale**2)).decompose().value
        self.data = self.data/self.ppbeam
        self.get_fits_info()
    
        # Set default dendrogram values
        self.default_min_value = 1.7*np.std(self.data)
        self.default_min_delta = 1.4*self.default_min_value
        self.default_min_npix = 7
    
    def get_fits_info(self):
        """
        Get information from FITS header.
        
        Supported Telescopes
        ----------
        ALMA
        """
        
        try:
            self.telescope = self.hdr['TELESCOP']
            
            if self.telescope == 'ALMA':
                
                # Get the frequency, either stored in CRVAL3 or CRVAL4
                self.nu = 'UNKNOWN'
                for i in range(len(self.hdr['CTYPE*'])):
                    if self.hdr['CTYPE*'][i] == 'FREQ':
                        self.nu = (self.hdr['CRVAL*'][i] 
                                   * u.Unit(self.hdr['CUNIT*'][i]))
                
                # Create a frequency identifier from nu
                if not self.freq_id:
                    self.freq_id = ('{:.0f}'.format(np.round(self.nu
                                            .to(u.GHz))).replace(' ', ''))
                
                # Get metadata - need to raise exceptions if data is missing
                self.metadata = {
                        'data_unit': u.Unit(self.hdr['BUNIT']),
                        'spatial_scale': self.pixel_scale,
                        'beam_major': self.beam.major,
                        'beam_minor': self.beam.minor,
                        'wavelength': self.nu,
                        'velocity_scale': u.km/u.s,
                        'wcs': self.wcs,
                            }
                            
        except KeyError:
            self.telescope = 'UNKNOWN'
    
    
    def to_dendrogram(self, parameters=None, save=True):
        """
        Calculates a dendrogram for the image.
        
        Parameters
        ----------
        parameters : list or tuple of floats, optional
            (min_value, min_delta, min_npix) to use when calculating the 
            dendrogram. Defaults will be used if not specified.
        save : bool, optional
            If enabled, the resulting dendrogram will be saved as an instance
            attribute. Default is True.
            
        Returns
        ----------
        ~astrodendro.dendrogram.Dendrogram object
            A dendrogram object calculated from the radio image.
        """              
        if parameters:
            min_val = parameters[0]
            min_delt = parameters[1]
            min_npix = parameters[2]
            
        else:
            min_val = self.default_min_value
            min_delt = self.default_min_delta
            min_npix = self.default_min_npix
        
        dend = Dendrogram.compute(self.data, 
                                  min_value=min_val, 
                                  min_delta=min_delt, 
                                  min_npix=min_npix, 
                                  wcs=self.wcs, 
                                  verbose=True)
        if save:
            self.dendrogram = dend
            
        return dend
    
    
    def to_cat(self, dendrogram=None):
        """
        Returns a position-position catalog of leaves in a dendrogram. If no
        dendrogram is specified, the dendrogram saved in the instance 
        attributes will be used if it exists.
        
        Parameters
        ----------
        dendrogram : ~astrodendro.dendrogram.Dendrogram object, optional
            The dendrogram object to extract sources from. If not specified,
            the 'dendrogram' instance attribute will be used if it exists. 
        """
        
        if not dendrogram:
            dendrogram = self.dendrogram
                       
        cat = pp_catalog(dendrogram.leaves, self.metadata)
        cat['_idx'] = range(len(cat))
    
        cat['major_sigma'] = cat['major_sigma']*np.sqrt(8*np.log(2))
        cat['minor_sigma'] = cat['minor_sigma']*np.sqrt(8*np.log(2))
        cat.rename_column('major_sigma', 'major_fwhm')
        cat.rename_column('minor_sigma', 'minor_fwhm')
        cat.rename_column('flux', 'dend_flux_{}'.format(self.freq_id))
        
        catobj = SourceCatalog(catalog=cat, imageobj=self, masked=True)
        catobj.__dict__.update(self.__dict__)
        
        return catobj
        
                                 
