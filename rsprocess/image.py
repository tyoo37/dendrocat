from astropy.io import fits
from astropy import wcs
import radio_beam
import numpy as np
import astropy.units as u
from astropy.table import Column
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
        freq_id : str, optional
            An identifier specifying the observation frequency (Ex: 226GHz). 
            If not specified, it will be generated from the FITS image header.
        """

        self.hdu = hdu
        self.header = hdu[0].header
        self.data = hdu[0].data.squeeze()
        self.region_id = region_id
        self.freq_id = freq_id
        
        self.wcs = wcs.WCS(self.header).celestial
        self.beam = radio_beam.Beam.from_fits_header(self.header)
        self.pixel_scale = (np.abs(self.wcs.pixel_scale_matrix.diagonal()
                            .prod())**0.5 * u.deg)
        self.ppbeam = (self.beam.sr/(self.pixel_scale**2)).decompose().value
        self.data = self.data/self.ppbeam
        self._get_fits_info()
    
        # Set default dendrogram values
        self.default_min_value = 1.7*np.std(self.data)
        self.default_min_delta = 1.4*self.default_min_value
        self.default_min_npix = 7
        
        # Set other default parameters
        self.default_threshold = 6.
        self.default_annulus_width = 1e-5 * u.deg
        self.default_annulus_padding = 1e-5 * u.deg
    
    def _get_fits_info(self):
        """
        Get information from FITS header.
        
        Supported Telescopes
        ----------
        ALMA
        """
        
        try:
            self.telescope = self.header['TELESCOP']
            
            if self.telescope == 'ALMA':
                
                # Get the frequency, either stored in CRVAL3 or CRVAL4
                self.nu = 'UNKNOWN'
                for i in range(len(self.header['CTYPE*'])):
                    if self.header['CTYPE*'][i] == 'FREQ':
                        self.nu = (self.header['CRVAL*'][i] 
                                   * u.Unit(self.header['CUNIT*'][i]))
                
                # Create a frequency identifier from nu
                if not self.freq_id:
                    self.freq_id = ('{:.0f}'.format(np.round(self.nu
                                            .to(u.GHz))).replace(' ', ''))
                
                # Get metadata - need to raise exceptions if data is missing
                self.metadata = {
                        'data_unit': u.Unit(self.header['BUNIT']),
                        'spatial_scale': self.pixel_scale,
                        'beam_major': self.beam.major,
                        'beam_minor': self.beam.minor,
                        'wavelength': self.nu,
                        'velocity_scale': u.km/u.s,
                        'wcs': self.wcs,
                            }
                            
        except KeyError:
            self.telescope = 'UNKNOWN'
    
    
    def to_dendrogram(self, min_value=None, min_delta=None, min_npix=None, 
                      save=True):
        """
        Calculates a dendrogram for the image.
        
        Parameters
        ----------
        min_value : float, optional
        
        min_delta : float, optional
        
        min_npix : float, optional
            
        save : bool, optional
            If enabled, the resulting dendrogram will be saved as an instance
            attribute. Default is True.
            
        Returns
        ----------
        ~astrodendro.dendrogram.Dendrogram object
            A dendrogram object calculated from the radio image.
        """              
      
        if not min_value:
            min_value = self.default_min_value
        
        if not min_delta:
            min_delta = self.default_min_delta
            
        if not min_npix:
            min_npix = self.default_min_npix
        
        dend = Dendrogram.compute(self.data, 
                                  min_value=min_value, 
                                  min_delta=min_delta, 
                                  min_npix=min_npix, 
                                  wcs=self.wcs, 
                                  verbose=True)
        if save:
            self.dendrogram = dend
            
        return dend
    
    
    def to_sourcecatalog(self, dendrogram=None, catalog=None):
        """
        Returns a rsprocess.SourceCatalog object containing a position-position 
        catalog of leaves in a dendrogram, with image information stored as
        instance attributes.
        
        Parameters
        ----------
        dendrogram : ~astrodendro.dendrogram.Dendrogram object, optional
            The dendrogram object to extract sources from.
        catalog : ~astropy.table.Table object, optional
            A position-position catalog, as output by astrodendro.pp_catalog.
        """
        
        if not dendrogram:
            try:
                dendrogram = self.dendrogram
            except AttributeError:
                dendrogram = self.to_dendrogram()
                
        if not catalog:
            try:
                cat = self.catalog
            except AttributeError:
                cat = pp_catalog(dendrogram.leaves, self.metadata)

        cat['_idx'] = range(len(cat))
    
        cat['major_sigma'] = cat['major_sigma']*np.sqrt(8*np.log(2))
        cat['minor_sigma'] = cat['minor_sigma']*np.sqrt(8*np.log(2))
        cat.rename_column('major_sigma', 'major_fwhm')
        cat.rename_column('minor_sigma', 'minor_fwhm')
        cat.rename_column('flux', 'dend_flux_{}'.format(self.freq_id))
        cat.add_column(Column(np.zeros(len(cat)), dtype=int), name='rejected')
        
        catobj = SourceCatalog(catalog=cat, imageobj=self, masked=True)
        catobj.__dict__.update(self.__dict__)
        
        return catobj
        
                                 
