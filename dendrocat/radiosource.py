from astropy.io import fits
from astropy import wcs
import radio_beam
import numpy as np
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D
from astropy.table import Column, Table
from astrodendro import Dendrogram, pp_catalog
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import regions
import warnings
warnings.filterwarnings('ignore')

if __package__ == '':
    __package__ = 'dendrocat'
from .utils import annulus, ellipse, rms

class RadioSource:
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
        self.default_min_value = 1.7*np.nanstd(self.data)
        self.default_min_delta = 1.4*self.default_min_value
        self.default_min_npix = 7
        
        # Set other default parameters
        self.default_threshold = 6.
        self.annulus_width = 1e-5 * u.deg
        self.annulus_padding = 1e-5 * u.deg
    
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
        Documentation needed
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
    
    
    def to_catalog(self, dendrogram=None, catalog=None):
        """
        Returns an astropy.table.Table object containing a position-position 
        catalog of leaves in a dendrogram.
        
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
        
        cat.add_column(Column(np.ones(len(cat)), dtype=int), 
                       name='detected_'+self.freq_id)
                       
        self.catalog = Table(cat, masked=True)
        
        return Table(cat, masked=True)


    def _make_cutouts(self, sidelength, catalog=None, save=True):
        """
        Make a cutout_data of cutout regions around all source centers in the 
        catalog.
        
        Parameters
        ----------
        sidelength : float
            Side length of the square (in degrees) to cut out of the image for 
            each source.
        save : bool, optional
            If enabled, the cutouts and cutout data will both be saved as 
            instance attributes. Default is True.
            
        Returns
        ----------
        List of astropy.nddata.utils.Cutout2D objects, list of cutout data
            
        """
        
        # remove this, replace names at some point (lazy coding)
        beam = self.beam
        wcs = self.wcs
        pixel_scale = self.pixel_scale
        data = self.data
        
        if catalog is None:   
            try:
                catalog = self.catalog
            except AttributeError:
                catalog = self.to_catalog()
        
        cutouts = []
        cutout_data = []
        
        for i in range(len(catalog)):
            x_cen = catalog['x_cen'][i] * u.deg
            y_cen = catalog['y_cen'][i] * u.deg
            
            position = coordinates.SkyCoord(x_cen, y_cen, frame='icrs',
                                            unit=(u.deg, u.deg))
            pixel_position = np.array(position.to_pixel(wcs))
            
            cutout = Cutout2D(data, position, sidelength, wcs, mode='partial')
            cutouts.append(cutout)
            cutout_data.append(cutout.data)
            
        cutouts = np.array(cutouts)
        cutout_data = np.array(cutout_data)
        
        if save:
            self._cutouts = cutouts
            self._cutout_data = cutout_data
        # NOTE: If 'sort' is called, the catalog's attributes also need to be
        # sorted accordingly. Might be tricky.
        
        return cutouts, cutout_data
    
    
    def get_pixels(self, aperture, cutouts=None, 
                   cutout_data=None, save=True, **kwargs):
        """
        Return a list of pixel arrays, each of which contains the pixels in
        an annulus of constant width and variable radius depending on the 
        major fwhm of the source.
        
        Parameters
        ----------
        Documentation needed!
        save : bool, optional
            If enabled, the pixel arrays and masks will both be saved as 
            instance attributes. Default is True.
            
        Returns
        ----------
        List of pixel arrays
        """
        try:
            catalog = self.catalog
        except AttributeError:
            catalog = self.to_catalog()
        
        if cutouts is None or cutout_data is None:
            try:
                cutouts = self._cutouts
                cutout_data = self._cutout_data
            except AttributeError:
                size = 2.2*(np.max(catalog['major_fwhm'])*u.deg 
                            + self.annulus_padding 
                            + self.annulus_width)
                cutouts, cutout_data = self._make_cutouts(size)
        
        pix_arrays = []
        masks = []
        
        for i in range(len(cutouts)):
            this_mask = aperture(self.catalog[i], 
                                 cutouts[i], 
                                 self,
                                 **kwargs)
                                    
            pix_arrays.append(cutouts[i].data[this_mask.astype('bool')])
            masks.append(this_mask)
        
        if save:
            self.__dict__['pixels_{}'.format(aperture.__name__)] = pix_arrays
            self.__dict__['mask_{}'.format(aperture.__name__)] = masks
        
        return pix_arrays
        
    
    def get_snr(self, pixels_in_source=None, pixels_in_background=None, 
                peak=True, save=True):
        
        """
        Return the SNR of all sources in the catalog.
        
        Parameters
        ----------
        peak : bool, optional
            Use peak flux of source pixels as 'signal'. Default is True.
        save : bool, optional
            If enabled, the snr will be saved as a column in the source catalog
            and as an instance attribute. Default is True.
        """
        
        if not pixels_in_background:
            try:
                pixels_in_background = self.pixels_annulus
            except:
                pixels_in_background = self.get_pixels(annulus)
                                                               
        if not pixels_in_source:
            try:
                pixels_in_source = self.pixels_ellipse
            except AttributeError:
                pixels_in_source = self.get_pixels(ellipse)
        
        snr_vals = []
        for i in range(len(self.catalog)):
            snr = (np.max(pixels_in_source[i])
                  / rms(pixels_in_background[i]))
            snr_vals.append(snr)
            
        if save:
            self.snr = np.array(snr_vals)
            self.catalog.add_column(Column(snr_vals), name='snr_'+self.freq_id)
            
        return np.array(snr_vals)
    
    
    def plot_grid(self, cutout_data=None, masks=None, snr_vals=None,
                  skip=False):
        """
        Plot a grid of sources with aperture mask overlays. Rejected sources
        are shown in gray.
        
        Parameters
        ----------
        Documentation needed
        """
        
        if not snr_vals:
            try:
                snr_vals = self.snr
            except AttributeError:
                snr_vals = self.get_snr()
                
        if not cutout_data:
            cutout_data = self._cutout_data
        
        if not masks:
            masks = np.array(list(zip(self.mask_ellipse, self.mask_annulus)))
            
        names = np.array(self.catalog['_idx'])
        rejected = np.array(self.catalog['rejected'])
        
        if skip:
            accepted_indices = np.where(self.catalog['rejected'] == 0)[0]
            snr_vals = snr_vals[accepted_indices]
            cutout_data = cutout_data[accepted_indices]
            masks = masks[accepted_indices]
            names = names[accepted_indices]
            rejected = rejected[accepted_indices]
        
        n_images = len(cutout_data)
        xplots = int(np.around(np.sqrt(n_images)))
        yplots = xplots + 1
        gs1 = gs.GridSpec(yplots, xplots, wspace=0.0, hspace=0.0)
        plt.figure(figsize=(9.5, 10))
        
        for i in range(n_images):
            image = cutout_data[i]
            plt.subplot(gs1[i])
            
            if rejected[i] == 1:
                plt.imshow(image, origin='lower', cmap='gray')
            else:
                plt.imshow(image, origin='lower')

            for j in range(len(masks[i])):
                plt.imshow(masks[i][j], origin='lower', cmap='gray', 
                           alpha=0.15)
                
            plt.text(0, 0, '{}  SN {:.1f}'.format(names[i], snr_vals[i]), 
                     fontsize=7, color='w')
            plt.xticks([])
            plt.yticks([])
        
        plt.tight_layout()
        plt.show()
    
    
    def autoreject(self, threshold=6.):
        """
        Reject noisy detections.
        
        Parameters
        ----------
        threshold : float, optional
            The signal-to-noise threshold below which sources are rejected
        """
            
        if not threshold:
            threshold = self.default_threshold
      
        try:
            snrs = self.snr
        except:
            snrs = self.get_snr()
            
        try:
            self.catalog['rejected'] = np.zeros(len(self.catalog), dtype=int)
        except KeyError:
            self.catalog.add_column(Column(np.zeros(len(self.catalog))), 
                                    name='rejected')
    
        for i in range(len(self.catalog)):
            if snrs[i] <= threshold:
                self.catalog['rejected'][i] = 1
                    

    def reject(self, rejected_list):
        
        for idx in rejected_list:
            self.catalog[np.where(self.catalog['_idx'] == idx)]['rejected'] = 1
            
            
    def accept(self, accepted_list):
        
        for idx in accepted_list:
            self.catalog[np.where(self.catalog['_idx'] == idx)]['rejected'] = 0

        
