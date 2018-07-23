from astropy.io import fits
from astropy import wcs
import radio_beam
import numpy as np
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D, NoOverlapError
from astropy.table import Column, Table
from astrodendro import Dendrogram, pp_catalog
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import regions
import warnings
warnings.filterwarnings('ignore')

if __package__ == '':
    __package__ = 'dendrocat'
from .aperture import ellipse, annulus
from .utils import rms

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
        self._get_fits_info()
    
        # Set default dendrogram values
        self.min_value = 1.7*np.nanstd(self.data)
        self.min_delta = 1.4*self.min_value
        self.min_npix = 7
        
        # Set other default parameters
        self.threshold = 6.
        self.annulus_width = 12 * self.pixel_scale
        self.annulus_padding = 12 * self.pixel_scale
    
        self.properties = {
            'min_value':self.min_value,
            'min_delta':self.min_delta,
            'min_npix':self.min_npix,
            'annulus_width':self.annulus_width,
            'annulus_padding':self.annulus_padding
                            }
    
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
                self.set_metadata()
                
            else:
                print('FITS info collection not currently supported for ' \
                      '{}. Please manually set the following instance' \
                      ' attributes:'.format(self.telescope))
                print(' nu\n', 'freq_id\n', 'metadata\n')
                
        except KeyError:
            self.telescope = 'UNKNOWN'
            print('Telescope not identified. Please manually set the ' \
                  'following instance attributes:')
            print(' telescope\n', 'nu\n', 'freq_id\n', 'metadata\n')
    
    
    def set_metadata(self):
        self.metadata = {
            'data_unit': u.Unit(self.header['BUNIT']),
            'spatial_scale': self.pixel_scale,
            'beam_major': self.beam.major,
            'beam_minor': self.beam.minor,
            'wavelength': self.nu,
            'velocity_scale': u.km/u.s,
            'wcs': self.wcs,
                        }
    
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
            min_value = self.min_value
        
        if not min_delta:
            min_delta = self.min_delta
            
        if not min_npix:
            min_npix = self.min_npix
        
        dend = Dendrogram.compute(self.data, 
                                  min_value=min_value, 
                                  min_delta=min_delta, 
                                  min_npix=min_npix, 
                                  wcs=self.wcs, 
                                  verbose=True)
        if save:
            self.dendrogram = dend
            
        return dend
    
    
    def to_catalog(self, dendrogram=None):
        """
        Returns an astropy.table.Table object containing a position-position 
        catalog of leaves in a dendrogram.
        
        Parameters
        ----------
        dendrogram : ~astrodendro.dendrogram.Dendrogram object, optional
            The dendrogram object to extract sources from.
        """
        
        if not dendrogram:
            try:
                dendrogram = self.dendrogram
            except AttributeError:
                dendrogram = self.to_dendrogram()
                
        cat = pp_catalog(dendrogram.leaves, self.metadata)
        for i, idx in enumerate(cat['_idx']):
            cat['_idx'][i] = int('{:.0f}{:03d}'.format(
                                       np.round(self.nu.to(u.GHz).value), idx))
    
        try:
            cat['major_sigma'] = cat['major_sigma']*np.sqrt(8*np.log(2))
            cat['minor_sigma'] = cat['minor_sigma']*np.sqrt(8*np.log(2))
            cat.rename_column('major_sigma', 'major_fwhm')
            cat.rename_column('minor_sigma', 'minor_fwhm')
            cat.rename_column('flux', 'dend_flux_{}'.format(self.freq_id))
        except KeyError:
            pass
        
        try:
            cat.remove_column('rejected')
            cat.remove_column('detected_'+self.freq_id)
        except KeyError:
            pass
            
        cat.add_column(Column(np.zeros(len(cat)), dtype=int), name='rejected')
        cat.add_column(Column(np.ones(len(cat)), dtype=int), 
                       name='detected_'+self.freq_id)
                       
        self.catalog = Table(cat, masked=True)
        
        return Table(cat, masked=True)


    def _make_cutouts(self, catalog=None, data=None, save=True):
        """
        Make a cutout_data of cutout regions around all source centers in the 
        catalog.
        
        Parameters
        ----------
        save : bool, optional
            If enabled, the cutouts and cutout data will both be saved as 
            instance attributes. Default is True.
            
        Returns
        ----------
        List of astropy.nddata.utils.Cutout2D objects, list of cutout data
            
        """
        
        if catalog is None:   
            try:
                catalog = self.catalog
            except AttributeError:
                catalog = self.to_catalog()
                
        if data is None:
            data = self.data

        size = 2.2*(np.max(catalog['major_fwhm'])*u.deg 
            + self.annulus_padding 
            + self.annulus_width)

        cutouts = []
        cutout_data = []
        
        for i in range(len(catalog)):
            x_cen = catalog['x_cen'][i] * u.deg
            y_cen = catalog['y_cen'][i] * u.deg
            
            position = coordinates.SkyCoord(x_cen, 
                                            y_cen, 
                                            frame='icrs',
                                            unit=(u.deg, u.deg))
                                            
            pixel_position = np.array(position.to_pixel(self.wcs))
            
            try:
                cutout = Cutout2D(data, 
                                  position, 
                                  size, 
                                  self.wcs, 
                                  mode='partial')
                cutouts.append(cutout)
                cutout_data.append(cutout.data)
                
            except NoOverlapError:
                catalog['rejected'][i] = 1
                cutouts.append(float('nan'))
                cutout_data.append(float('nan'))
        
        cutouts = np.array(cutouts)
        cutout_data = np.array(cutout_data)
        
        if save:
            self._cutouts = cutouts
            self._cutout_data = cutout_data
        # NOTE: If 'sort' is called, the catalog's attributes also need to be
        # sorted accordingly. Might be tricky.
        
        return cutouts, cutout_data
    
    
    def get_pixels(self, aperture, catalog=None, data=None, cutouts=None, 
                   save=True, **kwargs):
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
        List of pixel arrays, list of masks
        """
        if catalog is None:
            try:
                catalog = self.catalog
            except AttributeError:
                catalog = self.to_catalog()
                
        if data is None:
            data = self.data
        
        if cutouts is None:
            size = 2.2*(np.max(catalog['major_fwhm'])*u.deg 
                        + self.annulus_padding 
                        + self.annulus_width)
            cutouts, cutout_data = self._make_cutouts(catalog=catalog,
                                                      data=data)
        
        pix_arrays = []
        masks = []
        
        for i in range(len(cutouts)):
        
            if isinstance(cutouts[i], Cutout2D):
                pass
            else:
                pix_arrays.append(float('nan'))
                masks.append(float('nan'))
                continue
            
            this_mask = aperture(catalog[i], 
                                 cutouts[i], 
                                 self,
                                 **kwargs)
                                    
            pix_arrays.append(cutouts[i].data[this_mask.astype('bool')])
            masks.append(this_mask)
        
        if save:
            self.__dict__['pixels_{}'
                          .format(aperture.__name__)] = np.array(pix_arrays)
            self.__dict__['mask_{}'
                          .format(aperture.__name__)] = np.array(masks)
        
        return np.array(pix_arrays), np.array(masks)
        
    
    def get_snr(self, source=None, background=None, catalog=None, data=None, 
                cutouts=None, cutout_data=None, peak=True, save=True):
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
        
        if catalog is None:
                try:
                    catalog = self.catalog
                except AttributeError:
                    catalog = self.to_catalog()
        
        # Cascade check
        if source is None or background is None:    
                    
            if data is None:
                data = self.data

            if cutouts is None or cutout_data is None:
                size = 2.2*(np.max(catalog['major_fwhm'])*u.deg 
                                    + self.annulus_padding 
                                    + self.annulus_width)
                cutouts, cutout_data = self._make_cutouts(catalog=catalog, 
                                                          data=data)
        
            background = self.get_pixels(annulus, 
                                         catalog=catalog,
                                         data=data,
                                         cutouts=cutouts)[0]
                                          
            source = self.get_pixels(ellipse, 
                                     catalog=catalog,
                                     data=data,
                                     cutouts=cutouts)[0]
        
        snr_vals = []
        for i in range(len(catalog)):
            try:
                snr = np.max(source[i]) / rms(background[i])
            except (ZeroDivisionError, ValueError) as e:
                snr = 0.0
            snr_vals.append(snr)
            
        if save:
            self.snr = np.array(snr_vals)
            try:
                catalog.remove_column('snr_'+self.freq_id)
            except KeyError:
                pass
            catalog.add_column(Column(snr_vals), name='snr_'+self.freq_id)
            
        return np.array(snr_vals)
    
    
    def plot_grid(self, catalog=None, data=None, cutouts=None, cutout_data=None,
                  apertures=None, skip=True, outfile=None):
        """
        Plot sources in a grid.
        
        Parameters
        ----------
        catalog : astropy.table.Table object, optional
            The catalog used to extract source positions.
        data : numpy.ndarray, optional
            The image data displayed and used to make cutouts.
        cutouts : list of astropy.nddata.utils.Cutout2D objects, optional
            Image cutout regions to save computation time, if they have already
            been calculated.
        cutout_data : list of numpy.ndarrays, optional
            Image cutout region data to save on computation time, if it has 
            already been calculated.
        apertures : list of dendrocat.aperture functions, optional
            Apertures to plot over the image cutouts.
        skip : bool, optional
            If enabled, don't plot rejected sources. Default is True.
        """
        
        if catalog is None:
            try:
                catalog = self.catalog
            except AttributeError:
                catalog = self.to_catalog()
        
        if data is None:
            data = self.data
        
        # Make sure ellipse and annulus apertures are used
        if apertures is None:
            apertures = [ellipse, annulus]
        else:
            apertures = list(set(apertures+[ellipse, 
                                            annulus]))
        
        # Get cutouts
        if cutouts is None or cutout_data is None:
            cutouts, cutout_data = self._make_cutouts(catalog=catalog, 
                                                      data=data)
        
        # Get pixels and masks in each aperture
        ap_names = []
        pixels = []
        masks = []
        
        for aperture in apertures:
            some_pixels, a_mask = self.get_pixels(aperture,
                                                  catalog=catalog,
                                                  data=data,
                                                  cutouts=cutouts)
            ap_names.append(aperture.__name__)
            pixels.append(some_pixels)
            masks.append(a_mask)
        
        # Find SNR
        ellipse_pix = pixels[ap_names.index('ellipse')]
        annulus_pix = pixels[ap_names.index('annulus')]
        
        snr_vals = self.get_snr(source=ellipse_pix, background=annulus_pix,
                                catalog=catalog)
        
        names = np.array(catalog['_idx'])
        rejected = np.array(catalog['rejected'])
        
        if skip:
            accepted_indices = np.where(catalog['rejected'] == 0)[0]
            snr_vals = snr_vals[accepted_indices]
            cutout_data = cutout_data[accepted_indices]
            names = names[accepted_indices]
            rejected = rejected[accepted_indices]
            for k in range(len(masks)):
                masks[k] = masks[k][accepted_indices]
        
        an = np.ones(len(cutouts), dtype='bool')
        for i in range(len(cutouts)):
            try:
                np.isnan(cutouts[i])
                an[i] = 0
            except TypeError:
                pass
        snr_vals = snr_vals[an]
        cutout_data = cutout_data[an]
        for k in range(len(masks)):
            masks[k] = masks[k][an]
        names = names[an]
        rejected = rejected[an]
        
        n_images = len(cutout_data)
        xplots = int(np.around(np.sqrt(n_images)))
        yplots = xplots + 1
        gs1 = gs.GridSpec(yplots, xplots, wspace=0.0, hspace=0.0)
        plt.figure(figsize=(9.5, 10))
        
        for i in range(n_images):
            
            image = cutout_data[i]
            ax = plt.subplot(gs1[i])
            
            if rejected[i] == 1:
                plt.imshow(image, origin='lower', cmap='gray')
            else:
                plt.imshow(image, origin='lower')

            for k in range(len(masks)):
                plt.imshow(masks[k][i], origin='lower', cmap='gray', 
                           alpha=0.15)
                
            plt.text(0, 0, 'SN {:.1f}'.format(snr_vals[i]), fontsize=7, 
                     color='w', ha='left', va='bottom', 
                     transform=ax.transAxes)
            plt.text(0, 1, str(names[i]), fontsize=7, color='w', ha='left', 
                     va='top', transform=ax.transAxes)

            plt.xticks([])
            plt.yticks([])
        
        plt.tight_layout()
        
        if outfile is not None:
            plt.savefig(outfile, dpi=300, bbox_inches='tight')                              

    
    def autoreject(self, threshold=None):
        """
        Reject noisy detections.
        
        Parameters
        ----------
        threshold : float, optional
            The signal-to-noise threshold below which sources are rejected
        """
            
        if threshold is None:
            threshold = self.threshold
      
        try:
            snrs = self.snr
        except AttributeError:
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
            self.catalog['rejected'][np.where(self.catalog['_idx'] == idx)] = 1
            
            
    def accept(self, accepted_list):
        
        for idx in accepted_list:
            self.catalog['rejected'][np.where(self.catalog['_idx'] == idx)] = 0

    def reset(self):
        self.catalog['rejected'] = 0
        
