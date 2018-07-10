from astropy.io import fits
from astropy import wcs
import radio_beam
import numpy as np
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D
from astropy.table import Column, Table, MaskedColumn, vstack
from astrodendro import Dendrogram, pp_catalog
from copy import deepcopy
from func import mask, rms, commonbeam
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import regions
from astropy.utils.console import ProgressBar
import warnings
warnings.filterwarnings('ignore')


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
    
    
    def to_catalog(self, dendrogram=None, catalog=None):
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
        
        cat.add_column(Column(np.ones(len(cat)), dtype=int), 
                       name='detected_'+self.freq_id)
                       
        self.catalog = Table(cat, masked=True)
        
        return Table(cat, masked=True)


    def _make_cutouts(self, sidelength, save=True):
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
    
    
    def save_ds9_regions(self, outfile, hide_rejects=True):
    
        if hide_rejects:
            catalog = self.catalog[np.where(self.catalog['rejected'] == 0)]
        else:
            catalog = self.catalog
            
        with open(outfile, 'w') as fh:
            fh.write("icrs\n")
            for row in catalog:
                fh.write("ellipse({x_cen}, {y_cen}, {major_fwhm}, " \
                         "{minor_fwhm}, {position_angle}) # text={{{_idx}}}\n"
                         .format(**dict(zip(row.colnames, row))))
    
    
    def get_pixels_annulus(self, padding=None, width=None, save=True):
        """
        Return a list of pixel arrays, each of which contains the pixels in
        an annulus of constant width and variable radius depending on the 
        major fwhm of the source.
        
        Parameters
        ----------
        padding : astropy.units.deg, optional
            The additional spacing between the major fwhm of the source and
            the inner radius of the annulus.
        width : astropy.units.deg, optional
            The width of the annulus, in degrees.
        save : bool, optional
            If enabled, the pixel arrays and masks will both be saved as 
            instance attributes. Default is True.
            
        Returns
        ----------
        List of pixel arrays
        """
        
        if not padding:
            padding = self.default_annulus_padding
        
        if not width:
            width = self.default_annulus_width
        
        try:
            catalog = self.catalog
        except AttributeError:
            catalog = self.to_catalog()
        
        size = 2.2*(np.max(catalog['major_fwhm'])*u.deg + padding + width)
        cutouts, cutout_data = self._make_cutouts(size)
        
        pix_arrays = []
        masks = []
        
        for i in range(len(cutouts)):
            center = regions.PixCoord(cutouts[i].center_cutout[0], 
                                      cutouts[i].center_cutout[1])
            
            inner_r = self.catalog[i]['major_fwhm']*u.deg + padding
            outer_r = inner_r + width
            
            innerann_reg = regions.CirclePixelRegion(center, 
                                                     inner_r/self.pixel_scale)
            outerann_reg = regions.CirclePixelRegion(center, 
                                                     outer_r/self.pixel_scale)
            
            annulus_mask = (mask(outerann_reg, cutouts[i]) 
                            - mask(innerann_reg, cutouts[i]))
            
            pix_arrays.append(cutouts[i].data[annulus_mask.astype('bool')])
            masks.append(annulus_mask)
        
        if save:
            self.pixels_in_annulus = pix_arrays
            self.mask_annulus = masks
        
        return pix_arrays


    def get_pixels_ellipse(self, save=True):
        """
        Return a list of pixel arrays, each of which contains the pixels in
        the source ellipses.
        
        Parameters
        ----------
        save : bool, optional
            If enabled, the pixel arrays and masks will both be saved as 
            instance attributes. Default is True.
            
        Returns
        ----------
        List of pixel arrays
        """
        cutouts = self._cutouts
        # Currently, get_pixels_annulus needs to be run first to set the size
        # of the cutouts and save them as an instance attribute. Not ideal.
        
        pix_arrays = []
        masks = []
        
        for i in range(len(cutouts)):
            center = regions.PixCoord(cutouts[i].center_cutout[0], 
                                      cutouts[i].center_cutout[1])
            
            pix_major = self.catalog[i]['major_fwhm']*u.deg / self.pixel_scale
            pix_minor = self.catalog[i]['minor_fwhm']*u.deg / self.pixel_scale
            pa = self.catalog[i]['position_angle']*u.deg
            
            radius = self.catalog[i]['major_fwhm'] * u.deg
            reg = regions.EllipsePixelRegion(center, pix_major, pix_minor, 
                                             angle=pa)                 
            ellipse_mask = mask(reg, cutouts[i])
            pix_arrays.append(cutouts[i].data[ellipse_mask.astype('bool')])
            masks.append(ellipse_mask)
        
        if save:
            self.pixels_in_ellipse = pix_arrays
            self.mask_ellipse = masks
        
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
                pixels_in_background = self.pixels_in_annulus
            except:
                pixels_in_background = self.get_pixels_annulus(
                                          padding=self.default_annulus_padding,
                                          width=self.default_annulus_width
                                        )
                                                               
        if not pixels_in_source:
            try:
                pixels_in_source = self.pixels_in_ellipse
            except AttributeError:
                pixels_in_source = self.get_pixels_ellipse()
        
        snr_vals = []
        for i in range(len(self.catalog)):
            snr = np.max(pixels_in_source[i])/rms(pixels_in_background[i])
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
                           alpha=0.25)
                
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
    
    
    def match(self, other, verbose=True):
        """
        Find sources that match up between two radio images. 
        
        Parameters
        ----------
        other : rsprocess.RadioSource object or rsprocess.MasterCatalog object
            The catalog with which to compare radio sources.
            
        Returns
        ----------
        astropy.table.Table object
        """
      
        all_colnames = set(self.catalog.colnames + other.catalog.colnames)
        stack = vstack([self.catalog, other.catalog])
        
        all_colnames.add('_idy')
        stack.add_column(Column(range(len(stack)), name='_idy'))
        stack = stack[sorted(list(all_colnames))]
        
        rejected = np.where(stack['rejected'] == 1)[0]
        
        if verbose:
            print('Combining matches')
            pb = ProgressBar(len(stack) - len(rejected))
        
        i = 0
        while True:
            
            if i == len(stack) - 1:
                break
            
            if i in rejected:
                i += 1
                continue
            
            teststar = stack[i]
            delta_p = vstack([stack[:i], stack[i+1:]]
                             )['_idx', '_idy', 'x_cen', 'y_cen']
            delta_p['x_cen'] = np.abs(delta_p['x_cen'] - teststar['x_cen'])                
            delta_p['y_cen'] = np.abs(delta_p['y_cen'] - teststar['y_cen'])
            delta_p.sort('x_cen')
            
            threshold = 1e-5
            found_match = False
            
            dist_col = MaskedColumn(length=len(delta_p), name='dist', 
                                    mask=True)
            
            for j in range(10):
                dist_col[j] = np.sqrt(delta_p[j]['x_cen']**2. 
                                      + delta_p[j]['y_cen']**2)            
                if dist_col[j] <= threshold:
                    found_match = True
                    
            delta_p.add_column(dist_col)
            delta_p.sort('dist')
            
            if found_match:
                match_index = np.where(stack['_idy'] == delta_p[0]['_idy'])
                match = deepcopy(stack[match_index])
                stack.remove_row(match_index[0][0])
                
                # Find the common bounding ellipse
                new_x_cen = np.average([match['x_cen'], teststar['x_cen']])
                new_y_cen = np.average([match['y_cen'], teststar['y_cen']])
                
                # Find new ellipse properties
                new_maj, new_min, new_pa = commonbeam(
                                             float(match['major_fwhm']), 
                                             float(match['minor_fwhm']), 
                                             float(match['position_angle']),
                                             float(teststar['major_fwhm']),
                                             float(teststar['minor_fwhm']),
                                             float(teststar['position_angle'])
                                             )
                
                # Replace properties of test star
                stack[i]['x_cen'] = new_x_cen       
                stack[i]['y_cen'] = new_y_cen
                stack[i]['major_fwhm'] = new_maj.value
                stack[i]['minor_fwhm'] = new_min.value
                stack[i]['position_angle'] = new_pa.value
        
                # Replace masked data with available values from the match
                for k, masked in enumerate(stack.mask[i]):
                    colname = stack.colnames[k]
                    if masked:
                        stack[i][colname] = match[colname]
            i += 1
            if verbose:
                pb.update()
        
        # Fill masked detection column fields with 'False'
        for colname in stack.colnames:
            if colname.split('_')[0] == 'detected':
                stack[colname].fill_value = 0
        
        stack.remove_column('_idy')
        return MasterCatalog(self, other, catalog=stack)
        

class MasterCatalog:
    """
    An object to store combined data from two or more RadioSource objects.
    """
    
    def __init__(self, *args, catalog=None):
        """
        Create a new master catalog object.
        
        Parameters
        ----------
        catalog : astropy.table.Table object
            The master table from which to build the catalog object.
        *args : radiosource.RadioSource objects
            RadioSource objects from which the master table was built.
        """
        if catalog is not None:
            self.catalog = catalog
        
        obj_prefix = 'radiosource_'
        
        for obj in args:
        
            if isinstance(obj, MasterCatalog):
                for key in mc.__dict__.keys():
                    if key.split('_')[0] == obj_prefix:
                        self.__dict__[key] = obj[key]
            else:        
                self.__dict__[obj_prefix+obj.freq_id)] = obj
                
            # Now has support for n matches!
            
    
    
        
