from astropy.io import fits
from astropy import wcs
import radio_beam
import numpy as np
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D, NoOverlapError
from astropy.table import Column, Table, vstack
from astrodendro import Dendrogram, pp_catalog
import regions
import pickle
from copy import deepcopy
import warnings
warnings.filterwarnings('ignore')

if __package__ == '':
    __package__ = 'dendrocat'
from .aperture import Aperture, Ellipse, Circle, Annulus
from .utils import rms, ucheck

class UnknownApertureError(Exception):
    pass

class RadioSource:
    """
    An object to store radio image data.
    """

    def __init__(self, hdu,  name=None, freq_id=None):
        """
        Parameters
        ----------
        hdu : `~astropy.io.fits.hdu.image.PrimaryHDU`
            An astropy FITS HDU object containing the radio image data and
            header.
        region_id : str, optional
            An identifier specifying what sky object the radio image contains.
        freq_id : str, optional
            An identifier specifying the observation frequency (Ex: 226.0GHz).
            If not specified, it will be generated from the FITS image header.
        """
        self.hdu = hdu
        self.header = hdu[0].header
        self.data = hdu[0].data.squeeze()
        self.freq_id = freq_id

        self.__name__ = name

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
                    self.freq_id = ('{:.1f}'.format(self.nu
                                            .to(u.GHz)).replace(' ', ''))
                    if self.__name__ is None:
                        self.__name__ = 'Unknown_{}'.format(self.freq_id)
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
        """
        Sets RadioSource metadata using nu, WCS, and other FITS header data.
        """

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

        Parameters
        ----------
        min_value : float, optional
            Minimum detection level to be considered in the dendrogram.
        min_delta : float, optional
            How significant a dendrogram structure has to be in order to be
            considered a separate entity.
        min_npix : float, optional
            Minimum number of pixel needed for a dendrogram structure to be
            considered a separate entity.
        save : bool, optional
            If enabled, the resulting dendrogram will be saved as an instance
            attribute. Default is True.

        Returns
        ----------
        `~astrodendro.dendrogram.Dendrogram` object
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
        Creates a new position-position catalog of leaves in a dendrogram.
        This task will overwrite the existing catalog if there is one.

        Parameters
        ----------
        dendrogram : `~astrodendro.dendrogram.Dendrogram` object, optional
            The dendrogram object to extract sources from.

        Returns
        -------
        `~astropy.table.Table`
        """

        if not dendrogram:
            try:
                dendrogram = self.dendrogram
            except AttributeError:
                dendrogram = self.to_dendrogram()

        cat = pp_catalog(dendrogram.leaves, self.metadata)
        strarr=[str('{:.0f}{:04d}'.format(np.round(self.nu.to(u.GHz).value), idx)) for idx in cat['_idx']]  
        cat.add_column(Column(data=strarr), name='_name')
        cat.add_column(Column(data=range(len(cat))), name='_index')
        cat = cat[sorted(cat.colnames)]

        try:
            cat['major_sigma'] = cat['major_sigma']*np.sqrt(8*np.log(2))
            cat['minor_sigma'] = cat['minor_sigma']*np.sqrt(8*np.log(2))
            cat.rename_column('major_sigma', 'major_fwhm')
            cat.rename_column('minor_sigma', 'minor_fwhm')
            cat.rename_column('flux', '{}_dend_flux'.format(self.freq_id))
        except KeyError:
            pass

        try:
            cat.remove_column('rejected')
            cat.remove_column(self.freq_id+'_detected')
        except KeyError:
            pass

        cat.add_column(Column(np.zeros(len(cat)), dtype=int), name='rejected')
        cat.add_column(Column(np.ones(len(cat)), dtype=int),
                       name=self.freq_id+'_detected')

        self.catalog = Table(cat, masked=True)
        return Table(cat, masked=True)


    def add_sources(self, *args):
        """
        Adds external source entries to the existing catalog.

        Parameters
        ----------
        *args: `~astropy.table.Table`
            A source catalog containing the sources you wish to add to the
            existing catalog.
        """

        for sources in args:
            self.catalog = vstack([self.catalog, sources])
            self.catalog['_index'] = range(len(self.catalog))


    def _make_cutouts(self, catalog=None, data=None, save=True):
        """
        Make a cutout of cutout regions around all source centers in the
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

        size = 0.7*(np.max(catalog['major_fwhm'])*u.deg
                + self.annulus_padding
                + self.annulus_width)

        cutouts = []
        cutout_data = []

        for i in range(len(catalog)):
            x_cen = catalog['x_cen'][i] * u.deg
            y_cen = catalog['y_cen'][i] * u.deg

            position = coordinates.SkyCoord(x_cen,
                                            y_cen,
                                            frame=wcs.utils.wcs_to_celestial_frame(self.wcs).name,
                                            unit=(u.deg, u.deg))

            # commented out b/c not used
            # pixel_position = np.array(position.to_pixel(self.wcs))

            try:
                cutout = Cutout2D(data,
                                  position,
                                  size,
                                  wcs=self.wcs,
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
                   save=True):
        """
        Get pixels within an aperture for each entry in the specified catalog.

        Parameters
        ----------
        aperture: `~dendrocat.aperture.Aperture`
            The aperture determining which pixels to grab.
        catalog: `~astropy.table.Table`, optional
            A source catalog containing the center positions of each source.
        data: array-like
            Image data for the sources in the catalog.
        cutouts:
            For developer use

        Returns
        -------
        pixels, masks
        `~numpy.ndarray`, `~numpy.ndarray`
        """

        if catalog is None:
            try:
                catalog = self.catalog
            except AttributeError:
                catalog = self.to_catalog()

        if data is None:
            data = self.data

        if cutouts is None:
            cutouts, cutout_data = self._make_cutouts(catalog=catalog,
                                                      data=data)
        aperture_original = deepcopy(aperture)
        pix_arrays = []
        masks = []

        for i in range(len(cutouts)):

            if isinstance(cutouts[i], Cutout2D):
                pass
            else:
                pix_arrays.append(float('nan'))
                masks.append(float('nan'))
                continue

            frame = wcs.utils.wcs_to_celestial_frame(cutouts[i].wcs).name

            x_cen = catalog['x_cen'][i]
            y_cen = catalog['y_cen'][i]
            major = catalog['major_fwhm'][i]
            minor = catalog['minor_fwhm'][i]
            pa = catalog['position_angle'][i]

            if isinstance(aperture, Aperture):
                # If this is the case, then aperture has already been given
                # parameters. It should be 'fixed' dimensions. We just need to
                # replace the center value with the centers from the sources.

                if aperture.unit.is_equivalent(u.deg):
                    aperture.center = coordinates.SkyCoord(x_cen*u.deg,
                                                           y_cen*u.deg,
                                                           frame=frame)
                elif aperture.unit.is_equivalent(u.pix):
                    sky = coordinates.SkyCoord(x_cen*u.deg, y_cen*u.deg,
                                               frame=frame)
                    pixel = ucheck(sky.to_pixel(cutouts[i].wcs), u.pix)
                    aperture.center = pixel
                    aperture.x_cen, aperture.y_cen = pixel[0], pixel[1]

            elif issubclass(aperture, Aperture):
                # If this is the case, then the aperture type has been
                # specified and doesn't have any parameters associated to it.

                # DEFAULTS FOR VARIABLE APERTURES STORED HERE
                cen = [x_cen, y_cen]
                if aperture == Ellipse:
                    aperture = Ellipse(cen, major, minor, pa, unit=u.deg,
                                       frame=frame)

                elif aperture == Annulus:
                    inner_r = major*u.deg+self.annulus_padding
                    outer_r = major*u.deg+self.annulus_padding+self.annulus_width
                    aperture = Annulus(cen, inner_r, outer_r, unit=u.deg, frame=frame)

                elif aperture == Circle:
                    radius = major
                    aperture = Circle(cen, radius, unit=u.deg, frame=frame)

                else:
                    raise UnknownApertureError('Aperture not recognized. Pass'
                                               ' an instance of a custom aper'
                                               'ture instead.')

            this_mask = aperture.place(cutouts[i].data, wcs=cutouts[i].wcs)
            if this_mask.sum() == 0:
                raise ValueError("No pixels within aperture")
            pix_arrays.append(cutouts[i].data[this_mask])
            masks.append(this_mask)
            aperture = aperture_original # reset the aperture for the next source

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
        source: array-like
            Array of source fluxes to use in SNR calculation.
        background: array-like
            Array of background fluxes to use in SNR calculation.
        catalog: `~astropy.table.Table`
            The catalog of sources for which to calculate the SNR.
        data: array-like
            Image data for the sources in the catalog.
        cutouts:
            For debugging. Provides a specific set of cutouts instead of 
            letting the function generate them.
        cutout_data:
            For debugging. Provides a specific set of cutout data instead of 
            letting the function generate it.
        peak : bool, optional
            Use peak flux of source pixels as 'signal'. Default is True.
        save : bool, optional
            If enabled, the snr will be saved as a column in the source catalog
            and as an instance attribute. Default is True.

        Returns
        -------
        `~numpy.ndarray`

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
                #size = 2.2*(np.max(catalog['major_fwhm'])*u.deg
                #            + self.annulus_padding
                #            + self.annulus_width)
                cutouts, cutout_data = self._make_cutouts(catalog=catalog,
                                                          data=data)

            background = self.get_pixels(Annulus,
                                         catalog=catalog,
                                         data=data,
                                         cutouts=cutouts)[0]

            source = self.get_pixels(Ellipse,
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
                catalog.remove_column(self.freq_id+'_snr')
            except KeyError:
                pass
            catalog.add_column(Column(snr_vals), name=self.freq_id+'_snr')

        return np.array(snr_vals)


    def plot_grid(self, catalog=None, data=None, cutouts=None,
                  cutout_data=None, source_aperture=None, bkg_aperture=None,
                  skip_rejects=True, outfile=None, figurekwargs={}):
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
        skip_rejects : bool, optional
            If enabled, don't plot rejected sources. Default is True.
        """

        import matplotlib.gridspec as gs
        import matplotlib.pyplot as plt

        if catalog is None:
            try:
                catalog = self.catalog
            except AttributeError:
                catalog = self.to_catalog()

        if data is None:
            data = self.data

        if source_aperture is None:
            source_aperture = Ellipse

        if bkg_aperture is None:
            bkg_aperture = Annulus

        # Get cutouts
        if cutouts is None or cutout_data is None:
            cutouts, cutout_data = self._make_cutouts(catalog=catalog,
                                                      data=data)

        # Get pixels and masks in each aperture
        ap_names = []
        pixels = []
        masks = []

        for aperture in [source_aperture, bkg_aperture]:
            some_pixels, a_mask = self.get_pixels(aperture, catalog=catalog,
                                                  data=data, cutouts=cutouts)
            ap_names.append(aperture.__name__)
            pixels.append(some_pixels)
            masks.append(a_mask)

        # Find SNR
        ellipse_pix = pixels[0]
        annulus_pix = pixels[1]

        snr_vals = self.get_snr(source=ellipse_pix, background=annulus_pix,
                                catalog=catalog)

        names = np.array(catalog['_name'])
        rejected = np.array(catalog['rejected'])

        if skip_rejects:
            accepted_indices = np.where(catalog['rejected'] == 0)[0]
            snr_vals = snr_vals[accepted_indices]
            cutout_data = cutout_data[accepted_indices]
            cutouts = cutouts[accepted_indices]
            names = names[accepted_indices]
            rejected = rejected[accepted_indices]
            for k in range(len(masks)):
                masks[k] = masks[k][accepted_indices]

        an = np.ones(len(cutouts), dtype='bool')
        for i in range(len(cutouts)):
            try:
                # check whether cutouts[i] is a cutout or is NaN
                np.isnan(cutouts[i])
                an[i] = False
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
        plt.figure(figsize=(9.5, 10), **figurekwargs)

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
            plt.text(0, 1, names[i], fontsize=7, color='w', ha='left',
                     va='top', transform=ax.transAxes)

            plt.xticks([])
            plt.yticks([])

        plt.tight_layout()

        if outfile is not None:
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
        else:
            plt.show()

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

        snrs = self.get_snr()

        try:
            self.catalog['rejected'] = np.zeros(len(self.catalog), dtype=int)
        except KeyError:
            self.catalog.add_column(Column(np.zeros(len(self.catalog))),
                                    name='rejected')

        for i in range(len(self.catalog)):
            if snrs[i] <= threshold or np.isnan(snrs[i]):
                self.catalog['rejected'][i] = 1

        self.accepted = self.catalog[self.catalog['rejected']==0]
        self.rejected = self.catalog[self.catalog['rejected']==1]


    def reject(self, rejected_list):
        """
        Reject specific sources in the catalog.

        Parameters
        ----------
        rejected_list: list
            A list of ``_name``s, for which each corresponding entry will be
            marked rejected.
        """
        rejected_list = np.array(rejected_list, dtype=str)
        for nm in rejected_list:
            self.catalog['rejected'][np.where(self.catalog['_name'] == nm)] = 1
        self.accepted = self.catalog[self.catalog['rejected']==0]
        self.rejected = self.catalog[self.catalog['rejected']==1]

    def accept(self, accepted_list):
        """
        Accept specific sources in the catalog.

        Parameters
        ----------
        accpeted_list: list
            A list of ``_name``s, for which each corresponding entry will be
            marked accepted.
        """
        accepted_list = np.array(accepted_list, dtype=str)
        for nm in accepted_list:
            self.catalog['rejected'][np.where(self.catalog['_name'] == nm)] = 0
        self.accepted = self.catalog[self.catalog['rejected']==0]
        self.rejected = self.catalog[self.catalog['rejected']==1]

    def reset(self):
        """
        Reset all sources' rejection flags to 0 (all accepted).
        """

        self.catalog['rejected'] = 0
        self.accepted = self.catalog[self.catalog['rejected']==0]
        self.rejected = self.catalog[self.catalog['rejected']==1]

    def grab(self, name, skip_rejects=False):
        """
        Search the catalog for an entry matching a specific name, and return
         it.

        Parameters
        ----------
        name: tuple, list, or str
            The name or names of the sources to search for.
        skip_rejects: bool, optional
            If enabled, will only search accepted sources.
        """


        if skip_rejects:
            catalog = self.accepted
        else:
            catalog = self.catalog

        if type(name) == tuple or type(name) == list:
            name = np.array(name).astype(str)
            indices = []
            for i in range(len(catalog)):
                if catalog['_name'][i] in names:
                    indices.append(i)
            indices = np.array(indices)
            return catalog[indices]
        else:
            return self.catalog[self.catalog['_name']==str(name)]


    def dump(self, outfile):
        """
        Dump the `~dendrocat.RadioSource` object via pickle.

        Parameters
        ----------
        outfile : str
            Desired output file path.
        """
        outfile = outfile.split('.')[0]+'.pickle'
        with open(outfile, 'wb') as output:
            pickle.dump(obj, output, protocol=pickle.HIGHEST_PROTOCOL)
