from astropy.table import MaskedColumn, vstack, hstack, Table
from collections import OrderedDict
import numpy as np
import astropy.units as u
from copy import deepcopy
from astropy.coordinates import SkyCoord, Angle

if __package__ == '':
    __package__ = 'dendrocat'
from .utils import rms, specindex, ucheck
from .radiosource import RadioSource
from .aperture import Aperture

class ApertureError(Exception):
    pass

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
            self.accepted = catalog[catalog['rejected']==0]
        self.add_objects(*args)


    def grab(self, name, skip_rejects=False):
        """
        Grab a source or sources by name.

        Parameters
        ----------
        name : str or list
            String or list of strings to search the catalog "_name" header for.
        skip_rejects : bool, optional
            If enabled, rejected sources will not be queried. Disabled by
            default.
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


    def add_objects(self, *args):
        """
        Add a new `~dendrocat.RadioSource` object to the existing master
        catalog.

        Parameters
        ----------
        *args : '~dendrocat.RadioSource` objects
            RadioSource objects to add to the master catalog.
        """
        if not hasattr(self, 'other_catalogs'):
            self.other_catalogs = []

        for obj in args:
            if isinstance(obj, MasterCatalog):
                for key, value in obj.__dict__.items():
                    if isinstance(value, RadioSource):
                        self.__dict__[key] = value
            else:
                objname = [k for k, v in locals().items() if v is obj][0]
                self.other_catalogs.append(obj)
                setattr(self, objname, obj)


    def add_sources(self, *args):
        """
        Add source entries from another catalog.

        Parameters
        ----------
        *args : astropy.table.Table objects
            Source tables to vertically stack with the existing master catalog.
        """
        for sources in args:
            self.catalog = vstack([self.catalog, sources])
            self.catalog['_index'] = range(len(self.catalog))


    def photometer(self, *args, catalog=None):
        """
        Add photometry data columns to the master catalog.

        Parameters
        ----------
        args : `~dendrocat.Aperture` objects
            The apertures to use for photometry. Can be given as either
            instances or objects, to use fixed or variable aperture widths,
            respectively.

        catalog : astropy.table.Table object
            The catalog from which to extract source coordinates and ellipse
            parameters.
        """

        for aperture in args:
            if catalog is None:
                catalog = self.catalog

            rs_objects = []
            for i, obj in enumerate(self.__dict__.values()):
                if isinstance(obj, RadioSource):
                    rs_objects.append(obj)

            for i, rs_obj in enumerate(rs_objects):

                data = rs_obj.data
                cutouts, cutout_data = rs_obj._make_cutouts(catalog=catalog,
                                                            save=False)

                pix_in_aperture = rs_obj.get_pixels(
                                                    aperture,
                                                    catalog=catalog,
                                                    data=data,
                                                    cutouts=cutouts,
                                                    )[0]

                names = [
                    rs_obj.freq_id+'_'+aperture.__name__+'_peak',
                    rs_obj.freq_id+'_'+aperture.__name__+'_sum',
                    rs_obj.freq_id+'_'+aperture.__name__+'_rms',
                    rs_obj.freq_id+'_'+aperture.__name__+'_median',
                    rs_obj.freq_id+'_'+aperture.__name__+'_npix'
                ]

                peak_data = np.zeros(len(pix_in_aperture))
                for j in range(len(pix_in_aperture)):
                    try:
                        peak_data[j] = np.max(pix_in_aperture[j])
                    except ValueError:
                        peak_data[j] = float('nan')
                aperture_peak_col = MaskedColumn(data=peak_data,
                                                 name=names[0])

                sum_data = np.zeros(len(pix_in_aperture))
                for j in range(len(pix_in_aperture)):
                    ind = [pix_in_aperture[j] > 0.]
                    try:
                        sum_data[j] = np.sum(pix_in_aperture[j][ind])/rs_obj.ppbeam
                    except TypeError: # Catches single pixel apertures
                        sum_data[j] = float('nan')
                aperture_sum_col = MaskedColumn(data=sum_data,
                                                name=names[1])

                rms_data = np.zeros(len(pix_in_aperture))
                for j in range(len(pix_in_aperture)):
                    rms_data[j] = rms(pix_in_aperture[j])
                aperture_rms_col = MaskedColumn(data=rms_data,
                                                name=names[2])

                median_data = np.zeros(len(pix_in_aperture))
                for j in range(len(pix_in_aperture)):
                    median_data[j] = np.median(pix_in_aperture[j])
                aperture_median_col = MaskedColumn(data=median_data,
                                                   name=names[3])

                npix_data = np.zeros(len(pix_in_aperture))
                for j in range(len(pix_in_aperture)):
                    if np.isnan(pix_in_aperture[j]).any():
                        npix_data[j] = float('nan')
                    else:
                        npix_data[j] = len(pix_in_aperture[j])
                    aperture_npix_col = MaskedColumn(data=npix_data,
                                                     name=names[4])

                try:
                    self.catalog.remove_columns(names)
                except KeyError:
                    pass

                self.catalog.add_columns([
                    aperture_peak_col,
                    aperture_sum_col,
                    aperture_rms_col,
                    aperture_median_col,
                    aperture_npix_col
                ])

                # Mask NaN values
                for col in self.catalog.colnames:
                    try:
                        isnan = np.argwhere(np.isnan(list(self.catalog[col])))
                        self.catalog.mask[col][isnan] = True
                    except TypeError:
                        pass


    def ffplot(self, rsobj1, rsobj2, apertures=[], bkg_apertures=[],
               alphas=None, peak=False, label=False, log=True, outfile=None):

        """
        Produce a flux-flux plot for two `~dendrocat.RadioSource` objects.

        Parameters
        ----------
        rsobj1 : `~dendrocat.RadioSource` object
            One of two radio source objects from which to make a flux-flux
            plot.
        rsobj2 : `~dendrocat.RadioSource` object
            The other of two radio source objects from which to make a
            flux-flux plot.
        apertures : list
            List of `~dendrocat.Aperture` objects to use for source apertures.
        bkg_apertures : list
            List of `~dendrocat.Aperture` objects to use for background
            apertures.
        alphas : list, optional
            Spectral indices to overplot on top of the flux-flux data. 1, 2,
            and 3 will be used by default.
        peak : bool, optional
            If enabled, peak flux inside the aperture is used instead of
            aperture sum. Disabled by default.
        label : bool, optional
            If enabled, labels will be printed on the plot to identify sources.
            Disabled by default.
        log : bool, optional
            If enabled, results will be shown on log-log axes. Enabled by
            default.
        outfile : str, optional
            If provided, output plot will be saved to this file path.
        """

        import matplotlib.pyplot as plt

        if type(apertures) != list:
            apertures = list([apertures])

        if type(bkg_apertures) != list:
            bkg_apertures = list([bkg_apertures])

        if len(bkg_apertures) != len(apertures):
            raise ApertureError('Must give equal number of apertures and '
                                'background apertures')

        if rsobj1.nu > rsobj2.nu:
            rsobj1, rsobj2 = rsobj2, rsobj1

        catalog = deepcopy(self.catalog)
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        if alphas is None:
            alphas = [1, 2, 3]

        cols = []
        for aperture in apertures:

            if peak:
                cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_peak')
                cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_peak')
            else:
                cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_sum')
                cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_sum')

        for bkg_aperture in bkg_apertures:
            cols.append(rsobj1.freq_id+'_'+bkg_aperture.__name__+'_rms')
            cols.append(rsobj2.freq_id+'_'+bkg_aperture.__name__+'_rms')

        try:
            index = list(set(range(len(catalog)))^
                         set(np.nonzero(catalog.mask[cols])[0])
                         .union(set(np.where(catalog['rejected']==1)[0])))
        except KeyError:
            for aperture in apertures:
                self.photometer(aperture)
            catalog = deepcopy(self.catalog)
            index = list(set(range(len(catalog)))^
                         set(np.nonzero(catalog.mask[cols])[0])
                         .union(set(np.where(catalog['rejected']==1)[0])))

        catalog = catalog[index]

        flux1 = []
        flux2 = []
        err1 = []
        err2 = []

        for aperture in apertures:

            if peak:
                flux1.append(catalog['{}_{}_peak'.format(rsobj1.freq_id,
                                                         aperture.__name__)])
                flux2.append(catalog['{}_{}_peak'.format(rsobj2.freq_id,
                                                         aperture.__name__)])
            else:
                flux1.append(catalog['{}_{}_sum'.format(rsobj1.freq_id,
                                                         aperture.__name__)])
                flux2.append(catalog['{}_{}_sum'.format(rsobj2.freq_id,
                                                         aperture.__name__)])

        for bkg_aperture in bkg_apertures:
            err1.append(catalog[rsobj1.freq_id+bkg_aperture.__name__+'_rms'])
            err2.append(catalog[rsobj2.freq_id+bkg_aperture.__name__+'_rms'])

        marker_labels = catalog['_name']

        xflux = np.linspace(np.min(flux1), np.max(flux1), 10)
        yfluxes = []

        for alpha in alphas:
            yfluxes.append(specindex(rsobj1.nu, rsobj2.nu, xflux, alpha))

        n_images = len(apertures)
        xplots = int(np.around(np.sqrt(n_images)))
        yplots = xplots
        fig, axes = plt.subplots(ncols=yplots, nrows=xplots, figsize=(12, 12))

        for i in range(len(apertures)-1):
            ax = np.ndarray.flatten(np.array(axes))[i]
            ax.errorbar(flux1[i], flux2[i], xerr=err1[i], yerr=err2[i], ms=2,
                        alpha=0.75, elinewidth=0.5, color=colors[i], fmt='o',
                        label='{} Aperture Sums'.format(apertures[i].__name__))

            for j, yflux in enumerate(yfluxes):
                ax.plot(xflux, yflux, '--', color=colors[np.negative(j)],
                        label='Spectral Index = {}'.format(alphas[j]))

            ax.set_xticks([])
            ax.set_yticks([])

            if label:
                for j, label in enumerate(marker_labels):
                    ax.annotate(label,
                                (flux1[i][j], flux2[i][j]),
                                size=8)
            if peak:
                if log:
                    plt.xlabel('Log Peak Flux {}'.format(rsobj1.freq_id))
                    plt.ylabel('Log Peak Flux {}'.format(rsobj2.freq_id))
                    plt.xscale('log')
                    plt.yscale('log')
                else:
                    plt.xlabel('Peak Flux {}'.format(rsobj1.freq_id))
                    plt.ylabel('Peak Flux {}'.format(rsobj2.freq_id))
            else:
                if log:
                    plt.xlabel('Log Flux {}'.format(rsobj1.freq_id))
                    plt.ylabel('Log Flux {}'.format(rsobj2.freq_id))
                    plt.xscale('log')
                    plt.yscale('log')
                else:
                    plt.xlabel('Flux {}'.format(rsobj1.freq_id))
                    plt.ylabel('Flux {}'.format(rsobj2.freq_id))
            plt.suptitle('{} Flux v. {} Flux'.format(rsobj2.freq_id,
                                                     rsobj1.freq_id))
            plt.legend()

            if outfile is not None:
                plt.savefig(outfile, dpi=300, bbox_inches='tight')
