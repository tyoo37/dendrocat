from astropy.table import MaskedColumn, vstack, hstack, Table
from collections import OrderedDict
import numpy as np
import astropy.units as u
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from astropy.coordinates import SkyCoord, Angle

if __package__ == '':
    __package__ = 'dendrocat'
from .utils import rms, _matcher, specindex, findrow, check_units
from .radiosource import RadioSource
from .aperture import ellipse, annulus


def match(*args, verbose=True):
    """
    Wrapper for the _matcher method in .utils, which does the heavy lifting.
    """
    current_arg = args[0]
    for i in range(len(args)-1):
        current_arg = MasterCatalog(current_arg, args[i+1], 
                                    catalog=_matcher(current_arg, args[i+1]))
    return current_arg


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
        
        NOTE: Currently only supports multiple objects if they're of different
              frequencies (the freq_id must be unique).
        """
        if catalog is not None:
            self.catalog = catalog
        self.add_objects(*args)
        self.catalog

    def add_objects(self, *args):
        obj_prefix = 'radiosource_'
        for obj in args:
            if isinstance(obj, MasterCatalog):
                for key in obj.__dict__.keys():
                    if key.split('_')[0]+'_' == obj_prefix:
                        self.__dict__[key] = obj[key]
            else:        
                self.__dict__[obj_prefix+obj.freq_id] = obj
    
    
    def add_sources(self, *args):
        for sources in args:
            self.catalog = vstack([self.catalog, sources])
            self.catalog['_index'] = range(len(self.catalog))
    
    
    def photometer(self, aperture, catalog=None, **kwargs):
        """
        Add photometry data columns to the master catalog.
        
        Parameters
        ----------
        aperture : utils aperture mask function
            The function that creates a mask in the shape of the desired
            aperture, given source information and a cutout
        catalog : astropy.table.Table object
            The catalog from which to extract source coordinates and ellipse
            parameters.
        """
        
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
                                                **kwargs
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
                peak_data[j] = np.max(pix_in_aperture[j])
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
            
            
    def ffplot(self, rsobj1, rsobj2, apertures=[], alphas=None, peak=False,
               label=False, log=True, outfile=None):
           
        if rsobj1.nu > rsobj2.nu:
            rsobj1, rsobj2 = rsobj2, rsobj1
        
        catalog = deepcopy(self.catalog)
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
        apertures = list(set(apertures)|{ellipse, annulus})
        if alphas is None:
            alphas = [1, 2, 3]
        
        cols = []
        for aperture in apertures:
            
            cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_rms')
            cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_rms')
            
            if peak:
                cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_peak')
                cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_peak')
            else:
                cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_sum')
                cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_sum')
        
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
        
        for aperture in list(set(apertures)^set([annulus])):
        
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
                                     
        err1.append(catalog[rsobj1.freq_id+'_annulus_rms'])
        err2.append(catalog[rsobj2.freq_id+'_annulus_rms'])
        
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
            
            #ax.set_xlim([.6*np.min(flux1), 1.4*np.max(flux1)])
            #ax.set_ylim([.1*np.min(flux2), 1.9*np.max(flux2)])
            
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


    def match_external(self, cat, ra=None, dec=None, freq=None, flux_sum=None, 
                       flux_peak=None, err=None, shape=None,
                       tolerance=0.1*u.arcsec, skip_rejects=True):
        '''
        A large and frivolous method that serves a very specific purpose.
        Splits an external catalog by frequency and match each of them to 
        MasterCatalog sources, collects photometry data from the external 
        catalog, formats it to use the proper MasterCatalog table headers, and 
        then adds it (horizontally) to the existing MasterCatalog. Unmatched 
        entries are masked.
        
        Parameters
        ----------
        catalog : astropy.table.Table object
            The external catalog to use for source matching.
        ra : string
            The name of the column for the right ascension coordinate in the
            external catalog.
        dec : string
            The name of the column for the declination coordinate in the 
            external catalog.
        freq : str or float
            Either the column name of the frequencies in the catalog by which 
            the catalog will be split up, or a float specifying the frequency 
            of every entry in the catalog in GHz.
        flux_sum : str
            The name of the column for the aperture flux sum in the external
            catalog.
        flux_peak : str
            The name of the column for the peak flux in the external catalog.
        err : str
            The name of the column for the error in the external catalog.
        shape : str
            The name of the dendrocat.aperture shape to 'imitate' when adding
            the new catalog columns. Must match the aperture shape of any 
            existing photometry from the MasterCatalog you might want to use.           
        tolerance : astropy.units.quantity.Quantity
            The maximum sky distance between a master catalog source and an
            external catalog source for a match to be made.
        skip_rejects : bool, optional
            If enabled, rejected sources from the MasterCatalog will not be
            matched. Default is True.
            
        Returns
        ----------
        astropy.table.Table
            A copy of the original MasterCatalog.catalog, extended by columns
            from the external catalog. Non-matched entries under external
            columns are masked.
        ''' 
        catalogs = []
        freq_ids = []
        
        if shape is None:
            shape = 'ellipse'
        
        if type(freq) is float or type(freq) is int: # not yet debugged
            f_GHz = check_units(freq, u.GHz)
            freq_id = '{:.0f}'.format(np.round(f_GHz)).replace(' ', '')
            freq_ids.append(freq_id)
            new_cat = Table(masked=True)
            ### Finish coding this
            
        elif type(freq) is str:
            for f in set(list(cat[freq])):
                catalog = cat[cat[freq]==f]
                newcat = Table(masked=True)
                f_GHz = check_units(f, u.GHz)
                freq_id = '{:.1f}'.format(np.round(f_GHz)).replace(' ', '')
                freq_ids.append(freq_id)
                
                for col in catalog.colnames:
                    if flux_sum is not None:
                        if flux_sum in col:
                            newsum = MaskedColumn(
                                         data=catalog[col],
                                         name=freq_id+'_'+shape+'_sum')
                            newcat.add_column(newsum)
                    if flux_peak is not None:
                        if flux_peak in col:
                            newpeak = MaskedColumn(
                                          data=catalog[col],
                                          name=freq_id+'_'+shape+'_peak')
                            newcat.add_column(newpeak)
                    if err is not None:
                        if err in col:
                            newerr = MaskedColumn(
                                         data=catalog[col],
                                         name=freq_id+'_'+shape+'_err')
                            newcat.add_column(newerr)
                    if ra is not None:
                        if ra in col:
                            newra = MaskedColumn(data=catalog[col],
                                                  name='x_cen')
                            newcat.add_column(newra)
                    if dec is not None:
                        if dec in col:
                            newdec = MaskedColumn(data=catalog[col], 
                                                   name='y_cen')
                            newcat.add_column(newdec)
                            
                catalogs.append(newcat)
            
        current_table = deepcopy(self.catalog)
        
        for i, catalog in enumerate(catalogs):
            
            freq_id = freq_ids[i]
            
            ra = check_units(catalog['x_cen'])
            dec = check_units(catalog['y_cen'])
            
            selfra = check_units(current_table['x_cen'])
            selfdec = check_units(current_table['y_cen'])
            
            coords = SkyCoord(ra=ra, dec=dec)
            selfcoords = SkyCoord(ra=selfra, dec=selfdec)
            
            try:
                idx, d2d, d3d = selfcoords.match_to_catalog_sky(coords)
            except IndexError:
                continue
            matched = [d2d < tolerance]
                               
            ext_colnames = [c for c in catalog.colnames if freq_id in c]
            
            external_table = Table(catalog[ext_colnames][idx], masked=True)
            external_table.mask[~matched[0]] = True
            current_table = hstack([current_table, external_table])
            
        return current_table
            
    def plotsedgrid(self, row, alphas=None, path=None):
        row = Table(row, masked=True)
        apname = 'ellipse'
        method = 'sum'
        
        rsobjs = []
        for i, obj in enumerate(self.__dict__.values()):
            if isinstance(obj, RadioSource):
                rsobjs.append(obj)
        rsnus = [rsobj.nu for rsobj in rsobjs]
        rsnus, rsobjs = [list(s) for s in zip(*sorted(zip(rsnus, rsobjs)))]
        
        freq_ids = []
        nus = []
        fluxcols = []
        errcols = []

        for col in self.catalog.colnames:
            if 'GHz' in col:
                freq_id = col.split('_')[0]
                if row.mask[col][0] == False:
                    if apname in col and method in col:
                        freq_ids.append(freq_id)
                        fluxcols.append(col)
                    if 'annulus' in col and 'rms' in col:
                        errcols.append(col)
                    if 'ellipse' in col and 'err' in col:
                        errcols.append(col)
                        
        freq_ids = list(OrderedDict.fromkeys(freq_ids))
        fluxcols = list(OrderedDict.fromkeys(fluxcols))
        errcols = list(OrderedDict.fromkeys(errcols))
        
        nus = [float(s.split('GHz')[0]) for s in freq_ids]
        nus, sort = [list(s) for s in zip(*sorted(zip(nus, range(len(nus)))))]
        freq_ids = np.asarray(freq_ids, dtype=object)[sort]
        fluxcols = np.asarray(fluxcols, dtype=object)[sort]
        errcols = np.asarray(errcols, dtype=object)[sort]
        
        fluxes = [row[col][0] for col in fluxcols]
        errs = [row[errcol][0] for errcol in errcols]
        x = np.linspace(0.8*np.min(nus), 1.1*np.max(nus), 100)
        ys = []
        if alphas:
            if len(fluxes) <= 2:
                for a in alphas:
                    constant = fluxes[-1]/(nus[-1]**a)
                    ys.append(constant*(x**a))
            else:
                for a in alphas:
                    constant = np.median(fluxes)/(np.median(nus)**a)
                    ys.append(constant*(x**a))
        
        n_apertureplots = len(rsobjs)
        grid = gs.GridSpec(n_apertureplots, n_apertureplots, wspace=0.1, hspace=0.4)
        plt.figure(figsize=(8,10))
        
        for i, rsobj in enumerate(rsobjs):
            
            cutout, cutout_data = rsobj._make_cutouts(catalog=row, data=rsobj.data)
            cutout_data = cutout_data.squeeze()
            ax = plt.subplot(grid[i])
            ax.imshow(cutout_data, origin='lower')
            
            sidelength = np.shape(cutout_data)[1]
            beam = rsobj.beam.ellipse_to_plot(sidelength-7.5, 7.5, pixscale=rsobj.pixel_scale)
            beam.set(fill=False, ls='-', ec='darkred')
            plt.gca().add_artist(beam)
            
            apertures = [ellipse, annulus]
            ap_names = []
            pixels = []
            masks = []
            
            for aperture in apertures:
                some_pixels, a_mask = rsobj.get_pixels(aperture,
                                                       catalog=row,
                                                       data=rsobj.data,
                                                       cutouts=cutout)
                ap_names.append(aperture.__name__)
                pixels.append(some_pixels)
                masks.append(a_mask)
        
            ellipse_pix = pixels[ap_names.index('ellipse')]
            annulus_pix = pixels[ap_names.index('annulus')]
            
            snr_val = rsobj.get_snr(source=ellipse_pix, background=annulus_pix,
                                    catalog=row)[0]
            name = row['_name'][0]
            for k in range(len(masks)):
                if path is not None:
                    ax.imshow(masks[k].squeeze(), origin='lower', cmap='gray', alpha=0.4)
                else:
                    ax.imshow(masks[k].squeeze(), origin='lower', cmap='gray', alpha=0.15)
                
            plt.text(0, 0, 'SN {:.1f}'.format(snr_val), fontsize=7, color='w', 
                     ha='left', va='bottom', transform=ax.transAxes)
            plt.text(0, 1, name, fontsize=7, color='w', ha='left', va='top', 
                     transform=ax.transAxes)
            ax.set_title('Freq: {}'.format(freq_ids[i]))
            ax.set_xticks([])
            ax.set_yticks([])
            
        ax = plt.subplot(grid[n_apertureplots:])
        for i in range(len(fluxes)):
            if fluxes[i] < 3.*errs[i]:
                ax.scatter(nus[i], errs[i], marker='v', color='k', zorder=3, label=r'1 $\sigma$')
                ax.scatter(nus[i], 2.*errs[i], marker='v', color='darkred', zorder=3, label=r'2 $\sigma$')
                ax.scatter(nus[i], 3.*errs[i], marker='v', color='red', zorder=3, label=r'3 $\sigma$')
            else:
                ax.errorbar(nus[i], fluxes[i], yerr=errs[i], fmt='o', ms=2, 
                            elinewidth=0.75, color='k', zorder=3,
                            label='Aperture {}'.format(method))
            
        if ys:
            for i, y in enumerate(ys):                     
                ax.plot(x, y, '--',
                        label=r'$\alpha$ = {}'.format(alphas[i], zorder=2))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Log Frequency (GHz)')
        ax.set_ylabel('Log Flux (Jy)')
        ax.set_title('Spectral Energy Distribution for Source {}'.format(row['_name'][0]))
        handles, labels = plt.gca().get_legend_handles_labels()
        label = OrderedDict(zip(labels, handles))
        ax.legend(label.values(), label.keys())
        plt.suptitle('Name: {}'.format(row['_name'][0]))
        plt.tight_layout()
        plt.subplots_adjust(top=0.92,
                            bottom=0.075,
                            left=0.11,
                            right=0.94,
                            hspace=0.2,
                            wspace=0.2)
                            
        if path is not None:
            plt.savefig('{}SEDgrid_{}.pdf'.format(path, row['_name'][0]), dpi=300, bbox_inches='tight', overwrite=True)       
        
        return nus, fluxes, errs
        
        
