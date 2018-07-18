from astropy.table import MaskedColumn
import numpy as np
import astropy.units as u
from copy import deepcopy
import matplotlib.pyplot as plt

if __package__ == '':
    __package__ = 'dendrocat'
from .utils import rms, _matcher, specindex
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


    def add_objects(self, *args):
        obj_prefix = 'radiosource_'
        for obj in args:
            if isinstance(obj, MasterCatalog):
                for key in obj.__dict__.keys():
                    if key.split('_')[0]+'_' == obj_prefix:
                        self.__dict__[key] = obj[key]
            else:        
                self.__dict__[obj_prefix+obj.freq_id] = obj
    
    
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
                sum_data[j] = np.sum(pix_in_aperture[j])/rs_obj.ppbeam                   
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
                    npix_data[j] = None
                else:
                    npix_data[j] = len(pix_in_aperture[j])
                aperture_npix_col = MaskedColumn(data=npix_data,
                                                 name=names[4])
            
            try:
                self.catalog.remove_columns(names)
            except KeyError:
                self.catalog.add_columns([
                    aperture_peak_col,
                    aperture_sum_col,
                    aperture_rms_col,
                    aperture_median_col,
                    aperture_npix_col
                ])  
            
            # Mask NaN values        
            for col in self.catalog.colnames:
                isnan = np.argwhere(np.isnan(list(self.catalog[col])))
                self.catalog.mask[col][isnan] = True
            
            
    def ffplot(self, rsobj1, rsobj2, apertures=[], specs=[], peak=False,
               label=False, log=True, outfile=None):
           
        if rsobj1.nu > rsobj2.nu:
            rsobj1, rsobj2 = rsobj2, rsobj1
        
        catalog = deepcopy(self.catalog)
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
        apertures = list(set(apertures)|{ellipse, annulus})
        if specs is None:
            specs = [1, 2, 3]
        
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
        npix1 = []
        npix2 = []
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
                                     
            npix1.append(catalog[rsobj1.freq_id+'_'+aperture.__name__+'_npix'])
            npix2.append(catalog[rsobj2.freq_id+'_'+aperture.__name__+'_npix'])
        
        err1.append(catalog[rsobj1.freq_id+'_annulus_rms'])
        err2.append(catalog[rsobj2.freq_id+'_annulus_rms'])
        
        marker_labels = catalog['_idx']
        
        xflux = np.linspace(np.min(flux1), np.max(flux1), 10)
        yfluxes = []
        
        for spec in specs:
            yfluxes.append(specindex(rsobj1.nu, rsobj2.nu, xflux, spec))
        
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
                        label='Spectral Index = {}'.format(specs[j]))
            
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
            
