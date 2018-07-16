from astropy.table import MaskedColumn
import numpy as np
import astropy.units as u

if __package__ == '':
    __package__ = 'dendrocat'
from .utils import rms, _matcher
from .radiosource import RadioSource
from .aperture import ellipse, annulus


def match(*args, verbose=True):
    """
    Documentation needed
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
            
            data = rs_obj.data/rs_obj.ppbeam
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
                aperture.__name__+'_peak_'+rs_obj.freq_id,
                aperture.__name__+'_sum_'+rs_obj.freq_id,
                aperture.__name__+'_rms_'+rs_obj.freq_id,
                aperture.__name__+'_median_'+rs_obj.freq_id,
                aperture.__name__+'_npix_'+rs_obj.freq_id
            ]
            
            peak_data = np.zeros(len(pix_in_aperture))
            for j in range(len(pix_in_aperture)):
                peak_data[j] = np.max(pix_in_aperture[j])
            aperture_peak_col = MaskedColumn(data=peak_data,
                                             name=names[0])
            
            sum_data = np.zeros(len(pix_in_aperture))
            for j in range(len(pix_in_aperture)):
                sum_data[j] = np.sum(pix_in_aperture[j])                    
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
            
            
            
