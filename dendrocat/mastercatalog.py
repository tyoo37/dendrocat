from astropy.table import MaskedColumn
import numpy as np
import astropy.units as u

if __package__ == '':
    __package__ = 'dendrocat'
from .utils import rms
from .radiosource import RadioSource


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
                
        aperture_npix_col = MaskedColumn(length=len(catalog),
                                         name=aperture.__name__,
                                         mask=True)
        
        for i, rs_obj in enumerate(rs_objects):
        
            size = 2.2*(np.max(catalog['major_fwhm'])*u.deg 
                        + rs_obj.annulus_padding 
                        + rs_obj.annulus_width)
                        
            cutouts, cutout_data = rs_obj._make_cutouts(size, 
                                                        catalog=catalog, 
                                                        save=False)
        
            try:
                pix_in_aperture = rs_obj.__dict__['pix_'+aperture.__name__]
            except KeyError:                              
                pix_in_aperture = rs_obj.get_pixels(
                                                    aperture,
                                                    cutouts=cutouts,
                                                    cutout_data=cutout_data,
                                                    **kwargs
                                                    )
            
            names = [
                aperture.__name__+'_peak_'+rs_obj.freq_id,
                aperture.__name__+'_sum_'+rs_obj.freq_id,
                aperture.__name__+'_rms_'+rs_obj.freq_id,
                aperture.__name__+'_median_'+rs_obj.freq_id
            ]
            
            peak_data = np.amax(np.amax(pix_in_aperture, axis=1), axis=1)
            aperture_peak_col = MaskedColumn(data=peak_data,
                                             name=names[0],
            
                                             mask=True)
            sum_data = np.sum(np.sum(pix_in_aperture, axis=1), axis=1)                         
            aperture_sum_col = MaskedColumn(data=sum_data,
                                            name=names[1],
                                            mask=True)                  
                                            
            rms_data = np.zeros(len(pix_in_aperture))
            for i in range(len(pix_in_aperture)):
                rms_data[i] = rms(pix_in_aperture[i])                           
            aperture_rms_col = MaskedColumn(data=rms_data,
                                            name=names[2],
                                            mask=True)
            
            median_data = np.median(np.median(pix_in_aperture, axis=1), axis-1)                   
            aperture_median_col = MaskedColumn(data=median_data,
                                               name=names[3],
                                               mask=True)
            
            self.catalog.add_columns([
                aperture_peak_col,
                aperture_sum_col,
                aperture_rms_col,
                aperture_median_col
            ])
            
            
            
