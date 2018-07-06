from astropy.table import Table, MaskedColumn
import astropy.units as u
from astropy import coordinates
from astropy.nddata.utils import Cutout2D
import numpy as np

class SourceCatalog(Table):

    def __init__(self, catalog=None, imageobj=None, masked=None, names=None, 
                 dtype=None, meta=None, copy=True, rows=None, 
                 copy_indices=True, **kwargs):
                 
        Table.__init__(self, data=catalog, masked=masked, names=names,
                       dtype=dtype, copy=copy, rows=rows,
                       copy_indices=copy_indices, **kwargs)
                      
        if imageobj: 
            self.__dict__.update(imageobj.__dict__)
        self.catalog = catalog
        self.names = names
        self.meta = meta
        self.copy = copy
        self.rows = rows
        self.copy_indices = copy_indices


    def save_ds9_regions(self, outfile):
        with open(outfile, 'w') as fh:
            fh.write("icrs\n")
            for row in cat:
                fh.write("ellipse({x_cen}, {y_cen}, {major_fwhm}, " \
                         "{minor_fwhm}, {position_angle}) # text={{{_idx}}}\n"
                         .format(**dict(zip(row.colnames, row))))
    
    
    def make_cutouts(self, sidelength, imageobj=None, catalog=None, save=True):
        """
        Make a datacube of cutout regions around all source centers in the 
        catalog.
        
        Parameters
        ----------
        sidelength : float
            Side length of the square (in pixels) to cut out of the image for 
            each source.
        imageobj : rsprocess.image.Image object, optional
            Image object from which to make the cutouts. If unspecified, image 
            information from instance attributes will be used.
        catalog : astropy.table.Table object, optional
            Source catalog to use for cutout coordinates. If unspecified, 
            catalog stored in instance attributes will be used.
        save : bool, optional
            If enabled, the cutouts and cutout data will both be saved as 
            instance attributes. Default is True.
        Returns
        ----------
        List of astropy.nddata.utils.Cutout2D objects, list of cutout data
            
        """
    
        if imageobj:
            beam = imageobj.beam
            wcs = imageobj.wcs
            pixel_scale = imageobj.pixel_scale
            data = imageobj.data
        else:
            beam = self.beam
            wcs = self.wcs
            pixel_scale = self.pixel_scale
            data = self.data
            
        if not catalog:
            catalog = self.catalog
        
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
        
        if save:
            self.cutouts = cutouts
            self.cutout_data = cutout_data
        
        return cutouts, cutout_data
    
    
    def snr(self, cutouts=None):
        
        """
        Return the SNR of all sources in the catalog. Calculated using the peak 
        flux in the ellipse region over the background annulus RMS.
        """
        #self.snr = snr
        
        # return snr
    
    
    def reject(self, imgobj=None, catalog=None, threshold=6.):
        """
        Reject noisy detections.
        
        Parameters
        ----------
        imgobj : rsprocess.image.Image object, optional
            The radio image from which the source catalog was extracted.
        catalog : astropy.table.Table object, optional
            The source catalog 
        threshold : float, optional
            The signal-to-noise threshold below which sources are rejected
        """
        
        
    
    
    
if __name__ == '__main__':
    from astropy.io import fits
    from image import Image
    
    filename = '/lustre/aoc/students/bmcclell/w51/W51e2_cont_briggsSC_tclean.image.fits.gz'
    f = fits.open(filename)
    i = Image(f)
    i.to_dendrogram()
    c = i.to_cat()
